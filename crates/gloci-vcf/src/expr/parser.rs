//! Recursive-descent parser: tokens to a typed [`Node`] tree.

use super::lexer::{Op, Spanned, SyntaxError, Tok, err, lex};
use super::{Expr, Field, InfoTy, Node, Rhs, TypeKw};
use gloci_io::{UsageError, VcfHeader};
use std::collections::BTreeMap;

/// How deeply parentheses, `!` and `N_PASS(...)` may nest. The parser recurses
/// once per level, so an unbounded depth would overflow the stack.
const MAX_DEPTH: usize = 128;

const FUNCS: &str = "N_PASS";
const TYPE_KEYWORDS: &str = "snp, mnp, indel, other, ref";
const GT_KEYWORDS: [(&str, u8); 6] = [
    ("het", HET),
    ("hom", HOM_REF | HOM_ALT),
    ("ref", HOM_REF | HAP_REF),
    ("alt", HOM_ALT | HET | HAP_ALT),
    ("mis", MIS),
    ("hap", HAP_REF | HAP_ALT),
];

// Genotype classes (bit set); see `eval::gt_class`.
pub(crate) const HOM_REF: u8 = 1;
pub(crate) const HOM_ALT: u8 = 2;
pub(crate) const HET: u8 = 4;
pub(crate) const HAP_REF: u8 = 8;
pub(crate) const HAP_ALT: u8 = 16;
pub(crate) const MIS: u8 = 32;

/// What the header defines.
struct Schema {
    info: BTreeMap<Vec<u8>, InfoTy>,
    format: BTreeMap<Vec<u8>, ()>,
    filters: Vec<Vec<u8>>,
    n_samples: usize,
}

/// The value of `key` in a `##X=<ID=...,Key=Value,...>` line.
fn def_attr(line: &[u8], key: &[u8]) -> Option<Vec<u8>> {
    let start = line.iter().position(|&b| b == b'<')? + 1;
    let body = line[start..].strip_suffix(b">").unwrap_or(&line[start..]);
    let (mut i, mut in_q) = (0, false);
    let mut field_start = 0;
    let mut fields = Vec::new();
    while i <= body.len() {
        if i == body.len() || (body[i] == b',' && !in_q) {
            fields.push(&body[field_start..i]);
            field_start = i + 1;
        } else if body[i] == b'"' {
            in_q = !in_q;
        }
        i += 1;
    }
    fields.into_iter().find_map(|f| {
        let eq = f.iter().position(|&b| b == b'=')?;
        (&f[..eq] == key).then(|| f[eq + 1..].to_vec())
    })
}

impl Schema {
    fn new(header: &VcfHeader) -> Schema {
        let mut s = Schema {
            info: BTreeMap::new(),
            format: BTreeMap::new(),
            filters: Vec::new(),
            n_samples: header.samples.len(),
        };
        for ((kind, id), line) in header.definitions() {
            match kind {
                "INFO" => {
                    let ty = match def_attr(&line, b"Type").as_deref() {
                        Some(b"Integer" | b"Float") => InfoTy::Num,
                        Some(b"Flag") => InfoTy::Flag,
                        _ => InfoTy::Str,
                    };
                    s.info.insert(id, ty);
                }
                "FORMAT" => {
                    s.format.insert(id, ());
                }
                "FILTER" => s.filters.push(id),
                _ => {}
            }
        }
        s
    }
}

pub(crate) fn parse(src: &str, header: &VcfHeader) -> Result<Expr, UsageError> {
    let shown: String = if src.chars().count() > 200 {
        format!("{}...", src.chars().take(200).collect::<String>())
    } else {
        src.to_string()
    };
    let fail = |e: SyntaxError| {
        UsageError(format!(
            "invalid filter expression {shown:?}: {} (column {})",
            e.msg, e.col
        ))
    };
    let toks = lex(src).map_err(fail)?;
    let mut p = Parser {
        toks,
        pos: 0,
        schema: Schema::new(header),
        in_npass: false,
        depth: 0,
        flags: Flags::default(),
    };
    let root = p.parse_or().map_err(fail)?;
    let t = p.peek().clone();
    if t.tok != Tok::End {
        return err(
            t.col,
            format!(
                "unexpected {}; expected `&&`, `||` or the end of the expression",
                t.tok.describe()
            ),
        )
        .map_err(fail);
    }
    Ok(Expr {
        root,
        needs_stats: p.flags.stats,
        needs_genotypes: p.flags.genotypes || p.flags.stats,
        needs_info: p.flags.info,
    })
}

#[derive(Default)]
struct Flags {
    stats: bool,
    genotypes: bool,
    info: bool,
}

struct Parser {
    toks: Vec<Spanned>,
    pos: usize,
    schema: Schema,
    in_npass: bool,
    depth: usize,
    flags: Flags,
}

/// A parsed term before it is placed in a comparison.
enum Term {
    Num(f64),
    Str(String),
    Miss,
    Field(Field),
    Bool(Node),
}

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Num,
    Str,
}

fn kind_of(f: &Field) -> Kind {
    match f {
        Field::Chrom
        | Field::Id
        | Field::Ref
        | Field::Alt
        | Field::Filter
        | Field::Type
        | Field::Gt(_)
        | Field::Info {
            ty: InfoTy::Str, ..
        } => Kind::Str,
        _ => Kind::Num,
    }
}

fn is_sample_wise(f: &Field) -> bool {
    matches!(f, Field::FmtDp(_) | Field::FmtGq(_) | Field::Gt(_))
}

fn node_is_sample_wise(n: &Node) -> bool {
    match n {
        Node::Or(v) | Node::And(v) => v.iter().any(node_is_sample_wise),
        Node::Not(n) => node_is_sample_wise(n),
        Node::Cmp { field, .. } | Node::Truthy(field) => is_sample_wise(field),
    }
}

impl Parser {
    fn peek(&self) -> &Spanned {
        &self.toks[self.pos]
    }

    fn bump(&mut self) -> Spanned {
        let t = self.toks[self.pos].clone();
        if self.pos + 1 < self.toks.len() {
            self.pos += 1;
        }
        t
    }

    /// Count one more level of nesting; error when it is too deep.
    fn descend(&mut self) -> Result<(), SyntaxError> {
        if self.depth >= MAX_DEPTH {
            return err(
                self.peek().col,
                format!(
                    "expression nested too deeply (more than {MAX_DEPTH} levels of parentheses or `!`); \
simplify the expression"
                ),
            );
        }
        self.depth += 1;
        Ok(())
    }

    fn parse_or(&mut self) -> Result<Node, SyntaxError> {
        self.descend()?;
        let r = self.parse_or_list();
        self.depth -= 1;
        r
    }

    fn parse_or_list(&mut self) -> Result<Node, SyntaxError> {
        let first = self.parse_and()?;
        if self.peek().tok != Tok::Or {
            return Ok(first);
        }
        let mut v = vec![first];
        while self.peek().tok == Tok::Or {
            self.bump();
            v.push(self.parse_and()?);
        }
        Ok(Node::Or(v))
    }

    fn parse_and(&mut self) -> Result<Node, SyntaxError> {
        let first = self.parse_unary()?;
        if self.peek().tok != Tok::And {
            return Ok(first);
        }
        let mut v = vec![first];
        while self.peek().tok == Tok::And {
            self.bump();
            v.push(self.parse_unary()?);
        }
        Ok(Node::And(v))
    }

    fn parse_unary(&mut self) -> Result<Node, SyntaxError> {
        if self.peek().tok == Tok::Not {
            self.descend()?;
            self.bump();
            let inner = self.parse_unary();
            self.depth -= 1;
            return Ok(Node::Not(Box::new(inner?)));
        }
        self.parse_cmp()
    }

    fn parse_cmp(&mut self) -> Result<Node, SyntaxError> {
        let lcol = self.peek().col;
        let l = self.parse_term()?;
        let Tok::Cmp(op) = self.peek().tok else {
            // `QUAL+1>10` lexes as QUAL followed by the number +1
            if matches!(self.peek().tok, Tok::Num(_)) && matches!(l, Term::Field(_)) {
                return err(
                    self.peek().col,
                    "arithmetic is not supported; compare a field with a number, \
                     a string or `.`",
                );
            }
            return self.bare(l, lcol);
        };
        let opcol = self.bump().col;
        let rcol = self.peek().col;
        let r = self.parse_term()?;
        let (field, lit, op, fcol, lcol) = match (l, r) {
            (Term::Bool(_), _) => {
                return err(
                    lcol,
                    "a parenthesised condition can't be compared; compare fields or values instead",
                );
            }
            (_, Term::Bool(_)) => {
                return err(
                    rcol,
                    "a parenthesised condition can't be compared; compare fields or values instead",
                );
            }
            (Term::Field(f), lit @ (Term::Num(_) | Term::Str(_) | Term::Miss)) => {
                (f, lit, op, lcol, rcol)
            }
            (lit @ (Term::Num(_) | Term::Str(_) | Term::Miss), Term::Field(f)) => {
                if matches!(op, Op::Match | Op::NoMatch) {
                    return err(
                        opcol,
                        "the regular expression must be on the right of `~` / `!~`; write FIELD ~ \"pattern\"",
                    );
                }
                (f, lit, op.flipped(), rcol, lcol)
            }
            (Term::Field(_), Term::Field(_)) => {
                return err(
                    opcol,
                    "comparing two fields is not supported; compare a field with a number or a string",
                );
            }
            _ => {
                return err(
                    opcol,
                    "a comparison needs a field on one side; use QUAL, POS, INFO/KEY, FMT/KEY, ...",
                );
            }
        };
        let rhs = self.make_rhs(&field, op, lit, fcol, lcol, opcol)?;
        Ok(Node::Cmp { field, op, rhs })
    }

    /// A term used without an operator.
    fn bare(&mut self, t: Term, col: usize) -> Result<Node, SyntaxError> {
        match t {
            Term::Bool(n) => Ok(n),
            Term::Field(
                f @ Field::Info {
                    ty: InfoTy::Flag, ..
                },
            ) => Ok(Node::Truthy(f)),
            Term::Field(_) => err(
                col,
                "this field is used without a comparison; compare it, e.g. QUAL>30 \
                 (only INFO flags can stand alone)",
            ),
            _ => err(col, "a value is used without a comparison"),
        }
    }

    fn make_rhs(
        &self,
        field: &Field,
        op: Op,
        lit: Term,
        fcol: usize,
        lcol: usize,
        opcol: usize,
    ) -> Result<Rhs, SyntaxError> {
        let fk = kind_of(field);
        let is_re = matches!(op, Op::Match | Op::NoMatch);
        match (lit, field) {
            (Term::Str(s), Field::Type) => {
                if is_re || !matches!(op, Op::Eq | Op::Ne) {
                    return err(opcol, "TYPE supports only = and !=");
                }
                let k = match s.as_str() {
                    "snp" => TypeKw::Snp,
                    "mnp" => TypeKw::Mnp,
                    "indel" => TypeKw::Indel,
                    "other" => TypeKw::Other,
                    "ref" => TypeKw::Ref,
                    _ => {
                        return err(
                            lcol,
                            format!("{s:?} is not a TYPE; use one of {TYPE_KEYWORDS}"),
                        );
                    }
                };
                Ok(Rhs::TypeKw(k))
            }
            (Term::Miss, Field::Type) => err(
                lcol,
                format!("TYPE is compared with one of {TYPE_KEYWORDS}"),
            ),
            (Term::Str(s), Field::Gt(_)) if !is_re => {
                if !matches!(op, Op::Eq | Op::Ne) {
                    return err(opcol, "GT supports only =, !=, ~ and !~");
                }
                match GT_KEYWORDS.iter().find(|(k, _)| *k == s) {
                    Some((_, m)) => Ok(Rhs::GtKw(*m)),
                    None => Ok(Rhs::Str(s.into_bytes().into())),
                }
            }
            (Term::Str(s), _) if is_re => {
                if fk != Kind::Str {
                    return err(
                        fcol,
                        "~ and !~ match text; this field is numeric (use a string field such as ID, REF, ALT, FILTER, GT)",
                    );
                }
                regex::bytes::Regex::new(&s)
                    .map(Rhs::Re)
                    .map_err(|e| SyntaxError {
                        col: lcol,
                        msg: format!(
                            "{s:?} is not a valid regular expression ({}); fix the pattern",
                            e.to_string().lines().last().unwrap_or("syntax error").trim()
                        ),
                    })
            }
            (_, _) if is_re => err(lcol, "~ and !~ need a quoted pattern, e.g. ID~\"^rs\""),
            (Term::Str(s), _) => {
                if fk != Kind::Str {
                    return err(
                        lcol,
                        format!("{s:?} is text but this field is numeric; write the number without quotes"),
                    );
                }
                if !matches!(op, Op::Eq | Op::Ne) {
                    return err(opcol, "text fields support only =, !=, ~ and !~");
                }
                if matches!(field, Field::Filter)
                    && !matches!(s.as_str(), "PASS" | ".")
                    && !s.split(';').all(|t| {
                        matches!(t, "PASS" | ".")
                            || self.schema.filters.iter().any(|f| f == t.as_bytes())
                    })
                {
                    return err(
                        lcol,
                        format!("filter {s:?} is not defined in the header; add a ##FILTER line or fix the name"),
                    );
                }
                Ok(Rhs::Str(s.into_bytes().into()))
            }
            (Term::Num(n), _) => {
                if fk != Kind::Num {
                    return err(
                        lcol,
                        "this field is text but is compared with a number; quote the value, e.g. \"1\"",
                    );
                }
                Ok(Rhs::Num(n))
            }
            (Term::Miss, _) => {
                if !matches!(op, Op::Eq | Op::Ne) {
                    return err(opcol, "`.` (missing) can be compared only with = and !=");
                }
                Ok(if fk == Kind::Str {
                    Rhs::Str(b".".to_vec().into())
                } else {
                    Rhs::Miss
                })
            }
            _ => unreachable!("parse_cmp passes only literals"),
        }
        .and_then(|rhs| {
            if matches!(op, Op::Lt | Op::Le | Op::Gt | Op::Ge) && !matches!(rhs, Rhs::Num(_)) {
                return err(opcol, "< <= > >= compare numbers only");
            }
            Ok(rhs)
        })
    }

    fn parse_term(&mut self) -> Result<Term, SyntaxError> {
        let t = self.bump();
        match t.tok {
            Tok::Num(n) => Ok(Term::Num(n)),
            Tok::Str(s) => Ok(Term::Str(s)),
            Tok::Dot => Ok(Term::Miss),
            Tok::LParen => {
                let n = self.parse_or()?;
                let c = self.bump();
                if c.tok != Tok::RParen {
                    return err(
                        c.col,
                        format!(
                            "expected `)` but found {}; close the parenthesis opened at column {}",
                            c.tok.describe(),
                            t.col
                        ),
                    );
                }
                Ok(Term::Bool(n))
            }
            Tok::Ident(name) => {
                if self.peek().tok == Tok::LParen {
                    return self.parse_func(&name, t.col);
                }
                self.parse_field(&name, t.col).map(Term::Field)
            }
            Tok::End => err(
                t.col,
                "the expression ends here but a field or value is expected",
            ),
            other => err(
                t.col,
                format!(
                    "unexpected {}; expected a field, a number, a string or `(`",
                    other.describe()
                ),
            ),
        }
    }

    fn parse_func(&mut self, name: &str, col: usize) -> Result<Term, SyntaxError> {
        if name != "N_PASS" {
            return err(
                col,
                format!("function {name} is not supported; supported functions: {FUNCS}"),
            );
        }
        if self.in_npass {
            return err(col, "N_PASS can't be nested");
        }
        self.bump(); // (
        self.in_npass = true;
        let inner_col = self.peek().col;
        let inner = self.parse_or();
        self.in_npass = false;
        let inner = inner?;
        let c = self.bump();
        if c.tok != Tok::RParen {
            return err(
                c.col,
                format!(
                    "expected `)` to close N_PASS( but found {}",
                    c.tok.describe()
                ),
            );
        }
        if !node_is_sample_wise(&inner) {
            return err(
                inner_col,
                "N_PASS needs a per-sample condition such as FMT/GQ>30 or GT=\"het\"",
            );
        }
        self.flags.genotypes = true;
        Ok(Term::Field(Field::NPass(Box::new(inner))))
    }

    /// An optional `[INT]` after a field.
    fn parse_index(&mut self) -> Result<Option<usize>, SyntaxError> {
        if self.peek().tok != Tok::LBrack {
            return Ok(None);
        }
        self.bump();
        let t = self.bump();
        let idx = match t.tok {
            Tok::Num(n) if n >= 0.0 && n.fract() == 0.0 && n < 1e9 => n as usize,
            _ => {
                return err(
                    t.col,
                    "the index must be a whole number from 0, e.g. INFO/AF[0]",
                );
            }
        };
        let c = self.bump();
        if c.tok != Tok::RBrack {
            return err(c.col, "expected `]` after the index");
        }
        Ok(Some(idx))
    }

    /// An optional `[N]` sample subscript after a FORMAT field: N is numbered
    /// over the whole header and must exist.
    fn parse_sample_index(&mut self) -> Result<Option<usize>, SyntaxError> {
        let col = self.toks.get(self.pos + 1).map_or(1, |t| t.col);
        let Some(i) = self.parse_index()? else {
            return Ok(None);
        };
        let n = self.schema.n_samples;
        if i >= n {
            return err(
                col,
                format!(
                    "the sample index {i} is too large: the header has {n} sample{}; \
                     use an index from 0 to {}",
                    if n == 1 { "" } else { "s" },
                    n.saturating_sub(1)
                ),
            );
        }
        Ok(Some(i))
    }

    fn parse_field(&mut self, name: &str, col: usize) -> Result<Field, SyntaxError> {
        let no_index = |p: &mut Parser, f: Field| -> Result<Field, SyntaxError> {
            let c = p.peek().col;
            if p.peek().tok == Tok::LBrack {
                return err(c, format!("{name} is a single value and can't be indexed"));
            }
            Ok(f)
        };
        match name {
            "QUAL" => return no_index(self, Field::Qual),
            "POS" => return no_index(self, Field::Pos),
            "CHROM" => return no_index(self, Field::Chrom),
            "ID" => return no_index(self, Field::Id),
            "REF" => return no_index(self, Field::Ref),
            "ALT" => return no_index(self, Field::Alt),
            "FILTER" => return no_index(self, Field::Filter),
            "TYPE" => return no_index(self, Field::Type),
            "N_ALT" => return no_index(self, Field::NAlt),
            "N_SAMPLES" => {
                let f = Field::NSamples(self.schema.n_samples);
                return no_index(self, f);
            }
            "GT" => {
                self.flags.genotypes = true;
                return Ok(Field::Gt(self.parse_sample_index()?));
            }
            // bcftools reads a header-defined INFO/AC, AN or MAF here and
            // computes them only when the header lacks the key; F_MISSING
            // and N_MISSING are always computed.
            "MAF" | "AN" | "AC" if self.schema.info.contains_key(name.as_bytes()) => {}
            "MAF" | "AN" | "F_MISSING" | "N_MISSING" => {
                self.flags.stats = true;
                let f = match name {
                    "MAF" => Field::Maf,
                    "AN" => Field::An,
                    "F_MISSING" => Field::FMissing,
                    _ => Field::NMissing,
                };
                return no_index(self, f);
            }
            "AC" => {
                self.flags.stats = true;
                return Ok(Field::Ac(self.parse_index()?));
            }
            _ => {}
        }
        let (scope, key) = match name.split_once('/') {
            Some(("INFO", k)) => (Some(true), k),
            Some(("FMT" | "FORMAT", k)) => (Some(false), k),
            Some((p, _)) => {
                return err(
                    col,
                    format!("unknown prefix {p:?}; use INFO/, FMT/ or FORMAT/"),
                );
            }
            None => (None, name),
        };
        let kb = key.as_bytes();
        let is_info = self.schema.info.contains_key(kb);
        let is_fmt = self.schema.format.contains_key(kb) || key == "GT";
        let use_info = match scope {
            Some(true) => is_info,
            Some(false) => false,
            None if is_info && is_fmt && key != "GT" => {
                return err(
                    col,
                    format!(
                        "{key} is ambiguous: it is defined as both INFO and FORMAT in the header; \
                         write INFO/{key} or FMT/{key}"
                    ),
                );
            }
            None => is_info,
        };
        if use_info {
            let ty = self.schema.info[kb];
            self.flags.info = true;
            let idx = self.parse_index()?;
            return Ok(Field::Info {
                key: kb.into(),
                ty,
                idx,
            });
        }
        if scope != Some(true) && is_fmt {
            self.flags.genotypes = true;
            return match key {
                "GT" => Ok(Field::Gt(self.parse_sample_index()?)),
                "DP" => Ok(Field::FmtDp(self.parse_sample_index()?)),
                "GQ" => Ok(Field::FmtGq(self.parse_sample_index()?)),
                _ => err(
                    col,
                    format!("FMT/{key} is not supported; supported FORMAT fields: GT, DP, GQ"),
                ),
            };
        }
        match scope {
            Some(true) => err(
                col,
                format!(
                    "INFO/{key} is not defined in the header; add an ##INFO line or fix the name"
                ),
            ),
            Some(false) => err(
                col,
                format!(
                    "FORMAT/{key} is not defined in the header; add a ##FORMAT line or fix the name"
                ),
            ),
            None => err(
                col,
                format!(
                    "the tag {key:?} is not defined in the header (not an INFO or FORMAT key) \
                     and is not a built-in field; fix the name"
                ),
            ),
        }
    }
}
