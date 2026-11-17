//! Binds the names of a parsed query format to fields and to the tags a header
//! declares (D172).

use gloci_io::{UsageError, VcfHeader};

use super::parse::{Item, Spec};
use super::{Field, Piece, Query, Tag, TagType};

struct Declared {
    info: Vec<(Vec<u8>, TagType)>,
    format: Vec<(Vec<u8>, TagType)>,
}

fn tag_type(ty: Option<&[u8]>) -> TagType {
    match ty {
        Some(b"Integer") => TagType::Integer,
        Some(b"Float") => TagType::Float,
        Some(b"Flag") => TagType::Flag,
        _ => TagType::String,
    }
}

fn declared(header: &VcfHeader, kind: &'static [u8]) -> Vec<(Vec<u8>, TagType)> {
    crate::extract::declarations(&header.raw, kind)
        .map(|d| (d.id.to_vec(), tag_type(d.ty)))
        .collect()
}

pub(super) fn resolve(
    src: &str,
    items: &[Item],
    header: &VcfHeader,
    print_header: bool,
) -> Result<Query, UsageError> {
    let decl = Declared {
        info: declared(header, b"INFO"),
        format: declared(header, b"FORMAT"),
    };
    let cx = Cx {
        src,
        decl: &decl,
        print_header,
    };
    let mut pieces = cx.items(items, false)?;
    if !has_newline(&pieces) {
        pieces.push(Piece::Text(b"\n".to_vec()));
    }
    Ok(Query {
        pieces,
        warned: Default::default(),
    })
}

fn has_newline(pieces: &[Piece]) -> bool {
    pieces.iter().any(|p| match p {
        Piece::Text(t) => t.contains(&b'\n'),
        Piece::Block(inner) => has_newline(inner),
        Piece::Field(_) => false,
    })
}

struct Cx<'a> {
    src: &'a str,
    decl: &'a Declared,
    print_header: bool,
}

impl Cx<'_> {
    fn fail(&self, msg: &str, col: usize) -> UsageError {
        UsageError(format!(
            "invalid query format \"{}\": {msg} (column {col})",
            self.src
        ))
    }

    fn items(&self, items: &[Item], in_block: bool) -> Result<Vec<Piece>, UsageError> {
        let mut out = Vec::new();
        for item in items {
            out.push(match item {
                Item::Text(t) => Piece::Text(t.clone()),
                Item::Block(inner) => Piece::Block(self.items(inner, true)?),
                Item::Field(spec) => Piece::Field(self.field(spec, in_block)?),
            });
        }
        Ok(out)
    }

    fn field(&self, spec: &Spec, in_block: bool) -> Result<Field, UsageError> {
        let name = spec.name.as_str();
        let plain = |f: Field| -> Result<Field, UsageError> {
            match spec.sub {
                None => Ok(f),
                Some(_) => Err(self.fail(&format!("{name} takes no subscript"), spec.col)),
            }
        };
        let in_block_only = |f: Field| -> Result<Field, UsageError> {
            if !in_block {
                return Err(self.fail(&format!("%{name} is only valid inside [ ]"), spec.col));
            }
            plain(f)
        };
        match name {
            "CHROM" => plain(Field::Chrom),
            "POS" => plain(Field::Pos),
            "POS0" => plain(Field::Pos0),
            "END" => plain(Field::End),
            "END0" => plain(Field::End0),
            "ID" => plain(Field::Id),
            "REF" => plain(Field::Ref),
            "ALT" => Ok(Field::Alt(spec.sub)),
            "FIRST_ALT" => plain(Field::FirstAlt),
            "QUAL" => plain(Field::Qual),
            "FILTER" => plain(Field::Filter),
            "TYPE" => plain(Field::Type),
            "INFO" => plain(Field::Info),
            "LINE" => {
                if in_block {
                    return Err(self.fail("%LINE is only valid outside [ ]", spec.col));
                }
                if self.print_header {
                    return Err(self.fail("-H cannot be used with %LINE", spec.col));
                }
                plain(Field::Line)
            }
            "SAMPLE" => in_block_only(Field::Sample),
            "GT" => in_block_only(Field::Gt),
            "TGT" => in_block_only(Field::Tgt),
            _ => self.tag(spec, in_block),
        }
    }

    fn tag(&self, spec: &Spec, in_block: bool) -> Result<Field, UsageError> {
        let name = spec.name.as_str();
        if let Some((_, tag)) = name
            .strip_prefix("FMT/")
            .or_else(|| name.strip_prefix("FORMAT/"))
            .map(|t| (name, t))
        {
            return Err(self.fail(
                &format!("%{name} is not supported; write [%{tag}] instead"),
                spec.col,
            ));
        }
        let make = |tag: &str, ty| Tag {
            name: tag.as_bytes().to_vec(),
            ty,
            sub: spec.sub,
        };
        let find = |list: &[(Vec<u8>, TagType)], tag: &str| {
            list.iter().find(|(n, _)| n == tag.as_bytes()).map(|p| p.1)
        };
        if let Some(tag) = name.strip_prefix("INFO/") {
            return match find(&self.decl.info, tag) {
                Some(ty) => Ok(Field::InfoTag(make(tag, ty))),
                None => Err(self.fail(
                    &format!("%{name}: INFO tag {tag} is not declared in the header"),
                    spec.col,
                )),
            };
        }
        if in_block && let Some(ty) = find(&self.decl.format, name) {
            return Ok(Field::FormatTag(make(name, ty)));
        }
        if let Some(ty) = find(&self.decl.info, name) {
            return Ok(Field::InfoTag(make(name, ty)));
        }
        if find(&self.decl.format, name).is_some() {
            return Err(self.fail(
                &format!("%{name} is a FORMAT tag, which is only valid inside [ ]"),
                spec.col,
            ));
        }
        Err(self.fail(
            &format!("%{name}: {name} is not declared in the header"),
            spec.col,
        ))
    }
}
