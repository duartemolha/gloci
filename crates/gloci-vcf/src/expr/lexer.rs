//! Tokens of the filter expression language.

/// A comparison operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Op {
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
    Match,
    NoMatch,
}

impl Op {
    /// The operator with its operands swapped (`30 < QUAL` is `QUAL > 30`).
    pub(crate) fn flipped(self) -> Op {
        match self {
            Op::Lt => Op::Gt,
            Op::Le => Op::Ge,
            Op::Gt => Op::Lt,
            Op::Ge => Op::Le,
            o => o,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Tok {
    Num(f64),
    Str(String),
    Ident(String),
    /// A lone `.`: the missing value.
    Dot,
    LParen,
    RParen,
    LBrack,
    RBrack,
    Not,
    And,
    Or,
    Cmp(Op),
    End,
}

impl Tok {
    pub(crate) fn describe(&self) -> String {
        match self {
            Tok::Num(n) => format!("the number {n}"),
            Tok::Str(s) => format!("the string {s:?}"),
            Tok::Ident(s) => format!("`{s}`"),
            Tok::Dot => "`.`".into(),
            Tok::LParen => "`(`".into(),
            Tok::RParen => "`)`".into(),
            Tok::LBrack => "`[`".into(),
            Tok::RBrack => "`]`".into(),
            Tok::Not => "`!`".into(),
            Tok::And => "`&&`".into(),
            Tok::Or => "`||`".into(),
            Tok::Cmp(_) => "a comparison operator".into(),
            Tok::End => "the end of the expression".into(),
        }
    }
}

/// A token and its 1-based column (in characters).
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Spanned {
    pub tok: Tok,
    pub col: usize,
}

/// A lexing or parsing failure: what is wrong, where, and how to fix it.
#[derive(Debug)]
pub(crate) struct SyntaxError {
    pub col: usize,
    pub msg: String,
}

pub(crate) fn err<T>(col: usize, msg: impl Into<String>) -> Result<T, SyntaxError> {
    Err(SyntaxError {
        col,
        msg: msg.into(),
    })
}

fn is_ident_start(c: char) -> bool {
    c.is_ascii_alphabetic() || c == '_'
}

fn is_ident_char(c: char) -> bool {
    c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '/')
}

pub(crate) fn lex(src: &str) -> Result<Vec<Spanned>, SyntaxError> {
    let chars: Vec<char> = src.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    while i < chars.len() {
        let c = chars[i];
        let col = i + 1;
        let next = chars.get(i + 1).copied();
        let mut push = |tok: Tok, len: usize, i: &mut usize| {
            out.push(Spanned { tok, col });
            *i += len;
        };
        match c {
            _ if c.is_whitespace() => i += 1,
            '(' => push(Tok::LParen, 1, &mut i),
            ')' => push(Tok::RParen, 1, &mut i),
            '[' => push(Tok::LBrack, 1, &mut i),
            ']' => push(Tok::RBrack, 1, &mut i),
            '&' if next == Some('&') => push(Tok::And, 2, &mut i),
            '|' if next == Some('|') => push(Tok::Or, 2, &mut i),
            '&' | '|' => {
                return err(
                    col,
                    format!(
                        "`{c}` is not supported; use && / || (the sample-wise `&` and `|` are not available)"
                    ),
                );
            }
            '!' => match next {
                Some('=') => push(Tok::Cmp(Op::Ne), 2, &mut i),
                Some('~') => push(Tok::Cmp(Op::NoMatch), 2, &mut i),
                _ => push(Tok::Not, 1, &mut i),
            },
            '=' if next == Some('=') => push(Tok::Cmp(Op::Eq), 2, &mut i),
            '=' => push(Tok::Cmp(Op::Eq), 1, &mut i),
            '<' if next == Some('=') => push(Tok::Cmp(Op::Le), 2, &mut i),
            '<' => push(Tok::Cmp(Op::Lt), 1, &mut i),
            '>' if next == Some('=') => push(Tok::Cmp(Op::Ge), 2, &mut i),
            '>' => push(Tok::Cmp(Op::Gt), 1, &mut i),
            '~' => push(Tok::Cmp(Op::Match), 1, &mut i),
            '"' | '\'' => {
                let mut j = i + 1;
                let mut s = String::new();
                loop {
                    match chars.get(j) {
                        None => {
                            return err(
                                col,
                                format!("the string starting here has no closing {c}; add one"),
                            );
                        }
                        Some(&q) if q == c => break,
                        Some('\\') if chars.get(j + 1) == Some(&c) => {
                            s.push(c);
                            j += 1;
                        }
                        Some(&ch) => s.push(ch),
                    }
                    j += 1;
                }
                push(Tok::Str(s), j + 1 - i, &mut i);
            }
            _ if c.is_ascii_digit()
                || ((c == '.' || c == '-' || c == '+')
                    && next.is_some_and(|n| n.is_ascii_digit() || (c != '.' && n == '.'))) =>
            {
                let mut j = i + 1;
                while chars
                    .get(j)
                    .is_some_and(|d| d.is_ascii_digit() || *d == '.')
                {
                    j += 1;
                }
                if matches!(chars.get(j), Some('e' | 'E')) {
                    let mut k = j + 1;
                    if matches!(chars.get(k), Some('+' | '-')) {
                        k += 1;
                    }
                    if chars.get(k).is_some_and(|d| d.is_ascii_digit()) {
                        while chars.get(k).is_some_and(|d| d.is_ascii_digit()) {
                            k += 1;
                        }
                        j = k;
                    }
                }
                let text: String = chars[i..j].iter().collect();
                let Ok(n) = text.parse::<f64>() else {
                    return err(col, format!("{text:?} is not a valid number"));
                };
                push(Tok::Num(n), j - i, &mut i);
            }
            '.' => push(Tok::Dot, 1, &mut i),
            _ if is_ident_start(c) => {
                let mut j = i + 1;
                while chars.get(j).is_some_and(|&d| is_ident_char(d)) {
                    j += 1;
                }
                let text: String = chars[i..j].iter().collect();
                push(Tok::Ident(text), j - i, &mut i);
            }
            _ => {
                return err(
                    col,
                    format!("unexpected character {c:?}; remove it or quote it as a string"),
                );
            }
        }
    }
    out.push(Spanned {
        tok: Tok::End,
        col: chars.len() + 1,
    });
    Ok(out)
}
