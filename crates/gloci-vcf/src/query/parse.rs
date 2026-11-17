//! The syntax of a query format (D172): text, escapes, `%NAME{N}` fields and
//! `[ ]` blocks. Names are bound later, by `resolve`.

use gloci_io::UsageError;

/// A syntax item.
#[derive(Clone, Debug)]
pub(super) enum Item {
    Text(Vec<u8>),
    Field(Spec),
    /// A `[ ]` block.
    Block(Vec<Item>),
}

/// A `%NAME` or `%NAME{N}`; `col` is the column of its `%`.
#[derive(Clone, Debug)]
pub(super) struct Spec {
    pub name: String,
    pub sub: Option<usize>,
    pub col: usize,
}

fn err(src: &str, msg: &str, col: usize) -> UsageError {
    UsageError(format!(
        "invalid query format \"{src}\": {msg} (column {col})"
    ))
}

fn is_name(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_' || c == '.'
}

pub(super) fn parse(src: &str) -> Result<Vec<Item>, UsageError> {
    let chars: Vec<char> = src.chars().collect();
    let mut pos = 0;
    parse_items(src, &chars, &mut pos, None)
}

fn push_text(items: &mut Vec<Item>, bytes: &[u8]) {
    if let Some(Item::Text(t)) = items.last_mut() {
        t.extend_from_slice(bytes);
    } else {
        items.push(Item::Text(bytes.to_vec()));
    }
}

/// Items up to the end of the format, or up to the `]` that closes a block.
fn parse_items(
    src: &str,
    chars: &[char],
    pos: &mut usize,
    open: Option<usize>,
) -> Result<Vec<Item>, UsageError> {
    let mut items = Vec::new();
    while *pos < chars.len() {
        let c = chars[*pos];
        let col = *pos + 1;
        match c {
            '\\' => {
                let next = chars.get(*pos + 1).copied();
                let byte = match next {
                    Some('t') => b'\t',
                    Some('n') => b'\n',
                    Some('\\') => b'\\',
                    Some(o) => {
                        return Err(err(
                            src,
                            &format!("unknown escape \\{o}; the escapes are \\t, \\n or \\\\"),
                            col,
                        ));
                    }
                    None => {
                        return Err(err(src, "the format ends with a backslash", col));
                    }
                };
                push_text(&mut items, &[byte]);
                *pos += 2;
            }
            '[' => {
                if open.is_some() {
                    return Err(err(src, "nested [ ] is not allowed", col));
                }
                *pos += 1;
                let inner = parse_items(src, chars, pos, Some(col))?;
                // `parse_items` stops at the closing `]` (consumed) or errors.
                items.push(Item::Block(inner));
            }
            ']' => {
                if open.is_none() {
                    return Err(err(src, "] without [", col));
                }
                *pos += 1;
                return Ok(items);
            }
            '%' => {
                items.push(Item::Field(parse_field(src, chars, pos)?));
            }
            _ => {
                let mut buf = [0u8; 4];
                push_text(&mut items, c.encode_utf8(&mut buf).as_bytes());
                *pos += 1;
            }
        }
    }
    if let Some(col) = open {
        return Err(err(src, "missing ]", col));
    }
    Ok(items)
}

fn parse_field(src: &str, chars: &[char], pos: &mut usize) -> Result<Spec, UsageError> {
    let col = *pos + 1;
    *pos += 1;
    let name_end = |from: usize| {
        let mut i = from;
        while i < chars.len() && is_name(chars[i]) {
            i += 1;
        }
        i
    };
    let end = name_end(*pos);
    if end == *pos {
        return Err(err(src, "% must be followed by a field name", col));
    }
    let mut name: String = chars[*pos..end].iter().collect();
    *pos = end;
    if matches!(name.as_str(), "INFO" | "FMT" | "FORMAT") && chars.get(*pos) == Some(&'/') {
        let end = name_end(*pos + 1);
        name.push('/');
        name.extend(&chars[*pos + 1..end]);
        *pos = end;
    }
    let mut sub = None;
    if chars.get(*pos) == Some(&'{') {
        let brace = *pos + 1;
        let mut i = *pos + 1;
        while i < chars.len() && chars[i].is_ascii_digit() {
            i += 1;
        }
        let digits: String = chars[*pos + 1..i].iter().collect();
        let closed = chars.get(i) == Some(&'}') && !digits.is_empty();
        match closed.then(|| digits.parse::<usize>()) {
            Some(Ok(n)) => {
                sub = Some(n);
                *pos = i + 1;
            }
            Some(Err(_)) => {
                return Err(err(
                    src,
                    &format!("%{name}{{{digits}}}: the subscript is too large"),
                    brace,
                ));
            }
            None => {
                return Err(err(
                    src,
                    &format!("%{name} is followed by `{{` but not by `{{N}}` (N a number from 0)"),
                    brace,
                ));
            }
        }
    }
    Ok(Spec { name, sub, col })
}
