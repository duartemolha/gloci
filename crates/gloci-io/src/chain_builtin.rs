//! Built-in chains for `--chain` (D174, R13), resolved by D33's rule.

use crate::builtin::edit_distance;
use crate::chain::{ChainIndex, MtRule};
use crate::error::{IoError, UsageError};
use gloci_chains::BuiltinChain;
use std::path::Path;

fn all() -> impl Iterator<Item = &'static BuiltinChain> {
    gloci_chains::builtin_chains()
        .iter()
        .chain(gloci_chains_t2t::builtin_chains())
}

/// Every built-in chain id, comma-separated, in listing order.
pub fn builtin_chain_ids() -> String {
    all().map(|c| c.id).collect::<Vec<_>>().join(", ")
}

/// The built-in id closest to `name` (edit distance 2 or less, ignoring case).
fn closest(name: &str) -> Option<&'static str> {
    let name = name.to_lowercase();
    all()
        .map(|c| (edit_distance(&name, &c.id.to_lowercase()), c.id))
        .filter(|(d, _)| *d <= 2)
        .min_by_key(|(d, _)| *d)
        .map(|(_, id)| id)
}

/// Only GRCh37-to-GRCh38 refuses a mitochondrial name. Every other built-in chain names it as a
/// chain file does: GRCh38-to-GRCh37 writes GRCh37's `MT` as `chrMT` for a `chr` input, never as
/// hg19's `chrM`, and the T2T chains have no mitochondrial chain (R13).
fn mt_rule(c: &BuiltinChain) -> MtRule {
    match (c.from, c.to) {
        ("GRCh37", "GRCh38") => MtRule::Hg19,
        _ => MtRule::Plain,
    }
}

impl ChainIndex {
    /// Read the chain `arg` names. A value with a `/` in it, or one that names an existing file
    /// (a directory does not count), is a file; so are `-` (stdin) and a host's in-memory source
    /// (D84). Otherwise it is a built-in chain id, ignoring case (D33's rule).
    pub fn resolve(arg: &str) -> Result<ChainIndex, IoError> {
        let is_file = Path::new(arg).metadata().is_ok_and(|m| !m.is_dir());
        if arg == "-" || crate::memio::is_mem_path(arg) || arg.contains('/') || is_file {
            return Self::load(arg);
        }
        let Some(c) = all().find(|c| c.id.eq_ignore_ascii_case(arg)) else {
            let hint = closest(arg)
                .map(|c| format!(" (did you mean {c}?)"))
                .unwrap_or_default();
            return Err(IoError::Usage(UsageError(format!(
                "no such file or built-in chain: \"{arg}\"{hint}; --chain takes a chain file \
(plain, gzip or BGZF) or one of the built-in chains: {}",
                builtin_chain_ids()
            ))));
        };
        let mut chain = Self::parse(c.id, c.gz)?;
        chain.mt = mt_rule(c);
        Ok(chain)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_builtin_loads_with_its_rule() {
        let mut n = 0;
        for c in all() {
            let chain = ChainIndex::resolve(&c.id.to_uppercase()).unwrap();
            assert!(chain.chains() > 20 && chain.blocks() > 10_000, "{}", c.id);
            assert_eq!(chain.mt_rule(), mt_rule(c), "{}", c.id);
            assert!(chain.source("chr1").is_some() || chain.source("1").is_some());
            n += 1;
        }
        assert_eq!(n, 6);
    }

    #[test]
    fn rules_follow_the_source_assembly() {
        let rule = |id: &str| ChainIndex::resolve(id).unwrap().mt_rule();
        assert_eq!(rule("GRCh37-to-GRCh38"), MtRule::Hg19);
        for id in ["GRCh38-to-GRCh37", "CHM13v2-to-GRCh38", "GRCh37-to-CHM13v2"] {
            assert_eq!(rule(id), MtRule::Plain, "{id}");
        }
    }

    #[test]
    fn a_missing_file_is_a_usage_error_naming_the_ids() {
        let Err(IoError::Usage(UsageError(m))) = ChainIndex::resolve("nope.chain") else {
            panic!("expected a usage error");
        };
        assert!(
            m.contains("GRCh37-to-GRCh38") && !m.contains("did you mean"),
            "{m}"
        );
    }
}
