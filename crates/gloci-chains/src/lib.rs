//! Built-in human chain files for `gloci liftover` (D174, R13): GRCh37 and GRCh38 from Ensembl,
//! which imposes no restrictions on its data. The T2T-CHM13v2 chains are in `gloci-chains-t2t`
//! (each crate stays under crates.io's 10 MB). `data/SOURCES.toml` records where each file came
//! from, its hash and its licence.

/// One chain file built into gloci, kept as the gzip file it was committed as.
#[derive(Debug, Clone, Copy)]
pub struct BuiltinChain {
    /// What `--chain` takes, ignoring case: `GRCh37-to-GRCh38`.
    pub id: &'static str,
    /// The assembly the chain moves records from (its source, the chain's `t` side).
    pub from: &'static str,
    /// The assembly it moves them to.
    pub to: &'static str,
    /// How the source names its chromosomes: `1` or `chr1`.
    pub naming: &'static str,
    /// The licence its publisher states.
    pub licence: &'static str,
    /// The published chain file, gzip -9.
    pub gz: &'static [u8],
}

static CHAINS: [BuiltinChain; 2] = [
    BuiltinChain {
        id: "GRCh37-to-GRCh38",
        from: "GRCh37",
        to: "GRCh38",
        naming: "1 and chr1 (Ensembl lists both)",
        licence: "Ensembl: no restrictions",
        gz: include_bytes!("../data/GRCh37-to-GRCh38.chain.gz"),
    },
    BuiltinChain {
        id: "GRCh38-to-GRCh37",
        from: "GRCh38",
        to: "GRCh37",
        naming: "1 and chr1 (Ensembl lists both)",
        licence: "Ensembl: no restrictions",
        gz: include_bytes!("../data/GRCh38-to-GRCh37.chain.gz"),
    },
];

/// The chains of this crate, in listing order.
pub fn builtin_chains() -> &'static [BuiltinChain] {
    &CHAINS
}

#[cfg(test)]
mod tests {
    use super::*;
    use sha2::{Digest, Sha256};

    const SOURCES: &str = include_str!("../data/SOURCES.toml");

    #[test]
    fn sources_match_the_files() {
        let doc: toml::Table = SOURCES.parse().unwrap();
        let entries = doc["chain"].as_array().unwrap();
        assert_eq!(entries.len(), builtin_chains().len());
        for c in builtin_chains() {
            let e = entries
                .iter()
                .find(|e| e["id"].as_str() == Some(c.id))
                .unwrap_or_else(|| panic!("{}: no SOURCES.toml entry", c.id));
            let sha: String = Sha256::digest(c.gz)
                .iter()
                .map(|b| format!("{b:02x}"))
                .collect();
            assert_eq!(e["sha256"].as_str().unwrap(), sha, "{}: sha256", c.id);
            assert_eq!(
                e["licence"].as_str().unwrap(),
                c.licence,
                "{}: licence",
                c.id
            );
            for k in ["url", "retrieved", "naming"] {
                assert!(!e[k].as_str().unwrap().is_empty(), "{}: {k}", c.id);
            }
        }
    }

    #[test]
    fn the_files_stay_under_the_crate_limit() {
        let total: usize = builtin_chains().iter().map(|c| c.gz.len()).sum();
        assert!(total < 9 << 20, "{total} bytes; crates.io allows 10 MB");
    }

    #[test]
    fn every_chain_is_gzip() {
        for c in builtin_chains() {
            assert_eq!(&c.gz[..2], &[0x1f, 0x8b], "{}", c.id);
        }
    }
}
