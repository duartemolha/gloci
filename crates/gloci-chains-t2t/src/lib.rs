//! Built-in T2T-CHM13v2 chain files for `gloci liftover` (D174, R13), from the T2T consortium's
//! `chain/v1_nflo` release, which is CC0. The Ensembl chains are in `gloci-chains`; this crate
//! is separate so that each stays under crates.io's 10 MB. `data/SOURCES.toml` records where each
//! file came from, its hash and its licence.

pub use gloci_chains::BuiltinChain;

static CHAINS: [BuiltinChain; 4] = [
    BuiltinChain {
        id: "GRCh38-to-CHM13v2",
        from: "GRCh38",
        to: "CHM13v2",
        naming: "chr1",
        licence: "CC0",
        gz: include_bytes!("../data/GRCh38-to-CHM13v2.chain.gz"),
    },
    BuiltinChain {
        id: "CHM13v2-to-GRCh38",
        from: "CHM13v2",
        to: "GRCh38",
        naming: "chr1",
        licence: "CC0",
        gz: include_bytes!("../data/CHM13v2-to-GRCh38.chain.gz"),
    },
    BuiltinChain {
        id: "GRCh37-to-CHM13v2",
        from: "GRCh37",
        to: "CHM13v2",
        naming: "chr1",
        licence: "CC0",
        gz: include_bytes!("../data/GRCh37-to-CHM13v2.chain.gz"),
    },
    BuiltinChain {
        id: "CHM13v2-to-GRCh37",
        from: "CHM13v2",
        to: "GRCh37",
        naming: "chr1",
        licence: "CC0",
        gz: include_bytes!("../data/CHM13v2-to-GRCh37.chain.gz"),
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
