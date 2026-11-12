//! `cargo xtask paper bundle`: the Zenodo evidence bundle (D98), a gzipped
//! tar that is a pure function of the results and the repository files:
//! entries sorted by path, mtime = `finished_at`, uid and gid 0, fixed modes
//! (0644 files, 0755 directories) and a gzip header with mtime 0.
use std::collections::BTreeMap;
use std::io::Write;
use std::path::Path;

use anyhow::{Context, Result, anyhow};

use crate::diff::workspace_root;

/// Repository files that go in beside the results, by path.
const REPO_FILES: [&str; 5] = [
    "xtask/bench/datasets.toml",
    "paper/fetch-datasets.sh",
    "xtask/docker/Dockerfile",
    "xtask/docker/Dockerfile.bench",
    "CITATION.cff",
];

fn collect(dir: &Path, prefix: &str, out: &mut BTreeMap<String, Vec<u8>>) -> Result<()> {
    for e in std::fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
        let e = e?;
        let p = e.path();
        let name = format!("{prefix}{}", e.file_name().to_string_lossy());
        if p.is_dir() {
            collect(&p, &format!("{name}/"), out)?;
        } else {
            out.insert(name, std::fs::read(&p)?);
        }
    }
    Ok(())
}

fn readme(version: &str, commit: &str, finished: &str) -> String {
    format!(
        "gloci paper evidence bundle\n\
         ===========================\n\
         \n\
         gloci version: {version}\n\
         gloci commit:  {commit}\n\
         finished at:   {finished}\n\
         \n\
         results/<version>/   raw benchmark JSON (hyperfine and peak RSS per case),\n\
         \x20                    concordance.csv and concordance/, environment.json (host,\n\
         \x20                    tool versions, image digests, wheel) and repo-snapshot.json\n\
         tables/              T1 to T8 as CSV, computed from results/ by `cargo xtask paper`\n\
         numbers.tex, .json   every number the paper quotes, keyed (\\glnum{{key}})\n\
         xtask/bench/datasets.toml, paper/fetch-datasets.sh\n\
         \x20                    the datasets with their provider checksums, and the script\n\
         \x20                    that fetches and derives them\n\
         xtask/docker/        the Dockerfiles of the pinned reference tools\n\
         xtask/paper/*.toml   the benchmark and concordance cases\n\
         CITATION.cff         how to cite gloci\n\
         \n\
         To reproduce: check out the commit, run `cargo xtask paper fetch`, then\n\
         `cargo xtask paper run`, `cargo xtask paper` and `cargo xtask paper bundle`\n\
         on a machine like the reference one (8-core x86_64, NVMe, 32 GB). The\n\
         bundle is byte-reproducible from its inputs: entries are sorted, mtimes are\n\
         the finish time, owners are 0 and the gzip header carries no time.\n"
    )
}

/// Write the bundle of `results` (a `paper/results/<version>` directory) to
/// the file `out`.
pub fn bundle(results: &Path, out: &Path) -> Result<()> {
    let root = workspace_root();
    let r = super::data::load(results)?;
    let finished = r.environment["finished_at"]
        .as_str()
        .ok_or_else(|| anyhow!("environment.json has no finished_at (written by `paper run`)"))?;
    let mtime = super::data::parse_utc(finished)?;
    let commit = r.environment["gloci_commit"].as_str().unwrap_or("unknown");

    let mut files: BTreeMap<String, Vec<u8>> = BTreeMap::new();
    let mut res = BTreeMap::new();
    collect(results, "", &mut res)?;
    for (k, v) in res {
        files.insert(format!("results/{}/{k}", r.version), v);
    }
    for (k, v) in super::tables::render(results)? {
        files.insert(k, v.into_bytes());
    }
    for f in REPO_FILES {
        files.insert(
            f.to_string(),
            std::fs::read(root.join(f)).with_context(|| format!("reading {f}"))?,
        );
    }
    let mut tomls: Vec<_> = std::fs::read_dir(root.join("xtask/paper"))?
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .filter(|n| n.ends_with(".toml"))
        .collect();
    tomls.sort();
    for n in tomls {
        let rel = format!("xtask/paper/{n}");
        files.insert(rel.clone(), std::fs::read(root.join(&rel))?);
    }
    files.insert(
        "README.txt".into(),
        readme(&r.version, commit, finished).into_bytes(),
    );

    // Directory entries, then everything by path.
    let top = format!("gloci-paper-evidence-{}", r.version);
    let mut entries: BTreeMap<String, Option<Vec<u8>>> = BTreeMap::new();
    entries.insert(format!("{top}/"), None);
    for (k, v) in files {
        let mut at = 0;
        while let Some(i) = k[at..].find('/') {
            at += i + 1;
            entries.insert(format!("{top}/{}", &k[..at]), None);
        }
        entries.insert(format!("{top}/{k}"), Some(v));
    }

    let mut gz = flate2::GzBuilder::new()
        .mtime(0)
        .operating_system(3)
        .write(Vec::new(), flate2::Compression::new(9));
    {
        let mut tar = tar::Builder::new(&mut gz);
        for (path, data) in &entries {
            let mut h = tar::Header::new_gnu();
            h.set_mtime(mtime);
            h.set_uid(0);
            h.set_gid(0);
            match data {
                Some(d) => {
                    h.set_entry_type(tar::EntryType::Regular);
                    h.set_mode(0o644);
                    h.set_size(d.len() as u64);
                    tar.append_data(&mut h, path, d.as_slice())?;
                }
                None => {
                    h.set_entry_type(tar::EntryType::Directory);
                    h.set_mode(0o755);
                    h.set_size(0);
                    tar.append_data(&mut h, path, std::io::empty())?;
                }
            }
        }
        tar.finish()?;
    }
    let bytes = gz.finish()?;

    if let Some(d) = out.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(d)?;
    }
    let tmp = out.with_extension("tmp");
    std::fs::File::create(&tmp)?.write_all(&bytes)?;
    std::fs::rename(&tmp, out).with_context(|| format!("writing {}", out.display()))?;
    Ok(())
}
