//! Bioconda recipes (D93): `packaging/bioconda/<recipe>/meta.yaml` are templates. The
//! committed version line names the workspace version (`release-check` keeps it so) and
//! the checksum is a placeholder; `render` fills both in from the release.
use anyhow::{Context, Result, anyhow, bail};
use std::path::Path;

/// The recipes, in build order (the library recipe does not depend on the CLI's).
pub const RECIPES: [&str; 2] = ["gloci", "py-gloci"];
const SHA_PLACEHOLDER: &str = "@SHA256@";
const VERSION_LINE: &str = "{% set version = \"";

/// The version a recipe's `{% set version = "X" %}` line carries.
pub fn recipe_version(text: &str) -> Option<&str> {
    text.lines().find_map(|l| {
        l.strip_prefix(VERSION_LINE)?
            .split_once('"')
            .map(|(v, _)| v)
    })
}

fn is_version(v: &str) -> bool {
    !v.is_empty()
        && v.split('.').count() == 3
        && v.split('.')
            .all(|p| !p.is_empty() && p.bytes().all(|b| b.is_ascii_digit()))
}

/// Render one recipe. With `local_source`, the `url` and `sha256` of the source become
/// `path: <local_source>` (a directory), and `sha256` is not needed.
pub fn render(
    recipe: &Path,
    version: &str,
    sha256: &str,
    local_source: Option<&Path>,
) -> Result<String> {
    if !is_version(version) {
        bail!(
            "version `{version}` must be X.Y.Z (digits only); pass the release version without the leading v"
        );
    }
    if local_source.is_none()
        && !(sha256.len() == 64 && sha256.bytes().all(|b| b.is_ascii_hexdigit()))
    {
        bail!(
            "sha256 `{sha256}` must be 64 hex digits (the checksum of the recipe's source archive)"
        );
    }
    let text =
        std::fs::read_to_string(recipe).with_context(|| format!("reading {}", recipe.display()))?;
    if recipe_version(&text).is_none() {
        bail!(
            "{}: no `{{% set version = \"X.Y.Z\" %}}` line to render; add it at the top",
            recipe.display()
        );
    }
    let mut out = String::with_capacity(text.len() + 64);
    let mut in_source = false;
    for line in text.lines() {
        if line.starts_with(VERSION_LINE) {
            out.push_str(&format!("{VERSION_LINE}{version}\" %}}\n"));
            continue;
        }
        if !line.starts_with(' ') && !line.is_empty() {
            in_source = line == "source:";
        }
        let key = line.trim_start();
        if let (true, Some(dir)) = (in_source, local_source) {
            if key.starts_with("url:") {
                out.push_str(&format!("  path: {}\n", dir.display()));
                continue;
            }
            if key.starts_with("sha256:") {
                continue;
            }
        }
        out.push_str(&line.replace(SHA_PLACEHOLDER, sha256));
        out.push('\n');
    }
    if out.contains('@') {
        let at = out.lines().position(|l| l.contains('@')).unwrap() + 1;
        bail!(
            "{}: line {at} still has an `@` placeholder after rendering",
            recipe.display()
        );
    }
    Ok(out)
}

/// Render every recipe under `src` (the `packaging/bioconda` directory) into `out`,
/// copying the recipe's other files (build.sh, test/) alongside. With `local`, a
/// directory holding `gloci/` and `py-gloci/` source trees, no checksums are needed.
pub fn render_all(
    src: &Path,
    out: &Path,
    version: &str,
    sha256_cli: Option<&str>,
    sha256_py: Option<&str>,
    local: Option<&Path>,
) -> Result<()> {
    for (name, sha) in RECIPES.into_iter().zip([sha256_cli, sha256_py]) {
        let sha = match (sha, local) {
            (Some(s), _) => s,
            (None, Some(_)) => "",
            (None, None) => bail!(
                "no sha256 for the {name} recipe; pass --sha256-cli and --sha256-py (or --local DIR to build from source trees)"
            ),
        };
        let local_src = local.map(|d| d.join(name));
        let text = render(
            &src.join(name).join("meta.yaml"),
            version,
            sha,
            local_src.as_deref(),
        )?;
        let dest = out.join(name);
        std::fs::create_dir_all(&dest).with_context(|| format!("creating {}", dest.display()))?;
        std::fs::write(dest.join("meta.yaml"), text)?;
        copy_tree(&src.join(name), &dest)?;
    }
    Ok(())
}

/// Copy everything but `meta.yaml` (already rendered).
fn copy_tree(from: &Path, to: &Path) -> Result<()> {
    for e in std::fs::read_dir(from).with_context(|| format!("reading {}", from.display()))? {
        let e = e?;
        let (p, name) = (e.path(), e.file_name());
        if name == "meta.yaml" {
            continue;
        }
        let dest = to.join(&name);
        if e.file_type()?.is_dir() {
            std::fs::create_dir_all(&dest)?;
            copy_tree(&p, &dest)?;
        } else {
            std::fs::copy(&p, &dest)
                .map_err(|err| anyhow!("copying {} to {}: {err}", p.display(), dest.display()))?;
        }
    }
    Ok(())
}
