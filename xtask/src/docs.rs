//! Man page and shell completion generation.
use anyhow::Context;
use clap_complete::Shell;
use std::path::Path;

/// Write `gloci.1` and one `gloci-<sub>.1` per subcommand into `out`.
pub fn generate_man(out: &Path) -> anyhow::Result<()> {
    std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    clap_mangen::generate_to(gloci::command(), out).context("generating man pages")?;
    Ok(())
}

/// Write bash, zsh and fish completions into `out`.
pub fn generate_completions(out: &Path) -> anyhow::Result<()> {
    std::fs::create_dir_all(out).with_context(|| format!("creating {}", out.display()))?;
    for shell in [Shell::Bash, Shell::Zsh, Shell::Fish] {
        let mut cmd = gloci::command();
        clap_complete::generate_to(shell, &mut cmd, "gloci", out)
            .with_context(|| format!("generating {shell} completions"))?;
    }
    Ok(())
}
