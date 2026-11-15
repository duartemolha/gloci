//! `gloci links`: an HTML page of genome browser links, one per feature.
use super::igv::{open_features, write_text};
use crate::Ctx;
use crate::cli::LinksArgs;
use gloci_core::ChromTable;
use gloci_ops::{LinksOptions, links};

pub fn run(ctx: &Ctx, a: &LinksArgs) -> anyhow::Result<()> {
    let reader = open_features(ctx, &a.file, a.format.format.as_deref())?;
    let o = LinksOptions {
        base: a.base.clone(),
        org: a.org.clone(),
        db: a.db.clone(),
        title: if a.file == "-" {
            "stdin".to_string()
        } else {
            gloci_io::memio::display_name(&a.file)
        },
    };
    let mut chroms = ChromTable::first_appearance();
    write_text(ctx, a.output.as_deref(), |w| {
        links(reader, &mut chroms, &o, w)
    })
}
