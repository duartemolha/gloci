//! Developer tasks for gloci.
pub mod bench;
pub mod book;
pub mod conda;
pub mod diff;
pub mod docs;
pub mod genomes;
pub mod paper;
pub mod release;
// `gen` is a reserved keyword in edition 2024, so the module is named `synth`.
#[path = "gen.rs"]
pub mod synth;
