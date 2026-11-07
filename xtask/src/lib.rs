//! Developer tasks for gloci.
pub mod bench;
pub mod diff;
pub mod docs;
pub mod genomes;
// `gen` is a reserved keyword in edition 2024, so the module is named `synth`.
#[path = "gen.rs"]
pub mod synth;
