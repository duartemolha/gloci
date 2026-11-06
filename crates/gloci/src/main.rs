use clap::Parser;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// Genomic interval toolkit.
#[derive(Parser)]
#[command(name = "gloci", version, arg_required_else_help = true)]
struct Cli {}

fn main() {
    let _cli = Cli::parse();
}
