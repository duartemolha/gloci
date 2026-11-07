use std::path::Path;
use std::process::ExitCode;

fn main() -> ExitCode {
    let task = std::env::args().nth(1);
    let res = match task.as_deref() {
        Some("man") => xtask::docs::generate_man(Path::new("target/assets/man")),
        Some("completions") => {
            xtask::docs::generate_completions(Path::new("target/assets/completions"))
        }
        _ => {
            eprintln!("usage: cargo xtask <man|completions>");
            return ExitCode::from(2);
        }
    };
    match res {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e:#}");
            ExitCode::FAILURE
        }
    }
}
