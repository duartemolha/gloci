// The extension module links no libpython (the interpreter that loads it provides the
// symbols). maturin passes the linker arguments that allow that on macOS; this does the same
// for a plain `cargo build` or `cargo clippy --workspace`.
fn main() {
    pyo3_build_config::add_extension_module_link_args();
}
