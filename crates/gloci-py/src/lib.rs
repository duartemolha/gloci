//! `gloci._native`: the Python library's glue over [`gloci::embed`] (D83).

#[cfg(windows)]
compile_error!(
    "gloci does not build on Windows: use WSL (the Windows Subsystem for Linux) and install \
     gloci there with `pip install gloci`"
);

use pyo3::prelude::*;

#[pymodule]
fn _native(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
