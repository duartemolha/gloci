//! gloci-io: format detection and input opening.

pub mod atomic;
pub mod detect;
pub mod error;
pub mod input;
pub mod output;

pub use atomic::{
    AtomicFile, TempGuard, cleanup_registered, install_cleanup_handlers, register_temp,
};
pub use detect::{Compression, Format, detect_compression, detect_format};
pub use error::{IoError, Loc, UsageError};
pub use input::{Input, InputOptions, open_input, open_input_from};
pub use output::{
    OutFormat, Output, OutputOptions, OutputSpec, check_level, open_output, parse_output_spec,
    refuse_binary_to_tty,
};
