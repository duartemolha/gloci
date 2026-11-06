//! gloci-io: format detection and input opening.

pub mod detect;
pub mod error;
pub mod input;

pub use detect::{Compression, Format, detect_compression, detect_format};
pub use error::{IoError, Loc};
pub use input::{Input, InputOptions, open_input, open_input_from};
