//! gloci-io: format detection, input opening and record readers.

pub mod atomic;
pub mod bcf;
pub mod bed;
pub mod detect;
pub mod diag;
pub mod error;
pub mod gff;
pub mod input;
pub mod output;
pub mod record;
pub mod vcf;

pub use atomic::{
    AtomicFile, TempGuard, cleanup_registered, install_cleanup_handlers, register_temp,
};
pub use bcf::{BcfDicts, read_bcf_header, write_bcf_header};
pub use bed::bed12_blocks;
pub use detect::{Compression, Format, detect_compression, detect_format};
pub use diag::{Diagnostics, InvalidPolicy};
pub use error::{IoError, Loc, UsageError};
pub use input::{Input, InputOptions, open_input, open_input_from};
pub use output::{
    OutFormat, Output, OutputOptions, OutputSpec, check_level, open_output, parse_output_spec,
    refuse_binary_to_tty,
};
pub use record::{
    Header, ReaderOptions, RecordBatch, RecordMeta, RecordReader, open_reader, open_reader_from,
};
pub use vcf::{VcfHeader, shell_quote_args};
