//! gloci-io: format detection, input opening and record readers.

pub mod atomic;
pub mod bcf;
pub mod bed;
pub mod builtin;
pub mod convert;
pub mod detect;
pub mod diag;
pub mod error;
#[doc(hidden)]
pub mod fuzz_bodies;
pub mod genome;
pub mod gff;
pub mod index;
pub mod input;
pub mod output;
pub mod record;
pub mod sorted;
pub mod vcf;
pub mod writer;

pub use atomic::{
    AtomicFile, TempGuard, cleanup_registered, install_cleanup_handlers, register_temp,
};
pub use bcf::{BcfDicts, read_bcf_header, write_bcf_header};
pub use bed::{bed12_blocks, is_bed_record_line};
pub use builtin::{
    BuiltinGenome, Naming, builtin, builtin_genomes, resolve_genome, unknown_genome,
};
pub use convert::{HeaderMode, convert_variants};
pub use detect::{Compression, Format, detect_compression, detect_format};
pub use diag::{Diagnostics, InvalidPolicy};
pub use error::{IoError, Loc, OrderError, UsageError};
pub use genome::{chrom_table_for, chrom_table_lenient, load_genome};
pub use index::{
    IndexInfo, IndexKind, find_index, index_reference_names, open_region_reader, write_index,
};
pub use input::{Input, InputOptions, open_input, open_input_from};
pub use output::{
    OutFormat, Output, OutputOptions, OutputSpec, check_level, open_output, parse_output_spec,
    refuse_binary_to_tty,
};
pub use record::{
    Header, ReaderOptions, RecordBatch, RecordMeta, RecordReader, next_batch_size, open_reader,
    open_reader_from,
};
pub use sorted::{LEX_SORTED_NOTE, OrderGuard, sort_hint};
pub use vcf::{VcfHeader, shell_quote_args};
pub use writer::{
    RecordWriter, WriterOptions, check_family, family_mismatch, open_headerless_writer, open_record_writer,
};
