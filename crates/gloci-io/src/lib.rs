//! gloci-io: format detection, input opening and record readers.

pub mod atomic;
pub mod bcf;
pub mod bcf_encode;
pub mod bed;
pub mod builtin;
pub mod chrom_map;
pub mod convert;
pub mod detect;
pub mod diag;
pub mod error;
pub mod fields;
#[doc(hidden)]
pub mod fuzz_bodies;
pub mod genome;
pub mod gff;
pub mod index;
pub mod input;
pub mod naming;
pub mod output;
pub mod record;
pub mod sorted;
pub mod tsv;
pub mod vcf;
pub mod wig;
pub mod writer;

pub use atomic::{
    __test_with_registry_locked, AtomicFile, TempGuard, cleanup_registered,
    create_file_in_temp_dir, create_temp_dir_in, create_temp_file, install_cleanup_handlers,
    internal_error_line,
};
pub use bcf::{BcfDicts, read_bcf_header, write_bcf_header};
pub use bcf_encode::encode_gt;
pub use bed::{bed12_blocks, is_bed_record_line};
pub use builtin::{
    BuiltinGenome, Naming, builtin, builtin_alias, builtin_genomes, resolve_genome, unknown_genome,
};
pub use chrom_map::ChromMap;
pub use convert::{BcfText, HeaderMode, TextToBcf, convert_variants};
pub use detect::{Compression, Format, detect_compression, detect_format};
pub use diag::{Diagnostics, InvalidPolicy};
pub use error::{IoError, Loc, OrderError, SORTED_INPUTS, UsageError};
pub use fields::{column_count, field, split_fields};
pub use genome::{chrom_table_for, chrom_table_lenient, load_genome};
pub use index::{
    IndexInfo, IndexKind, IndexedFile, TsvColumns, file_chrom_name, find_index,
    index_non_empty_reference_names, index_reference_names, open_region_reader, write_index,
    write_tsv_index,
};
pub use input::{
    BGZF_EOF, Input, InputOptions, InputSpool, is_regular_file, open_input, open_input_from,
};
pub use naming::{
    chr_prefix_hint, format_chrom_list, format_disjoint_chroms_warning, nary_naming_warning,
    unknown_chrom_err,
};
pub use output::{
    OutFormat, Output, OutputOptions, OutputSpec, check_level, open_output, parse_output_spec,
    refuse_binary_to_tty,
};
pub use record::{
    Header, ReaderOptions, RecordBatch, RecordMeta, RecordReader, next_batch_size, open_reader,
    open_reader_from,
};
pub use sorted::{LEX_SORTED_NOTE, OrderGuard, sort_hint};
pub use tsv::{TsvHeader, open_tsv_writer, tsv_output_ok};
pub use vcf::{VcfHeader, shell_quote_args};
pub use writer::{
    RecordWriter, WriterOptions, check_family, family_mismatch, open_headerless_writer,
    open_record_writer,
};
