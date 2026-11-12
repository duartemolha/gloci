//! gloci-io: format detection, input opening and record readers.

pub mod align;
pub mod atomic;
pub mod bcf;
pub mod bcf_encode;
mod bcf_fields;
pub mod bed;
pub mod bedpe;
pub mod bgzf_copy;
pub mod builtin;
pub mod cancel;
pub mod chrom_map;
pub mod convert;
pub mod detect;
pub mod diag;
pub mod error;
pub mod fasta;
pub mod fields;
#[doc(hidden)]
pub mod fuzz_bodies;
pub mod genome;
pub mod gff;
pub mod index;
pub mod input;
pub mod memio;
pub mod naming;
pub mod output;
pub(crate) mod pool;
pub use pool::WORKER_STACK;
mod psl;
pub mod record;
pub mod rmsk;
pub mod seekable;
pub mod sorted;
pub mod staged;
pub mod tsv;
pub mod vcf;
pub mod wig;
pub mod writer;

pub use align::{
    Accept, AlignOptions, Aln, AlnBed, AlnBedLayout, AlnBedOptions, AlnBlocks, AlnHeader, AlnRef,
    AlnSortOrder, AlnSource, CramPrelude, CramReader, FlagFields, FlagFieldsVisitor, RecordVisitor,
    SAM_FLAGS, SamAln, SamToBam, SortOrder, alignment_query_span, aln_ref, aln_view,
    bam_header_blocks, blocks, blocks_within, check_cram_version, check_sq_compatible,
    check_sq_compatible_named, cram_eof_present, cram_first_container, flag_names,
    open_alignment_writer, open_aln_source, open_aln_source_from, open_aln_source_regions,
    parse_flags, read_bam_file_header, read_bam_header_with_extent, read_cram_header,
    read_cram_prelude, read_new_header, regions_source_from, sam_line_to_bam,
    warn_unused_read_filters, warn_unused_reference, write_cram_header_container,
};
pub use atomic::{
    __test_with_registry_locked, AtomicFile, TempGuard, cleanup_registered,
    create_file_in_temp_dir, create_temp_dir_in, create_temp_file, install_cleanup_handlers,
    internal_error_line,
};
pub use bcf::{BcfDicts, check_record as check_bcf_record, read_bcf_header, write_bcf_header};
pub use bcf_encode::encode_gt;
pub use bed::{bed12_blocks, is_bed_record_line};
pub use bedpe::{BedpeBatch, BedpeEnd, BedpeMeta, BedpeReader, bedpe_name};
pub use builtin::{
    BuiltinGenome, Naming, builtin, builtin_alias, builtin_genomes, resolve_genome, unknown_genome,
};
pub use chrom_map::ChromMap;
pub use convert::{BcfText, HeaderMode, TextToBcf, convert_variants};
pub use detect::{Compression, Format, detect_compression, detect_format};
pub use diag::{Diagnostics, InvalidPolicy};
pub use error::{GateKind, IoError, Loc, OrderError, SORTED_INPUTS, UsageError};
pub use fasta::{
    FETCH_CHUNK, FaiEntry, FastaLine, FastaStream, IndexedFasta, build_fai, fai_bytes, fai_path,
    find_fai, gzi_path, parse_fai, write_fasta_index,
};
pub use fields::{column_count, field, split_fields};
pub use genome::{chrom_table_for, chrom_table_lenient, load_genome};
pub use index::{
    AlnIndexCounts, AlnRegionCursor, CraiEntry, CramIndex, IndexInfo, IndexKind, IndexedFile,
    RefCounts, TsvColumns, crai_path, file_chrom_name, find_crai, find_index,
    index_non_empty_reference_names, index_reference_names, open_region_reader,
    write_alignment_index, write_crai, write_index, write_tsv_index,
};
pub use input::{
    BGZF_EOF, Input, InputOptions, InputSpool, is_regular_file, open_input, open_input_from,
    open_raw,
};
pub use naming::{
    SeenChroms, chr_prefix_hint, format_chrom_list, format_disjoint_chroms_warning,
    nary_naming_warning, unknown_chrom_err,
};
pub use output::{
    OutFormat, Output, OutputOptions, OutputSpec, check_level, open_output, open_raw_output,
    parse_output_spec, plain_text_output_spec, refuse_binary_to_tty,
};
pub use record::{
    Header, Opened, ReaderOptions, RecordBatch, RecordMeta, RecordReader, check_map_against_index,
    next_batch_size, open_reader, open_reader_from, open_reader_noting_empty,
    open_reader_or_alignments, open_reader_or_alignments_from,
    open_reader_or_alignments_from_nonempty,
};
pub use rmsk::rmsk_fields;
pub use seekable::{ChromRange, SeekableFile, chrom_ranges, open_range_reader};
pub use sorted::{LEX_SORTED_NOTE, OrderGuard, sort_hint};
pub use staged::StagedSet;
pub use tsv::{TsvHeader, open_tsv_writer, tsv_output_ok};
pub use vcf::{VcfHeader, shell_quote_args};
pub use writer::{
    RecordWriter, WriterOptions, check_family, family_mismatch, open_headerless_writer,
    open_record_writer,
};
