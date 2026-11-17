//! bigWig and bigBed (D171), read through each file's own index.

mod read;

pub use read::bbi_chrom_sizes;

/// The refusal of `gloci index` and of tabix/CSI reads for a bigWig or bigBed, which carries
/// its own index (D171).
pub fn bbi_unindexed(name: &str, format: crate::Format) -> crate::UsageError {
    crate::UsageError(format!(
        "{name} is {}, which carries its own index; gloci reads it directly",
        format.name()
    ))
}
pub(crate) use read::open_bbi;
