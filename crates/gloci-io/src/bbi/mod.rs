//! bigWig and bigBed (D171), read through each file's own index.

mod read;

pub use read::bbi_chrom_sizes;
pub(crate) use read::open_bbi;
