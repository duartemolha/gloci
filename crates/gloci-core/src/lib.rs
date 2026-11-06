//! gloci-core: pure rules (coordinates, chromosomes, genome, regions).

pub mod chrom;
pub mod error;
pub mod fmt;
pub mod genome;
pub mod interval;
pub mod region;

pub use chrom::{ChromId, ChromTable, OrderSource};
pub use error::{CoordError, GenomeError, RegionError, UnknownChrom};
pub use fmt::thousands;
pub use genome::Genome;
pub use interval::{Interval, Strand, one_based_closed, parse_coord};
pub use region::Region;
