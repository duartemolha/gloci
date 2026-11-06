//! gloci-core: pure rules (coordinates, chromosomes, genome, regions).

pub mod chrom;
pub mod distance;
pub mod error;
pub mod fmt;
pub mod genome;
pub mod interval;
pub mod overlap;
pub mod region;

pub use chrom::{ChromId, ChromTable, OrderSource};
pub use distance::{SignFrame, distance, signed_distance};
pub use error::{CoordError, FractionError, GenomeError, RegionError, UnknownChrom};
pub use fmt::thousands;
pub use genome::Genome;
pub use interval::{Interval, Strand, one_based_closed, parse_coord};
pub use overlap::{Fraction, OverlapCriteria, StrandRule, Window, overlap_len, overlaps};
pub use region::Region;
