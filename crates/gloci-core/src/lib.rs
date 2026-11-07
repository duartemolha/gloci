//! gloci-core: pure rules (coordinates, chromosomes, genome, regions).

pub mod chrom;
pub mod distance;
pub mod error;
pub mod fmt;
pub mod footprint;
pub mod genome;
pub mod interval;
pub mod order;
pub mod overlap;
pub mod region;

pub use chrom::{ChromId, ChromOrder, ChromTable, OrderSource, karyotype_cmp};
pub use distance::{SignFrame, distance, signed_distance};
pub use error::{
    CoordError, FootprintError, FractionError, GenomeError, RegionError, UnknownChrom,
};
pub use fmt::thousands;
pub use footprint::{Variant, VcfSpan, alt_footprints, footprint};
pub use genome::Genome;
pub use interval::{Interval, Strand, one_based_closed, parse_coord};
pub use order::{OrderChecker, OrderViolation};
pub use overlap::{Fraction, OverlapCriteria, StrandRule, Window, overlap_len, overlaps};
pub use region::Region;
