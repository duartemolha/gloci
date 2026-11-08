//! Sample selection for `vcf filter`.

use gloci_io::{IoError, UsageError};

/// Which samples of the header to keep, in header order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleSelection {
    pub keep: Vec<bool>,
}

impl SampleSelection {
    /// True if every sample is kept.
    pub fn all_kept(&self) -> bool {
        self.keep.iter().all(|&k| k)
    }

    /// Count of kept samples.
    pub fn count(&self) -> usize {
        self.keep.iter().filter(|&&k| k).count()
    }
}

fn unknown_sample_err(unknown: &str, header_samples: &[Vec<u8>]) -> UsageError {
    let first_5: Vec<_> = header_samples
        .iter()
        .take(5)
        .map(|s| String::from_utf8_lossy(s).into_owned())
        .collect();
    if first_5.is_empty() {
        UsageError(format!(
            "unknown sample \"{unknown}\"; the header has no samples"
        ))
    } else {
        UsageError(format!(
            "unknown sample \"{unknown}\"; first samples in header: {}",
            first_5.join(", ")
        ))
    }
}

/// Compute which samples to keep given inclusion and exclusion lists.
///
/// An unknown sample name in either list is a usage error naming the first
/// unknown sample and listing up to five sample names from the header.
pub fn select_samples(
    header_samples: &[Vec<u8>],
    include: Option<&[String]>,
    exclude: &[String],
) -> Result<SampleSelection, UsageError> {
    if let Some(inc) = include {
        for name in inc {
            if !header_samples
                .iter()
                .any(|s| s.as_slice() == name.as_bytes())
            {
                return Err(unknown_sample_err(name, header_samples));
            }
        }
    }
    for name in exclude {
        if !header_samples
            .iter()
            .any(|s| s.as_slice() == name.as_bytes())
        {
            return Err(unknown_sample_err(name, header_samples));
        }
    }

    let mut keep: Vec<bool> = match include {
        Some(inc) => header_samples
            .iter()
            .map(|s| inc.iter().any(|name| name.as_bytes() == s.as_slice()))
            .collect(),
        None => vec![true; header_samples.len()],
    };

    for (i, s) in header_samples.iter().enumerate() {
        if exclude.iter().any(|name| name.as_bytes() == s.as_slice()) {
            keep[i] = false;
        }
    }

    Ok(SampleSelection { keep })
}

/// Read sample names from a file, one per line. Blank lines and `#` comments
/// are skipped.
pub fn read_samples(path: &str) -> Result<Vec<String>, IoError> {
    let mut out = Vec::new();
    super::lists::each_line(path, |_, line| {
        out.push(String::from_utf8_lossy(line).into_owned());
        Ok(())
    })?;
    Ok(out)
}
