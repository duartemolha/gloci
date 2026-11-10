//! All-or-nothing publication of several output files.

use crate::{AtomicFile, IoError};
use std::ffi::OsString;
use std::io;
use std::path::PathBuf;

/// Files that have been written and flushed but not yet renamed into place.
/// A writer "commits" by handing its file here; [`StagedSet::publish`]
/// renames them all once every statistic has finished.
#[derive(Default)]
pub struct StagedSet {
    files: Vec<(AtomicFile, PathBuf)>,
}

impl StagedSet {
    /// Stage `file` to be renamed onto `dest` by [`StagedSet::publish`].
    pub fn push(&mut self, file: AtomicFile, dest: PathBuf) {
        self.files.push((file, dest));
    }

    /// Rename every staged file onto its destination. A destination that
    /// already exists (the previous run's output) is first moved aside to a
    /// registered temp name in the same directory. When every file is in
    /// place the backups are deleted; if a step fails, the files published in
    /// this run are removed and the backups put back, so a failed run leaves
    /// the previous outputs as they were and no partial set of new ones.
    pub fn publish(self) -> Result<Vec<PathBuf>, IoError> {
        // A directory on a destination cannot be replaced; say so before
        // anything is moved.
        for (_, dest) in &self.files {
            if std::fs::symlink_metadata(dest).is_ok_and(|m| m.is_dir()) {
                let d = dest.display().to_string();
                return Err(IoError::Detect {
                    path: d.clone(),
                    msg: format!(
                        "cannot write {d}: a directory with that name exists; \
remove or rename it, or choose another --out-prefix"
                    ),
                });
            }
        }
        let mut done: Vec<PathBuf> = Vec::new();
        // (destination, where its previous content waits)
        let mut backups: Vec<(PathBuf, PathBuf)> = Vec::new();
        for (file, dest) in self.files {
            let step = (|| -> io::Result<()> {
                if std::fs::symlink_metadata(&dest).is_ok() {
                    let backup = back_up(
                        &dest,
                        |a, b| std::fs::hard_link(a, b),
                        |a, b| std::fs::rename(a, b),
                    )?;
                    backups.push((dest.clone(), backup));
                }
                file.commit()
            })();
            if let Err(e) = step {
                let kept = roll_back(&done, &backups);
                let err = IoError::Write {
                    path: dest.display().to_string(),
                    source: e,
                };
                return Err(if kept.is_empty() {
                    err
                } else {
                    IoError::Detect {
                        path: dest.display().to_string(),
                        msg: format!(
                            "{err}; the previous outputs could not all be put back: {kept}"
                        ),
                    }
                });
            }
            done.push(dest);
        }
        for (_, backup) in &backups {
            let _ = std::fs::remove_file(backup);
        }
        Ok(done)
    }
}

/// Keep the previous content of `dest` under a registered name next to it.
/// `dest` stays where it is, as a hard link, so the rename that publishes the
/// new file replaces it atomically and a signal never finds it missing. When
/// the link cannot be made (a file system without hard links) the file is
/// moved aside instead; if that fails too, the reserved name is removed.
fn back_up(
    dest: &std::path::Path,
    link: impl Fn(&std::path::Path, &std::path::Path) -> io::Result<()>,
    rename: impl Fn(&std::path::Path, &std::path::Path) -> io::Result<()>,
) -> io::Result<PathBuf> {
    let backup = reserve_backup_name(dest)?;
    // The name is registered, so it stays reserved for this run; a link needs
    // it free.
    std::fs::remove_file(&backup)?;
    if link(dest, &backup).is_ok() {
        return Ok(backup);
    }
    rename(dest, &backup).inspect_err(|_| {
        let _ = std::fs::remove_file(&backup);
    })?;
    Ok(backup)
}

/// A new, registered, empty file next to `dest` to move its previous content
/// onto (the move replaces it).
fn reserve_backup_name(dest: &std::path::Path) -> io::Result<PathBuf> {
    let name = dest.file_name().unwrap_or_default();
    for n in 0..1000u32 {
        let mut b: OsString = name.to_owned();
        b.push(format!(".bak-{}-{n}", std::process::id()));
        let path = dest.with_file_name(b);
        match crate::create_temp_file(&path) {
            Ok(_) => return Ok(path),
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {}
            Err(e) => return Err(e),
        }
    }
    Err(io::Error::other("cannot create a backup file name"))
}

/// Remove the files published in this run and move the backups back. Returns
/// a description of the previous outputs that could not be restored (empty
/// when all were).
fn roll_back(done: &[PathBuf], backups: &[(PathBuf, PathBuf)]) -> String {
    for p in done {
        let _ = std::fs::remove_file(p);
    }
    let mut kept = Vec::new();
    for (dest, backup) in backups {
        if std::fs::rename(backup, dest).is_err() {
            kept.push(format!(
                "{} is kept as {}",
                dest.display(),
                backup.display()
            ));
        } else {
            // A hard-linked backup of a dest that was never replaced is the
            // same file: the rename keeps both names.
            let _ = std::fs::remove_file(backup);
        }
    }
    kept.join(", ")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn staged(dir: &std::path::Path, name: &str, body: &str) -> (AtomicFile, PathBuf) {
        let dest = dir.join(name);
        let mut f = AtomicFile::create(&dest).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        (f, dest)
    }

    fn names(dir: &std::path::Path) -> Vec<String> {
        let mut v: Vec<_> = std::fs::read_dir(dir)
            .unwrap()
            .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
            .collect();
        v.sort();
        v
    }

    #[test]
    fn publish_puts_the_previous_outputs_back_when_a_later_rename_fails() {
        let d = tempfile::tempdir().unwrap();
        std::fs::write(d.path().join("o.frq"), "previous frq").unwrap();
        std::fs::write(d.path().join("o.hwe"), "previous hwe").unwrap();
        let first = staged(d.path(), "o.frq", "new frq");
        let second = staged(d.path(), "o.hwe", "new hwe");
        // The second temp file vanishes, so its rename fails after the first
        // file was published.
        std::fs::remove_file(second.0.temp_path()).unwrap();
        let err = StagedSet {
            files: vec![first, second],
        }
        .publish()
        .unwrap_err();
        assert!(err.to_string().contains("o.hwe"), "{err}");
        assert_eq!(names(d.path()), vec!["o.frq", "o.hwe"]);
        assert_eq!(
            std::fs::read_to_string(d.path().join("o.frq")).unwrap(),
            "previous frq"
        );
        assert_eq!(
            std::fs::read_to_string(d.path().join("o.hwe")).unwrap(),
            "previous hwe"
        );
    }

    #[test]
    fn back_up_hard_links_and_leaves_the_destination_in_place() {
        let d = tempfile::tempdir().unwrap();
        let dest = d.path().join("o.frq");
        std::fs::write(&dest, "previous").unwrap();
        let backup = back_up(
            &dest,
            |a, b| std::fs::hard_link(a, b),
            |a, b| std::fs::rename(a, b),
        )
        .unwrap();
        assert_eq!(std::fs::read_to_string(&dest).unwrap(), "previous");
        assert_eq!(std::fs::read_to_string(&backup).unwrap(), "previous");
    }

    #[test]
    fn back_up_falls_back_to_a_move_and_cleans_up_when_that_fails() {
        let d = tempfile::tempdir().unwrap();
        let dest = d.path().join("o.frq");
        std::fs::write(&dest, "previous").unwrap();
        let no_link = |_: &std::path::Path, _: &std::path::Path| -> io::Result<()> {
            Err(io::Error::other("no links"))
        };
        let backup = back_up(&dest, no_link, |a, b| std::fs::rename(a, b)).unwrap();
        assert!(!dest.exists());
        assert_eq!(std::fs::read_to_string(&backup).unwrap(), "previous");
        std::fs::rename(&backup, &dest).unwrap();
        let no_move = |_: &std::path::Path, _: &std::path::Path| -> io::Result<()> {
            Err(io::Error::other("no move"))
        };
        back_up(&dest, no_link, no_move).unwrap_err();
        assert_eq!(names(d.path()), vec!["o.frq"]);
    }

    #[test]
    fn publish_removes_what_it_published_when_there_was_no_previous_output() {
        let d = tempfile::tempdir().unwrap();
        let first = staged(d.path(), "o.frq", "new frq");
        let second = staged(d.path(), "o.hwe", "new hwe");
        std::fs::remove_file(second.0.temp_path()).unwrap();
        StagedSet {
            files: vec![first, second],
        }
        .publish()
        .unwrap_err();
        assert!(names(d.path()).is_empty(), "{:?}", names(d.path()));
    }
}
