#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// mimalloc's `mi_option_purge_delay`, which `libmimalloc-sys` leaves out of its constants.
/// mimalloc keeps deprecated options in its enum so that the numbers never move; the test
/// below checks the neighbours `libmimalloc-sys` does name.
const MI_OPTION_PURGE_DELAY: libmimalloc_sys::mi_option_t = 15;

/// How long mimalloc keeps freed memory before returning it to the OS, in milliseconds:
/// mimalloc's default of 1000 is longer than most runs, so every buffer freed during a run
/// stayed resident to its end (8-thread BGZF intersect 118 -> 81 MiB, flagstat of a CRAM
/// 20 -> 16 MiB at 10 ms, at the same speed; at 0 the CRAM decoder's churn cost 57 %).
/// `MIMALLOC_PURGE_DELAY` in the environment still wins.
const PURGE_DELAY_MS: std::ffi::c_long = 10;

fn main() -> std::process::ExitCode {
    if std::env::var_os("MIMALLOC_PURGE_DELAY").is_none() {
        // SAFETY: mi_option_set is not thread safe; this is the first thing main does, before
        // any thread is started.
        unsafe { libmimalloc_sys::mi_option_set(MI_OPTION_PURGE_DELAY, PURGE_DELAY_MS) };
    }
    gloci::run(std::env::args_os().collect())
}

#[cfg(test)]
mod tests {
    /// In mimalloc's enum, `purge_delay` follows `reserve_os_memory` (9) and five deprecated
    /// options, and precedes `use_numa_nodes`.
    #[test]
    fn purge_delay_sits_between_the_options_libmimalloc_sys_names() {
        assert_eq!(libmimalloc_sys::mi_option_reserve_os_memory, 9);
        assert_eq!(
            libmimalloc_sys::mi_option_use_numa_nodes - 1,
            super::MI_OPTION_PURGE_DELAY
        );
    }
}
