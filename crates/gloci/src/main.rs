#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

/// mimalloc's `mi_option_purge_delay`, which `libmimalloc-sys` leaves out of its constants.
/// mimalloc keeps deprecated options in its enum so that the numbers never move; the tests
/// below check it.
const MI_OPTION_PURGE_DELAY: libmimalloc_sys::mi_option_t = 15;

/// How long mimalloc keeps freed memory before returning it to the OS, in milliseconds:
/// mimalloc's default of 1000 is longer than most runs, so every buffer freed during a run
/// stayed resident to its end (8-thread BGZF intersect 118 -> 80 MiB, flagstat of a CRAM
/// 20 -> 17 MiB at 10 ms, at the same speed; at 0 the CRAM decoder's churn cost 57 %).
const PURGE_DELAY_MS: std::ffi::c_long = 10;

/// Whether the environment sets the purge delay, under either name mimalloc reads and in any
/// case, as mimalloc matches them. mimalloc reads its options when the process loads, before
/// `main`, so `mi_option_set_default` could not tell an environment value from its default.
fn purge_delay_in_env() -> bool {
    std::env::vars_os().any(|(k, _)| {
        k.to_str().is_some_and(|k| {
            k.eq_ignore_ascii_case("MIMALLOC_PURGE_DELAY")
                || k.eq_ignore_ascii_case("MIMALLOC_RESET_DELAY")
        })
    })
}

fn main() -> std::process::ExitCode {
    if !purge_delay_in_env() {
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

    /// Option 15 of the mimalloc that is compiled in reads back mimalloc's default purge delay
    /// of 1000 ms (a test binary never runs `main`).
    #[test]
    fn option_15_is_the_purge_delay() {
        if super::purge_delay_in_env() {
            return;
        }
        // SAFETY: reading an option is thread safe.
        assert_eq!(
            unsafe { libmimalloc_sys::mi_option_get(super::MI_OPTION_PURGE_DELAY) },
            1000
        );
    }
}
