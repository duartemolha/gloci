# Stability

gloci is young. Its commands are in two tiers, so that you know what you can build a pipeline on today.

## Stable core

`genomes`, `index`, `sort`, `view`, `intersect`, `merge`, `subtract`, `closest` and `coverage`.

`gloci --help` lists them first, under the heading `Stable core`.

For these commands:

- options and output change only with a version bump and an entry in the [changelog](changelog.md)
  (before 1.0, a breaking change bumps the minor version);
- every documented behaviour has a golden case against the tool it replaces, checked on every commit;
- the nightly workflow, which runs weekly, covers them with multi-GB inputs, memory ceilings and a
  live comparison with the pinned reference tools.

## Experimental

Every other command. Its help starts with a line that says so:

```text
$ gloci cluster --help
Experimental: checked against the tool it replaces, but its options and output may still change
between releases.
```

`gloci --help` puts every other command in a group whose heading ends in `(experimental)`.

Experimental does not mean untested: these commands have golden cases too, and most are compared with
bedtools, bcftools or samtools output byte for byte. It means that they have had less use on real data,
so an option may be renamed or an output column changed in a later release without a deprecation
period. Compare the output with the tool you are replacing before you rely on one, as described below.

## Running gloci beside the tool it replaces

Before moving a pipeline step to gloci, run both tools on your own data for a while and compare.
`contrib/gloci-shadow.sh` does that for one step: it runs your current command and the gloci command,
passes your current command's output through unchanged, and appends one line to a log saying whether
the two outputs were identical, with the time each took.

```sh
contrib/gloci-shadow.sh --name merge-peaks --log shadow.tsv \
    --ref 'bedtools merge -i peaks.bed' \
    --new 'gloci merge peaks.bed' > merged.bed
```

The step keeps producing the old tool's output, so nothing downstream changes while you collect
evidence. Both outputs are held in `$TMPDIR` until they are compared, so it needs room for them. A run that differs keeps the first differing lines next to the log. Before reporting a
difference, look it up in [Differences](differences.md): the intentional ones are listed there with the
option that gives the old behaviour (`--vcf-span ref`, `--inclusive-points`).

`--unordered` compares the outputs as sets of lines, for commands whose row order is allowed to differ.

## Promoting a command

A command moves to the stable core when:

1. it has been run beside the tool it replaces on real data, and every difference found is either fixed
   or listed in [Differences](differences.md);
2. it has a large-input test with a memory ceiling in the nightly workflow;
3. its options have been reviewed once more for names that should change, because after promotion they
   cannot change freely.

The list lives in one place, `CORE_COMMANDS` in `crates/gloci/src/lib.rs`; a test fails when a command
is neither in it nor marked experimental.
