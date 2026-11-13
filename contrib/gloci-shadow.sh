#!/usr/bin/env bash
# Run a pipeline step with the tool it uses today and with gloci, and record whether the two
# outputs are the same. The step keeps the old tool's output and exit status, so nothing
# downstream changes while the evidence is collected.
#
#   gloci-shadow.sh --name merge-peaks --log shadow.tsv \
#       --ref 'bedtools merge -i peaks.bed' \
#       --new 'gloci merge peaks.bed --chrom-order lex' > merged.bed
#
# Both commands run under `bash -o pipefail` and must write their result to standard output.
# One tab-separated line is appended to the log:
#
#   date  name  verdict  ref_seconds  new_seconds  ref_lines  new_lines  ref_status  new_status
#
# where verdict is `same`, `differs`, `ref-failed` or `new-failed`. For a run that differs, the
# first differing lines are kept in <log>.d/<name>.<date>.diff.
set -u

usage() {
    sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
    echo
    echo "Options:"
    echo "  --name NAME    a label for the log (required)"
    echo "  --ref CMD      the command in use today (required)"
    echo "  --new CMD      the gloci command that should give the same output (required)"
    echo "  --log FILE     the log to append to (default: gloci-shadow.tsv)"
    echo "  --unordered    compare the outputs as sorted lines (row order may differ)"
    echo "  --keep N       differing lines to keep per run (default: 20)"
}

name="" ref="" new="" log="gloci-shadow.tsv" unordered=0 keep=20
while [ $# -gt 0 ]; do
    case "$1" in
        --name) name="${2-}"; shift 2 ;;
        --ref) ref="${2-}"; shift 2 ;;
        --new) new="${2-}"; shift 2 ;;
        --log) log="${2-}"; shift 2 ;;
        --keep) keep="${2-}"; shift 2 ;;
        --unordered) unordered=1; shift ;;
        -h | --help) usage; exit 0 ;;
        *) echo "gloci-shadow: unknown option $1 (see --help)" >&2; exit 2 ;;
    esac
done
if [ -z "$name" ] || [ -z "$ref" ] || [ -z "$new" ]; then
    echo "gloci-shadow: --name, --ref and --new are required (see --help)" >&2
    exit 2
fi
case "$name" in
    */* | *$'\t'*) echo "gloci-shadow: --name may not contain a slash or a tab" >&2; exit 2 ;;
esac

work="$(mktemp -d "${TMPDIR:-/tmp}/gloci-shadow.XXXXXX")" || exit 1
trap 'rm -rf "$work"' EXIT

# Seconds as bash's own `time` reports them, so no GNU time or date is needed.
run() { # run OUT ERR TIMEFILE CMD
    local out="$1" err="$2" tf="$3" cmd="$4" status
    TIMEFORMAT=%R
    { time bash -o pipefail -c "$cmd" > "$out" 2> "$err"; } 2> "$tf"
    status=$?
    return $status
}

run "$work/ref.out" "$work/ref.err" "$work/ref.time" "$ref"
ref_status=$?
run "$work/new.out" "$work/new.err" "$work/new.time" "$new"
new_status=$?

# The step's own result: the old tool's output, messages and status, untouched.
cat "$work/ref.out"
cat "$work/ref.err" >&2

if [ "$unordered" -eq 1 ]; then
    LC_ALL=C sort "$work/ref.out" > "$work/ref.cmp"
    LC_ALL=C sort "$work/new.out" > "$work/new.cmp"
else
    ln "$work/ref.out" "$work/ref.cmp"
    ln "$work/new.out" "$work/new.cmp"
fi

stamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
if [ "$ref_status" -ne 0 ]; then
    verdict="ref-failed"
elif [ "$new_status" -ne 0 ]; then
    verdict="new-failed"
elif cmp -s "$work/ref.cmp" "$work/new.cmp"; then
    verdict="same"
else
    verdict="differs"
fi

if [ "$verdict" = differs ] || [ "$verdict" = new-failed ]; then
    mkdir -p "$log.d"
    {
        echo "# $stamp $name: $verdict"
        echo "# ref: $ref"
        echo "# new: $new"
        if [ "$verdict" = new-failed ]; then
            echo "# gloci exited $new_status and wrote:"
            head -n "$keep" "$work/new.err"
        else
            diff "$work/ref.cmp" "$work/new.cmp" | head -n "$keep"
        fi
    } > "$log.d/$name.${stamp//:/}.diff"
fi

lines() { wc -l < "$1" | tr -d ' '; }
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$stamp" "$name" "$verdict" \
    "$(tail -n 1 "$work/ref.time")" "$(tail -n 1 "$work/new.time")" \
    "$(lines "$work/ref.out")" "$(lines "$work/new.out")" \
    "$ref_status" "$new_status" >> "$log"

exit "$ref_status"
