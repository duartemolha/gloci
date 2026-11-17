#!/bin/sh
# Rebuilds the committed mosdepth 0.3.10 outputs under this directory, one directory per case.
# Run it from this directory, with docker; the image is the comparison image (xtask/docker):
#   sh regen.sh                    # uses gloci-diff:p0
#   IMAGE=gloci-diff:other sh regen.sh
# The cases need mosdepth 0.3.10 in the image: rebuild it after the Dockerfile change
# (docker build -t gloci-diff:p0 -f xtask/docker/Dockerfile xtask/docker), or point IMAGE at
# an image that has it. The `gloci/` directories (pinned files of documented differences) are kept.
# The case table is shared with depth_files.rs (CASES there): keep the two in step. Each case
# holds mosdepth's files without the prefix, with `regions.bed` and `thresholds.bed` decompressed.
set -eu
cd "$(dirname "$0")"
IMAGE="${IMAGE:-gloci-diff:p0}"
IN=../../../../xtask/diff/inputs/samtools

# Check the image before anything is deleted.
if ! docker run --rm "$IMAGE" mosdepth --version >/dev/null 2>&1; then
  echo "regen.sh: the image $IMAGE has no mosdepth: rebuild it with docker build -t gloci-diff:p0 -f xtask/docker/Dockerfile xtask/docker" >&2
  exit 1
fi

case_run() {
  name="$1"
  input="$2"
  shift 2
  mkdir -p "$name"
  find "$name" -maxdepth 1 -type f -delete
  # `-n`: no per-base file. Usage: case_run NAME INPUT MOSDEPTH-ARGS...
  docker run --rm -u "$(id -u):$(id -g)" -v "$PWD:/work" -v "$(cd "$IN" && pwd):/in:ro" \
    -w /work "$IMAGE" sh -c '
      set -eu
      name="$1"; input="$2"; shift 2
      mosdepth -n "$@" "$name/o" "$input"
      cd "$name"
      for f in o.mosdepth.*; do mv "$f" "${f#o.}"; done
      for f in o.*.bed.gz; do [ -e "$f" ] || continue; b="${f#o.}"; zcat "$f" > "${b%.gz}"; rm "$f"; done
      rm -f o.*.csi
    ' sh "$name" "$input" "$@"
}

case_run depth-window-100 /in/depth.bam --by 100
case_run depth-targets-thresholds /in/depth.bam --by /in/targets.bed -T 1,2,10
case_run pairs-window-thresholds /in/mosdepth.bam --by 100 -T 1,2,3
case_run pairs-bed-mapq /in/mosdepth.bam -Q 20 --by /in/mosdepth.bed -T 2,1
case_run pairs-no-by /in/mosdepth.bam
case_run pairs-flags-zero /in/mosdepth.bam -F 0 --by 500
case_run means-base-by-base /in/mosdepth-means.bam --by 200
case_run empty /in/empty.bam --by 100
case_run mates-contained /in/mates.bam --by 100
case_run stale-pair-entry /in/mosdepth-stale.bam -F 1540 --by 100
case_run bed-blank-line /in/mosdepth.bam --by /in/mosdepth-blank.bed
