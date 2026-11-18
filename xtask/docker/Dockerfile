# Reference tools for differential tests and benchmarks (pinned).
FROM mambaorg/micromamba:1.5.10
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        bedtools=2.31.1 bedops=2.4.41 vcftools=0.1.16 htslib=1.21 bcftools=1.21 samtools=1.21 \
        mosdepth=0.3.10 \
        python=3.12 pysam=0.23.3 bx-python=0.16.0 pybigwig=0.3.29 numpy=2.5.5 pip \
        hyperfine=2.0.0 coreutils \
    && micromamba clean --all --yes
# CrossMap 0.7.4 (GPL-3.0; test-only, never shipped) from PyPI without its dependencies: every
# bioconda crossmap build depends on a UCSC wig-to-bigWig tool package, which this image must
# not contain.
# `CrossMap --version` prints 0.7.3 for 0.7.4 (its version string was not bumped).
RUN printf 'crossmap==0.7.4 --hash=sha256:2c976861eb0d64e3afd8cc03413ffe98c496e8fd4e9d60c7913593cff4e62648\n' > /tmp/crossmap.txt \
    && micromamba run -n base python -m pip install --no-cache-dir --no-deps --require-hashes -r /tmp/crossmap.txt \
    && rm /tmp/crossmap.txt \
    && ! micromamba list -n base | grep -qi ucsc
