# Reference tools for differential tests and benchmarks (pinned).
FROM mambaorg/micromamba:1.5.10
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        bedtools=2.31.1 bedops=2.4.41 vcftools=0.1.16 htslib=1.21 bcftools=1.21 \
        hyperfine coreutils \
    && micromamba clean --all --yes
