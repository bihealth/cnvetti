# CNVetti CHANGELOG

## HEAD (v0.1.0; unreleased)

- Adding `cnvetti visualize` sub commands:
    - `cov-to-igv` (extract coverage from BCF to IGV format).
- Adding `cnvetti quick` sub commands:
    - `wis-build-model` (build model for within-sample CNV calling from BAM files)
    - `wis-call` (within-sample CNV calling from BAM and model BCF file)
- Adding `cnvetti cmd` sub commands:
    - `coverage` (genome-wide fragments/coverage, or target-wise fragments)
    - `normalize` (normalize with total coverage or by GC-wise median)
    - `merge-cov` (merge coverage BCF files)
    - `build-model-wis` (build model for the within-sample CNV calling approach)
    - `mod-cov` (compute model-based coverage)
- Starting out with CLI skelleton for `cnvetti cmd *`.
