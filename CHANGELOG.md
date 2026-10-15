# nf-core/tumourevo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [date]

### `Changed`

- `SUBCLONAL_INTERPRETATION`: revised expansion score. The evidence score Φ (driver, tail, signature) is multiplied by a reliability factor R, the minimum of read support, cluster size, own density mode and posterior predictive fit (purity and copy number are left to CNAqc). Tiers are assigned on the final score; the clonal cluster is re-anchored on CCF; φ_driver is graded by tail status and corrected for chance; φ_signature is read against a same-size null; non-estimable components take a neutral value instead of being dropped.
- `GENOME_INTERPRETER` passes the clustering fits and the CNAqc joint table to `SUBCLONAL_INTERPRETATION`; `PLOT_CLONE_TREE` uses the tiers computed there.

### `Added`

- `--interpreter_*` parameters controlling the expansion score.
- `*_scores.tsv` and `*_sample_quality.tsv` outputs of `SUBCLONAL_INTERPRETATION`, and nf-test stub tests for the module.

## v1.0.0 - [date]

Initial release of nf-core/tumourevo, created with the [nf-core](https://nf-co.re/) template.
