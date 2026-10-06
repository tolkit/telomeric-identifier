# Changelog

## 0.3.0

### Breaking changes

- `tidk explore` output has new column names and a new column order:
  `canonical_repeat_unit`, `copies`, `copies_as_unit`, `copies_as_revcomp`.
  `copies` replaces `count_repeat_runs_gt_<threshold>`. It is the number of copies of the
  unit in runs longer than `--threshold`.
- `tidk explore` now allows sequencing errors (substitutions and indels) within repeat runs by
  default. Use `--exact` for the previous behaviour, where any error breaks a run.
  `-e`/`--error-tolerant` is still accepted and does nothing.
- Counts from `tidk explore` are not comparable with earlier versions (see Fixed).

### Added

- `tidk length`: measures telomere lengths from long reads, from telomeres at read ends. It
  reports both the whole telomere, including telomere variant repeats (TVRs), and the canonical
  tract at its outer edge, which matches Telogator2's per-read lengths on HG002. Telomeres are
  extended over degraded read ends. Lengths are
  summarised separately for G-rich and C-rich telomeres, and only for reads that span the whole
  telomere.
- `tidk ends`: calls a telomere at both ends of each sequence (`present`, `wrong_strand`,
  `not_terminal` or `absent`) and reports which sequences are telomere-to-telomere, as TSV, BED
  and a JSON summary. The repeat unit is given with `--string`, or discovered from the sequence
  ends.
- `copies_as_unit` and `copies_as_revcomp` columns in `tidk explore`: copies reading as the
  canonical unit (e.g. `CCCTAA` for `AACCCT`) and as its reverse complement (`TTAGGG`). `NA` for
  units that are a rotation of their own reverse complement.
- A warning when one of the top five units has under 10% of its copies on one strand. In reads
  this usually means strand-specific basecalling errors, as seen for telomeres in older ONT data.
- `benchmarks/accuracy`: the paper's error simulation, with an optional indel error model, for
  comparing tidk versions.

### Changed

- `tidk explore` collapses rotations, reverse complements and exact multimers (e.g.
  `AAACCCTAAACCCT`) into one primitive canonical unit, instead of reporting multimers as
  separate rows.
- Runs of the same unit found at several kmer lengths in a `--minimum`/`--maximum` range are
  merged, so each telomere is counted once.
- The final aggregation step in `tidk explore` is linear instead of quadratic in the number of
  candidate repeats.
- Error-tolerant runs are extended both ways from their seed, so a repeat is measured the same
  whichever strand it reads on.
- Refactored how `tidk plot` groups records by sequence, and the removal of overlapping matches
  in `tidk search`/`tidk find`.

### Fixed

- The first run on each chromosome arm in `tidk explore` started at position 0, which inflated
  counts and let non-telomeric repeats through. In the paper's error simulations, the true
  repeat is now the top unit in 17/21 conditions with `--exact` (7/21 before) and 21/21 by
  default.
- `tidk explore` aggregation could count some runs twice and drop others.
- Runs on the right-hand chromosome arm now use whole-sequence coordinates.
