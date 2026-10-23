# Changelog

## Unreleased

### Added

- `tidk caps`: what is at the chromosome ends, with or without a canonical telomere. Per end
  and interior window: telomere call, tandem arrays from a built-in finder grouped into
  families across the assembly, shared non-tandem sequence, assembly gaps near the tip; per
  assembly: class (`telomeric`, `terminal_satellite`, `shared_complex`, `none`), whether the
  two ends of each sequence share their top family, and segments copied from the tip between
  ends. With `--reads`, long reads past each end (open/closed), the repeats their overhangs
  carry, and read depth near the tip; reads that diverge from the assembly before the tip
  are not counted. Use the default window for the class and `-w 100000` for end pairing and
  tip copies.
  Without `--unit`, the unit discovered at the ends must have a G run (GG/CC), so AT-rich
  microsatellites are not taken as the telomere.
- `tidk find --min-species N` (`-m`): only search repeats found in at least N species of the clade
  (also applies to `--print`). The database has a new `Species count` column; rows without one
  count as one species, and older tidk versions ignore the column.

### Changed

- Clade database (`clades/curated.csv`, fetched by `tidk build`): rebuilt from the confident calls
  of the ToL-wide tidk 0.3.0 run (2,650 species, 206 orders, 301 order-repeat rows, up from 52
  orders). One row per order and repeat. Old entries kept only for orders with no new call
  (Symphypleona); old calls that disagree were dropped (e.g. Crassiclitellata AAGGAC, now AACCCT).
- `tidk find -c`: each repeat in a clade is searched once, and the clade list has no duplicates
  wherever they sit in the database.

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
