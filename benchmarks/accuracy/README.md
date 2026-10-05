# explore accuracy simulation

The error simulation from the tidk paper (Supplementary Table S3). Each condition is
1000 sequences made entirely of `TTAGGG` (canonical unit `AACCCT`), with
substitution errors at 0–10%, at lengths of 600, 12000 and 30000 bp. `tidk explore` runs with
`--minimum 5 --maximum 12 --distance 0.5`, and `-t 10` (600 bp) or `-t 100` (otherwise).

```bash
cargo build --release
./benchmarks/accuracy/run.bash new=target/release/tidk > results.tsv
```

Pass several `name=binary` pairs to compare versions on the same simulated fastas.
Extra `explore` flags can follow the binary (`"exact=target/release/tidk --exact"`).
Set `ERROR_MODEL=mixed` to simulate substitutions, insertions and deletions in equal proportion
(the default, `substitution`, matches the paper).
The generator is not seeded, so counts vary slightly between runs.

## Results

`results_e0aafc0_vs_18c1e44.tsv` compares the primitive unit / run merging rewrite
(`new`, e0aafc0) with the previous version (`old`, 18c1e44).

- `AACCCT` is the top unit in 17/21 conditions for `new`, and 7/21 for `old`.
- `new` counts are exact (e.g. 100000 copies for 600 bp × 1000 with no errors).
- At 5–10% error on 12 kb and 30 kb sequences, `new` reports nothing, because no exact run
  reaches `-t 100`. `old` reported junk units in these conditions. Error-tolerant run detection (roadmap phase 2)
  is meant to fix this.

### Error-tolerant mode

`results_error_tolerant_substitution.tsv` and `results_error_tolerant_mixed.tsv` compare
exact runs with `--error-tolerant`, which became the default in 0.3.0 (`--exact` gives the old behaviour).

| | substitution | mixed (indels) |
|---|---|---|
| `AACCCT` top, exact | 16/21 | 16/21 |
| `AACCCT` top, tolerant | 21/21 | 21/21 |
| Copies recovered, tolerant, 5% error | 98–100% | 95–100% |
| Copies recovered, tolerant, 10% error | 93–100% | 71–78% |
| Units reported, tolerant (max) | 1 | 4 |
