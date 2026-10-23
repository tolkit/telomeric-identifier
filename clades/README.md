# tidk database

The curated telomeric repeat database used by `tidk find`. Fetch it with `tidk build` (tidk >= 0.2.6).

`curated.csv` has one row per order and telomeric repeat. Since October 2026 it is built from the
confident calls of a tidk 0.3.0 run over one curated assembly for each of 3,422 Darwin Tree of Life
species (2,650 species, 206 orders): the repeat is at chromosome ends (`tidk ends`), and/or the
telomerase RNA template encodes it. The species named on each row is an example, and `Notes` gives
the number of species and the evidence. `Species count` is what
`tidk find --min-species` filters on (empty for older rows: one species). Species with no NCBI order are filed under their class.
Orders without a confident call in that run keep their earlier entry (tidk 0.2, 2021-2023).

Per-species calls are in [a-telomeric-repeat-database](https://github.com/tolkit/a-telomeric-repeat-database).
