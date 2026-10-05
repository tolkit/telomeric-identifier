# tidk roadmap

Where tidk goes after the 2025 Bioinformatics applications note. Phases 1–3 are tool
work; they feed the methods of the phase 4 comparative paper.

## Phase 0: housekeeping

- [x] Finish the `plot.rs` / `utils.rs` refactor (grouping in `generate_plot_data`, `remove_overlapping_indexes`)
- [ ] Commit it
- [ ] Fix clippy: `flush_group` has too many arguments (`plot.rs`), `is_multiple_of` (`utils.rs`)

## Phase 1: `explore` output quality

- [x] Collapse rotations, reverse complements and exact multimers to one primitive canonical unit (`utils::primitive_telomere_unit`)
- [x] Replace the O(m²) pairwise aggregation with a hashmap keyed on the primitive unit
- [x] Merge overlapping runs per sequence, so a telomere found at several k in a `--minimum/--maximum` range is counted once
- [x] Fix: the first run on each arm started at position 0, which inflated counts and let junk repeats through
- [x] Fix: right-arm runs now use whole-record coordinates
- [ ] Rename the count column (it is now copies of the unit in merged runs, not "runs > threshold"), with a note in the changelog
- [ ] Optionally report the run length distribution and which k found each unit
- [x] Re-run the paper's error simulations (0–10% substitution error) against the previous version: the true unit is top in 17/21 conditions (was 7/21). The other 4 (12 kb/30 kb at 5–10%) report nothing, because no exact run reaches the threshold; the old version reported junk there

## Phase 2: error tolerance and compound repeats

- [x] Error-tolerant run detection (`explore --error-tolerant`): runs are seeded on exact adjacent chunks, then extended with a scoring rule that tolerates substitutions, with a penalised rotation switch for indels and consensus labelling. True unit is top in 21/21 conditions for both substitution-only and mixed (indel) errors, with ≥95% of copies recovered at 5% error
- [x] Real-data test (oak dhQueRobu3.1, meadow brown ilManJurt1.1; assemblies plus ~3 Gbp HiFi/ONT subsamples): same top unit as exact on assemblies, at the same speed; 3–50× more telomeric copies from reads; old R9 ONT gives nothing exact, but the tolerant top unit is a basecalling artefact (below)
- [ ] Decide whether `--error-tolerant` becomes the default (cost: more secondary units from imperfect microsatellites, e.g. ACAG/AAAG)
- [x] Report strand balance for each unit (`count_as_unit`/`count_as_revcomp` columns, plus a warning below 10% on the minor strand). On real data, assemblies and HiFi are balanced (32–47% minor strand); every unit in the 2019 meadow brown ONT is 0%. That covers satellites too (ACAG: long arrays only basecalled cleanly on one strand), not just telomeres. Both ONT sets miscall the G-rich telomere strand (meadow brown TTAGG→TTGGG, oak TTTAGGG→AGGG), while the C-rich strand is called correctly. True telomeric units should appear in both orientations across reads, so a strong strand bias would flag basecalling artefacts
- [ ] Remaining tolerant mode artefacts: small 11-mer units (≤0.6%) at 10% indel error, where both seed copies share an indel
- [ ] Describe compound / HOR telomeres (e.g. *Bombus* AACCT + AACCCG, mixed plant TTTAGGG/TTAGGG): unit composition and alternation per run, not a single winner
- [x] Benchmark on simulated reads (`benchmarks/accuracy`, `ERROR_MODEL=mixed` adds indels)
- [ ] Benchmark on known-repeat genomes and real reads

## Phase 3: telomere QC for T2T assemblies

- [ ] New subcommand (working name `tidk ends`): per-chromosome-end presence, length, strand and T2T status
- [ ] BED + JSON output, plus a whole-assembly summary ("18/20 T2T")
- [ ] MultiQC module; integrate with Tree of Life assembly/curation pipelines
- [ ] Per-read telomere length from HiFi/ONT reads, for any species (no known repeat needed)

## Phase 4: comparative telomere evolution paper

- [ ] Run tidk over all chromosome-level Tree of Life / ERGA / EBP assemblies (thousands, not the original 500)
- [ ] Map repeat-unit transitions onto a phylogeny: frequency, clades, rates
- [ ] **Plants:** pair repeat units with existing telomerase RNA (TR) and TERT data to test whether TR template changes explain repeat transitions
- [ ] Interstitial telomeric sequences as markers of chromosome fusion and karyotype evolution
- [ ] Survey taxa that lack canonical telomeres
- [ ] Target journal: MBE, Genome Research or PNAS, with tidk v1.0 released alongside
