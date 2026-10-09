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
- [x] Rename the count column to `copies` (0.3.0)
- [ ] Optionally report the run length distribution and which k found each unit
- [x] Re-run the paper's error simulations (0–10% substitution error) against the previous version: the true unit is top in 17/21 conditions (was 7/21). The other 4 (12 kb/30 kb at 5–10%) report nothing, because no exact run reaches the threshold; the old version reported junk there

## Phase 2: error tolerance and compound repeats

- [x] Error-tolerant run detection (`explore --error-tolerant`): runs are seeded on exact adjacent chunks, then extended with a scoring rule that tolerates substitutions, with a penalised rotation switch for indels and consensus labelling. True unit is top in 21/21 conditions for both substitution-only and mixed (indel) errors, with ≥95% of copies recovered at 5% error
- [x] Real-data test (oak dhQueRobu3.1, meadow brown ilManJurt1.1; assemblies plus ~3 Gbp HiFi/ONT subsamples): same top unit as exact on assemblies, at the same speed; 3–50× more telomeric copies from reads; old R9 ONT gives nothing exact, but the tolerant top unit is a basecalling artefact (below)
- [x] Make error-tolerant runs the default (0.3.0), with `--exact` for the old behaviour; `-e` is kept as a hidden no-op. Cost: more secondary units from imperfect microsatellites, e.g. ACAG/AAAG
- [x] Report strand balance for each unit (`copies_as_unit`/`copies_as_revcomp` columns, plus a warning below 10% on the minor strand). On real data, assemblies and HiFi are balanced (32–47% minor strand); every unit in the 2019 meadow brown ONT is 0%. That covers satellites too (ACAG: long arrays only basecalled cleanly on one strand), not just telomeres. Both ONT sets miscall the G-rich telomere strand (meadow brown TTAGG→TTGGG, oak TTTAGGG→AGGG), while the C-rich strand is called correctly. True telomeric units should appear in both orientations across reads, so a strong strand bias would flag basecalling artefacts
- [ ] Remaining tolerant mode artefacts: small 11-mer units (≤0.6%) at 10% indel error, where both seed copies share an indel
- [ ] Describe compound / HOR telomeres (e.g. *Bombus* AACCT + AACCCG, mixed plant TTTAGGG/TTAGGG): unit composition and alternation per run, not a single winner
- [x] Benchmark on simulated reads (`benchmarks/accuracy`, `ERROR_MODEL=mixed` adds indels)
- [ ] Benchmark on known-repeat genomes and real reads

## Phase 3: telomere QC for T2T assemblies

- [x] `tidk ends`: presence, length, strand and T2T status at each chromosome end, with `wrong_strand` and `not_terminal` flags. Oak dhQueRobu3.1 is 7/12 T2T. Meadow brown ilManJurt1.1 is 12/30, with 10 `not_terminal` ends: correctly oriented telomeres 1.5–8 kb in from the end (possibly TRAS/SART retrotransposon insertions; not yet checked)
- [x] BED + JSON output, plus a whole-assembly summary ("18/20 T2T")
- [ ] MultiQC module; integrate with Tree of Life assembly/curation pipelines
- [x] `tidk length`: per-read telomere length from HiFi/ONT reads, summarised by strand and using only anchored reads. Oak HiFi gives median 4.3 kb (G-rich) and 6.7 kb (C-rich) from about 4× coverage; the 2019 meadow brown ONT gives almost no anchored telomeres
- [x] Compare `tidk length` with Telogator2 on HG002 (HiFi and ONT telomere reads). tidk measures 93–98% of Telogator2's reads, and per-read whole-telomere lengths correlate with r = 0.91. Telogator2 measures the canonical tract only, so tidk now reports `canonical_length` (median difference +1 / −40 bp, r = 0.78–0.81) as well as the whole telomere including TVRs (`telomere_length`, about +1 kb)
- [x] Degraded read ends: telomeres starting up to 2 kb from the read end are extended to it when at least 20% of the positions in between start an exact copy of the unit. HG002 HiFi reads measured out of Telogator2's went from 93.4% to 95.6%, with agreement unchanged
- [ ] The remaining 4.4% of HiFi reads Telogator2 measures have read ends of heavily degenerate repeat (e.g. `GGGTTGGGG`, `CCCTACCCC`) with almost no exact copies. A G/C composition test catches about two thirds of these, but also 11% of subtelomere next to telomeres, so it isn't used. A variant-kmer set like Telogator2's might separate them better
- [ ] Validate against a non-sequencing method (e.g. TRF Southern blots) and in another species
- [ ] Understand why C-rich telomeres are longer than G-rich in the same read sets (meadow brown HiFi p = 0.01; oak HiFi and ONT in the same direction). Reverse complementing the reads swaps the result, so it is in the data, not the algorithm
- [ ] Interrupted telomeres (e.g. Lepidoptera TRAS/SART insertions): report the whole telomeric region as well as the uninterrupted terminal tract

## Phase 4: comparative telomere evolution paper

- [x] Run tidk over all Darwin Tree of Life assemblies: 3,422 species (2026-10-08; `/lustre/scratch122/tol/teams/blaxter/users/mb39/tol_telomeres`). A telomere call at ≥25% of chromosome ends in 2,292 species, after also checking known telomeric units at the ends (the top 5 explore units alone missed telomeres ranked just below them). Merged into `a-telomeric-repeat-database` (branch `tidk-0.3`) with NCBI taxonomy. ERGA / EBP still to add
- [ ] Map repeat-unit transitions onto a phylogeny: frequency, clades, rates. About 30 independent animal shifts with ≥2 species (e.g. bees/crabronids AACCCAGACCT, Hemiptera AACCATCCCT, bryozoans AAACCCC)
- [ ] **Plants:** pair repeat units with existing telomerase RNA (TR) and TERT data to test whether TR template changes explain repeat transitions
- [x] **TERT in the shift lineages** (`tol_telomeres/shifts`): found in 196/198 shift species and 105/108 controls, so shifts are telomerase products; absent from Diptera (0/6, expected) and spiders (0/6, with no telomere at any end: telomerase loss)
- [x] **TR templates validate the telomere calls**: where a TR is found (Hymenoptera models from Fajkus et al. 2023, Lepidoptera models from their Table S4, Rfam vertebrate CM), its template encodes the tidk repeat in 1,087/1,096 strongly called species; most exceptions are tidk miscalls
- [ ] TR in the remaining shift lineages: model walking (Hymenoptera 109 → 146/188), covariance models, synteny and template co-variation (within *Bombus* it finds the published TR blind, rank 1 of 2)
- [ ] **Neat result: TR genes move in Hymenoptera but not in Lepidoptera.** Projecting a known TR's flanks (nucleotide or BUSCO-gene anchors) into a relative finds the relative's TR in place across Lepidoptera from genus to superfamily (7/9 pairs) and across core sawfly families, but between pimpline genera, crabronid subfamilies, *Tenthredo*/*Euura* and *Cephus*/*Urocerus* the orthologous interval is conserved and empty, with the TR on another chromosome. TR relocation (retroposition?) looks frequent in Hymenoptera. Survey of 250 random pairs (`tol_telomeres/shifts/denovo`): TR position conserved in Lepidoptera at genus 25/25, family 22/25, superfamily 14/21; in Hymenoptera at genus 16/25, family 4/24, superfamily 8/23. Hymenoptera genomes carry more TR-like copies (mean 3.3 vs 1.2, ≥3 copies 29% vs 3%, p = 3e-22), and pimplines have 3–7 copies with usually one template-bearing, the ancestral position empty rather than decayed. To do: the mechanism (are the copies retrocopies: poly-A tails, target-site duplications, flanking TEs?), when the switch in the functional copy happens, and whether it coincides with telomeric-repeat shifts
- [ ] Interstitial telomeric sequences as markers of chromosome fusion and karyotype evolution
- [ ] Survey taxa that lack canonical telomeres (started: Diptera, spiders, mayflies)
- [ ] Target journal: MBE, Genome Research or PNAS, with tidk v1.0 released alongside
