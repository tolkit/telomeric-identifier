# tidk roadmap

Where tidk goes after the 2025 Bioinformatics applications note. Phases 1–3 are tool
work; they feed the methods of the phase 4 comparative paper.

## Phase 0: housekeeping

- [x] Finish the `plot.rs` / `utils.rs` refactor (grouping in `generate_plot_data`, `remove_overlapping_indexes`)
- [x] Commit it (23a7e38)
- [ ] Fix clippy: `flush_group` has too many arguments (9/7, `plot.rs`), manual `is_multiple_of` (`utils.rs`); both still present at 0.3.0

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

## Phase 3b: chromosome caps without canonical telomeres (`tidk caps`)

Why: the paper's loss lineages (Diptera, spiders and 10 others) were checked with one-off
scripts in `tol_telomeres/shifts/d2_ends/` (`ends.py`: terminal vs interior tandem repeats
and shared sequence; `end_reads.py`: HiFi reads past the tip). Every new clade means redoing
them. D3 showed the read check is what makes a "no cap" call defensible. Spider assemblies
stop short of the chromosome ends (HiFi overhangs of 3–5 kb at 11–31 ends in 8 species).
Lepidoptera controls also look "open", because assemblies trim telomeres, but their
overhangs are the telomeric unit. So the call must rest on what the overhang is, not on
whether one exists.

One command per assembly (`tidk caps -f asm.fa [--reads hifi.fa.gz] [--unit AACCT]`); a
TSV per chromosome end, a species summary and JSON. Sections:

- [x] **Telomere call** at each end: reuse `ends.rs` (present, length, `wrong_strand`, `not_terminal`)
- [x] **End-enriched tandem repeats** (`src/caps.rs`, `src/tandem.rs`):
  - per end, the dominant tandem array in the terminal window (default 20 kb), with unit, period, array length and orientation toward the tip;
  - arrays grouped into families across ends (rotation-aware canonical k-mers);
  - enrichment: the fraction of ends where a family is top vs interior windows of the same chromosomes;
  - species class: `telomeric`, `terminal_satellite` (e.g. Diptera 160–200 bp caps), `none`.

  Built-in tandem finder (k-mer lags, periods 1 bp–2 kb, tolerant of divergent copies). Families compare consensus and repeated k-mers; with a unit, arrays ≥ 30% unit k-mers form one telomeric family (variant-rich plant telomeres). `--include`/`--exclude` regexes choose chromosomes (ToL: `'^SUPER_'`, `unloc`).
  Validation against the D2 TRF scripts on 600 species (same chromosomes): same class in 133/149 controls and 404/451 losses; telomere top at ends in 130/149 controls (scripts 132/150). Most disagreements sit at the 50% end-fraction threshold with the same period. Three interior windows per chromosome (scripts: one) catch satellites that are also interior (e.g. *Atherix*). 2–30 s and < 1 GB per assembly
- [x] **Is the end real?** (`--reads`, `src/anchor.rs`):
  - map reads to the terminal ~50 kb (minimap2 crate or a k-mer anchor; reuse the `length.rs` anchoring);
  - per end: tip reads, reads running ≥ 1 kb past the tip, median overhang, open/mixed/closed;
  - the top tandem family in the overhangs, counted per end. Pooled reads mislead: one collapsed array gave 13,220 of 14,230 *Anyphaena* overhang reads;
  - [ ] flag ends with N gaps or low coverage near the tip.

  Done with a k-mer anchor rather than minimap2, so tidk stays a pure Rust binary. 25-mers unique across the terminal 30 kb of all ends are looked up at every 8th read position, and a consistent diagonal places the read; the tip is extrapolated from the anchor nearest it. Overhangs (≥ 1 kb, up to 50 per end) go through the same tandem finder and families as the windows. On the 14 D3 species, the end states match `end_reads.py`. Lepidoptera overhangs fall in the telomeric family at 8–20 open ends (5–6 with the per-end script). Spider overhangs are 0–9% tandem, with 7–25 families across open ends. *Palloptera* overhangs are its 1,015 bp satellite. 10–13 min and < 200 MB per species (2 HiFi files, 8 threads)
- [ ] **Shared non-tandem sequence** at ends (e.g. *Drosophila* HeT-A/TART arrays): end windows aligned to each other with tandem repeats masked, vs interior; head-to-tail orientation bias.
- [ ] Later: classify end-enriched families with rep (satellite vs TRAS/SART/HeT-A-type TE vs unknown); this also covers the interrupted Lepidoptera telomeres (phase 3).
- [ ] **Validation set:** the D2/D3 species. Controls should come out `telomeric` (132/150 with the scripts); Diptera `terminal_satellite` (154/402); spiders open ends with non-telomeric overhangs (8/8); *Drosophila* shared non-tandem.
- [ ] Then rerun D2/D3 over all 12 loss lineages with `tidk caps` and cite it in the paper's methods.

## Phase 4: comparative telomere evolution paper

How the repositories fit together:

| repository | role in phase 4 |
|---|---|
| this one (tidk) | the tool; this roadmap is the project plan |
| [tolkit/a-telomeric-repeat-database](https://github.com/tolkit/a-telomeric-repeat-database) | telomere calls for every DToL species; `tol_run/` runs tidk over all assemblies |
| [ARU-life-sciences/tr_tert](https://github.com/ARU-life-sciences/tr_tert) (renamed from `tr_tert_plants`: its plant question is answered) | telomerase: the plant pipeline (top level), `tol/` (TERT and TR across the tree, built on the database's calls), `santalales/` (mistletoe request) |

Flow: tidk → `tol_run/` (calls) → `tr_tert/tol/` (TERT, TR templates against
the calls) → the paper. Both `tol_run/` and `tol/` were consolidated from the
`tol_telomeres` workspace (2026-10-10; its on-disk outputs move once the
running ToL TR round ends).

**Paper plan:** the outline anchoring the paper (thesis, sections, numbers so far,
to-dos) is `tol/notes/paper_outline.md` in tr_tert, with the literature behind it in
`tol/notes/literature_synthesis.md`.

### Telomeric repeats across the tree

- [x] tidk over all Darwin Tree of Life assemblies: 3,422 species (2026-10-08). A telomere call at ≥25% of chromosome ends in 2,292, after also checking 28 known telomeric units at the ends (the top 5 explore units alone missed telomeres ranked just below them). In `a-telomeric-repeat-database` with NCBI taxonomy
- [ ] Refresh the assembly list (one assembly has since vanished from disk: *Eloeophila maculata*) and add ERGA / EBP assemblies
- [ ] Map repeat-unit transitions onto a phylogeny: frequency, clades, rates. About 30 independent animal shifts with ≥2 species (e.g. bees/crabronids AACCCAGACCT, Hemiptera AACCATCCCT, bryozoans AAACCCC); 27 genera have a within-genus change
- [ ] Compound / mixed telomeres (phase 2 item) matter here: several calls are compounds of a simpler unit

### Telomerase: TERT

- [x] TERT in the shift lineages: found in 196/198 shift species and 105/108 controls (miniprot + Telomerase_RBD, 3 rounds of panel walking, then a profile-HMM search of six-frame ORFs). Shifts are telomerase products
- [x] Absent from Diptera (0/6, expected) and spiders (0/6, which also have no telomere at any end): telomerase loss in spiders
- [ ] TERT sequence changes associated with repeat shifts (TRBD, motif T) across the 196 shift and 105 control proteins
- [ ] **Santalales / mistletoe (collaborator request):** locate TERT (and TR, POT1, TRB, CST, dyskerin, Ku) in non-*Viscum* Santalales - *Santalum*, *Malania* (RefSeq-annotated), Loranthaceae mistletoes *Taxillus*, *Scurrula*, *Dendrophthoe* (an independent origin of aerial parasitism), *Thesium*, *Balanophora*, ... - then their synteny into *Viscum* (94 Gb; no TERT or TR found genome-wide, no telomere at any end) and a targeted search of that interval for a gene or pseudogene. Set up in `tr_tert/santalales/`; download paused

### Telomerase: TR

- [x] TR templates validate the telomere calls: where a TR is found (Hymenoptera models from Fajkus et al. 2023 Table S1, Lepidoptera from their Table S4, Rfam vertebrate CM), its template encodes the tidk repeat in 1,087/1,096 strongly called species; most exceptions are tidk miscalls
- [x] Hymenoptera: model walking (6 rounds, telomere-validated seeds), per-lineage models, covariance models, synteny and copy-aware checks: strongly called species with a TR 109 → 157/188. Lessons: lineage-specific models beat pooled ones; genomes can hold several TR-like copies, only one with a template
- [ ] Hymenoptera still without a TR: Argidae (AACCC; every method fails, including template conservation), *Xiphydria*, *Perithous albicinctus*
- [ ] Template co-variation as a de novo method: works within a genus (*Bombus*: the published TR blind, rank 1 of 2), not beyond; needs distinctive unit changes
- [ ] **ToL-wide TR search** (`tr_tert/tol/tr_tol`): stage 1, every model on every genome (52 HMMs + 3 Rfam CMs): 1,572 validated, all in clades with existing models, plus one beetle (*Otiorhynchus rugosostriatus*, via a Lepidoptera model, E 5.7e-42: check). Stage 2 (walking per lineage, 169 seed models) running. Then: the 7 large or failed genomes, synteny where TRs stay put, co-variation in the 27 genera with a repeat change. Groups with no seed (Hemiptera, Odonata, molluscs, fungi, mosses) need literature TRs or de novo discovery
- [ ] **Plants:** pair repeat units with the TR/TERT data in `tr_tert` to test whether TR template changes explain repeat transitions

### TR genomic position: a neat result

- [x] **TR genes move in Hymenoptera but not in Lepidoptera.** Survey of 250 random pairs (TR flanks projected into a relative by nucleotide or BUSCO-gene anchors): TR position conserved in Lepidoptera at genus 25/25, family 22/25, superfamily 14/21; in Hymenoptera at genus 16/25, family 4/24, superfamily 8/23. Between pimpline genera, crabronid subfamilies, *Tenthredo*/*Euura* and *Cephus*/*Urocerus*, the orthologous interval is conserved and empty (no relic even at E ≤ 10), with the TR on another chromosome
- [x] Prior work: Klapproth et al. 2025 (bioRxiv 10.1101/2025.04.29.651168, Stadler/Findeiß, with J. Chen) found 8 TR copies in *Andrena dorsata* with two functional template types making two telomeric repeats (5' RACE in *A. hattorfiana*), TR paralogs diverging within *Andrena*, and single-copy TERT; TR paralogs are also known in plants (Závodník et al. 2023). So multi-copy TRs in *Andrena* are known; what is new here is the scale (all Hymenoptera vs Lepidoptera), the positional survey, the pimpline one-functional-copy pattern, and testing their "duplication drives repeat diversity" hypothesis tree-wide (do repeat shifts track multi-copy lineages?). It also explains our odd *Andrena* calls: two functional TRs, two repeats
- [x] Hymenoptera genomes carry more TR-like copies (mean 3.3 vs 1.2; ≥3 copies 29% vs 3%, p = 3e-22). Pimplines: 3–7 copies, usually one with a template; halictids up to 52 (*Lasioglossum calceatum*)
- [ ] Mechanism: no retrocopy signature (poly-A, target-site duplications at background). Two patterns: dispersed TR-sized insertions with unique flanks (pimplines, *Lasioglossum*) and copies inside larger shared repeats (*Andrena*, *Seladonia*). Next: full-length copy alignments (5' ends), TE annotation of flanks, copy ages; when the functional copy switches, and whether switches coincide with repeat shifts

### Other

- [ ] Interstitial telomeric sequences as markers of chromosome fusion and karyotype evolution
- [ ] Survey taxa that lack canonical telomeres: Diptera (satellite caps in 154/402), spiders (assemblies open at the ends; caps unresolved), 10 other losses, *Viscum*. Moving to `tidk caps` (phase 3b)
- [ ] Target journal: MBE, Genome Research or PNAS, with tidk v1.0 released alongside
