[<img alt="github" src="https://img.shields.io/badge/github-tolkit/tidk-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/tolkit/telomeric-identifier)
[<img alt="crates.io" src="https://img.shields.io/crates/v/tidk.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/tidk)
[<img alt="bioconda" src="https://img.shields.io/badge/bioconda-tidk-44A833?style=for-the-badge&labelColor=555555&logo=Anaconda" height="20">](https://bioconda.github.io/recipes/tidk/README.html)
[![DOI](https://zenodo.org/badge/DOI/tidk%20citation.svg)](https://doi.org/10.1093/bioinformatics/btaf049)


# A Telomere Identification toolKit (`tidk`)

`tidk` is a toolkit to identify and visualise telomeric repeats for the Darwin Tree of Life genomes. `tidk` works especially well on chromosomal genomes, but can also work on PacBio HiFi reads as well (see <a href="https://github.com/tolkit/a-telomeric-repeat-database">the telomeric repeat database</a> for many examples). There are a few modules in the tool, which may be useful to anyone investigating telomeric repeat sequences in a genome.

1. `explore` - tries to find the telomeric repeat unit in the genome.
2. `find` and `search` are essentially the same. They identify a repeat sequence in windows across the genome. `find` uses an in-built table of telomeric repeats, in `search` you supply your own.
3. `plot` does what is says on the tin, and plots the csv output of `find` or `search` as an SVG.
4. `build` builds the telomeric repeat database and saves on your local machine for use in `tidk find`.
5. `ends` calls telomeres at both ends of each sequence, for telomere-to-telomere (T2T) assessment of an assembly.
6. `length` measures telomere lengths from long reads (PacBio HiFi or ONT).

## Install

The easiest way to install is through conda:

```bash
conda install -c bioconda tidk
```

Otherwise please see the releases page here. It supports all major platforms.

If you really want to, you can complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, `cd` into it, and then run:

`cargo install --path=.`

To install into `$PATH` as `tidk`.

## Usage

Below is some usage guidance. From 0.2.3 onwards there have been breaking changes to the CLI interface. They will be pointed out below, and in the release changelog.

### Build

Before using `tidk find`, you will need to fetch the data using `tidk build`. You can do this from version 0.2.6 onwards.

### Explore 

`tidk explore` will attempt to find the simple telomeric repeat unit in the genome provided. It will report this repeat in its canonical form (e.g. TTAGG -> AACCT). A simple TSV is printed to STDOUT. Use the `distance` parameter to search only in a proportion of the chromosome arms. The default is 1% of the length of the chromosome either side, but feel free to change this. In particular with raw reads (PacBio), I'd recommend setting the distance flag to 0.5 (`--distance 0.5` or `--distance=0.5`), to process the full length of each read.

The output has a row per canonical repeat unit with the number of copies found in runs longer than `--threshold` (the `copies` column; before 0.3.0 this was named `count_repeat_runs_gt_<threshold>`), split into copies reading as the unit (e.g. `CCCTAA` for `AACCCT`) and copies reading as its reverse complement (`TTAGGG`). A real telomeric repeat reads both ways round, at both ends of a chromosome or across both strands in reads. A strong one-sided bias in reads usually comes from strand-specific basecalling errors (common for telomeres in older ONT data), and `explore` prints a warning for such units among the top five. Units that are a rotation of their own reverse complement have `NA` for both.

Since version 0.3.0, runs of a repeat may contain sequencing errors (substitutions and indels), so raw reads and noisy assemblies still give the true repeat unit. In simulations it is recovered at up to 10% error (see `benchmarks/accuracy`). Use `--exact` for the previous behaviour, where any error breaks a run. `-e`/`--error-tolerant` is still accepted and does nothing.

For example:
`tidk explore --minimum 5 --maximum 12 fastas/iyBomHort1_1.20210303.curated_primary.fa > out.tsv` searches the genome for repeats from length 5 to length 12 sequentially on the <a href="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/905/332/935/GCA_905332935.1_iyBomHort1.1/"><i>Bombus hortorum</i> genome</a>.

```
Use a range of kmer sizes to find potential telomeric repeats.
One of either length, or minimum and maximum must be specified.

Usage: tidk explore [OPTIONS] <FASTA>

Arguments:
  <FASTA>  The input fasta file

Options:
  -l, --length [<LENGTH>]        Length of substring
  -m, --minimum [<MINIMUM>]      Minimum length of substring [default: 5]
  -x, --maximum [<MAXIMUM>]      Maximum length of substring [default: 12]
  -t, --threshold [<THRESHOLD>]  Positions of repeats are only reported if they occur sequentially in a greater number than the threshold [default: 100]
      --distance [<DISTANCE>]    The distance from the end of the chromosome as a proportion of chromosome length. Must range from 0-0.5. [default: 0.01]
      --exact                    Only count runs of exactly identical repeats, as before version 0.3.0. By default runs may contain sequencing errors (substitutions and indels).
  -v, --verbose                  Print verbose output.
      --log                      Output a log file.
  -h, --help                     Print help
  -V, --version                  Print version
```

### Ends

`tidk ends` checks both ends of each sequence in an assembly for a telomere, and reports which sequences are telomere-to-telomere (T2T). Give the repeat unit with `--string` (any rotation or orientation), or leave it out to discover it from the sequence ends with the `explore` algorithm.

For each end, `tidk ends` searches a window (`--window`, default 10kb) for runs of the repeat that are at least `--min-length` bp long (default 200), allowing sequencing errors unless `--exact` is given. Each end gets one of these statuses:

- `present`: a telomere within `--max-offset` bp of the end (default 1000), on the expected strand. That is C-rich (e.g. `CCCTAA`) reading along the sequence at its start, and G-rich (`TTAGGG`) at its end.
- `wrong_strand`: a telomere at the end, but on the other strand. Often a misjoin or an inverted end.
- `not_terminal`: a telomere on the expected strand, but further than `--max-offset` from the end, so there is sequence beyond it to check.
- `absent`: none of the above.

Only `present` ends count towards T2T. Three files are written:

- `<output>.ends.tsv`: every end with its status, telomere coordinates, length, distance from the end and strand. `window_limited` is `true` when the telomere reaches the inner edge of the window, so it may be longer than reported.
- `<output>.ends.bed`: the telomeres found (0-based, half open), named by end and status.
- `<output>.ends.json`: a summary (T2T, one end, no ends, and the flagged ends) for pipelines and reports.

For example:
`tidk ends --min-sequence-length 1000000 -o ilManJurt1 ilManJurt1.1.fa.gz` calls telomeres on the chromosomes of the meadow brown (*Maniola jurtina*) assembly, skipping unplaced scaffolds.

```
Call telomeres at both ends of each sequence, for telomere-to-telomere (T2T) assessment of an assembly.

Usage: tidk ends [OPTIONS] --output <OUTPUT> <FASTA>

Arguments:
  <FASTA>  The input fasta file

Options:
  -s, --string [<STRING>]
          The telomeric repeat unit, in any rotation or orientation. Discovered from the sequence ends if not given.
  -w, --window [<WINDOW>]
          Length of sequence to search at each end, in bp [default: 10000]
      --min-length [<MIN_LENGTH>]
          Minimum telomere length to call, in bp [default: 200]
      --max-offset [<MAX_OFFSET>]
          Maximum distance of a telomere from the sequence end, in bp. Telomeres further in are flagged as not_terminal [default: 1000]
      --min-sequence-length [<MIN_SEQUENCE_LENGTH>]
          Skip sequences shorter than this, in bp (e.g. unplaced scaffolds) [default: 0]
      --exact
          Only count runs of exactly identical repeats. By default runs may contain sequencing errors (substitutions and indels).
  -o, --output <OUTPUT>
          Output filename prefix
  -d, --dir [<DIR>]
          Output directory to write files to [default: .]
  -h, --help
          Print help
  -V, --version
          Print version
```

### Length

`tidk length` measures telomere lengths from long reads, using reads that reach a chromosome end. The repeat unit must be given with `--string`; use `tidk explore` on the reads to find it if needed.

A read from a chromosome start begins with the telomere, reading C-rich (e.g. `CCCTAA`). A read from a chromosome end finishes with it, reading G-rich (`TTAGGG`). Telomeres at a read end within `--max-offset` bp (default 200) and at least `--min-length` bp long (default 200) are measured. Read ends are often degraded, so a telomere starting up to 2kb from the read end is also measured, and extended to the read end, if at least 20% of the positions in between start an exact copy of the unit. A telomere on the other strand (G-rich at a read start, or C-rich at a read end) is reported as `wrong_strand` and left out of the lengths: it can come from chimeric reads, interstitial repeats or basecalling artefacts.

A telomere is `anchored` when the read has non-telomeric sequence inward of it, and no more telomeric repeat on the same strand further in. Only anchored telomeres are used for the length summary. Otherwise the read may start or end within the telomere, or the telomere may be interrupted (e.g. by telomeric retrotransposons in some insects), so the length is a lower bound.

Two lengths are reported for each telomere:

- `telomere_length`: the whole telomere, found with error-tolerant runs. This includes most of the telomere variant repeats (TVRs) at its inner edge.
- `canonical_length`: the canonical tract at the outer edge, made of exact copies of the unit, before the TVRs. It is measured inward from the read end in windows of 20 copies of the unit, and ends where under 80% of a window's positions start an exact copy, for two windows in a row. Degraded repeat at the read end counts as part of the tract. `tvr_length` is the rest of the telomere.

On HG002 reads, `canonical_length` matches the per-read lengths from [Telogator2](https://github.com/zstephens/telogator2), which measures the canonical tract: median difference +1bp (ONT) and -40bp (HiFi). `telomere_length` is about 1kb longer, roughly the TVR length.

Lengths are summarised separately for G-rich and C-rich telomeres, as one strand can be basecalled much worse than the other (e.g. the G-rich strand in older ONT data). A warning is printed if one strand has under 10% of the telomeric reads. Two files are written:

- `<output>.length.tsv`: every telomere found at a read end, with its strand, status, `telomere_length`, `canonical_length`, `tvr_length` and whether it is anchored.
- `<output>.length.json`: for each strand, the number of reads and the median, mean, 10th and 90th percentile and maximum of anchored `telomere_length` and `canonical_length`.

Lengths have been compared with Telogator2 on human reads, but not with non-sequencing methods or other species, so treat them as estimates. The distributions can be wide and bimodal, so look at the per-read lengths rather than only the medians.

```
Measure telomere lengths from long reads (PacBio HiFi or ONT), using telomeres at read ends.

Usage: tidk length [OPTIONS] --string <STRING> --output <OUTPUT> <FASTA>

Arguments:
  <FASTA>  The input reads, as fasta

Options:
  -s, --string <STRING>            The telomeric repeat unit, in any rotation or orientation. Find it with `tidk explore` if unknown.
      --min-length [<MIN_LENGTH>]  Minimum telomere length to call, in bp [default: 200]
      --max-offset [<MAX_OFFSET>]  Maximum distance of a telomere from the read end, in bp (e.g. for untrimmed adapters) [default: 200]
      --exact                      Only count runs of exactly identical repeats. By default runs may contain sequencing errors (substitutions and indels).
  -o, --output <OUTPUT>            Output filename prefix
  -d, --dir [<DIR>]                Output directory to write files to [default: .]
  -h, --help                       Print help
  -V, --version                    Print version
```

### Find

`tidk find` will take an input clade, and match the known or putative telomeric repeat for that clade (or repeats plural) and search the genome. Now uses a custom curated telomeric repeat database. As more telomeric repeats are found and added, the dictionary of sequences used will increase.

```
Supply the name of a clade your organsim belongs to, and this submodule will find all telomeric repeat matches for that clade.

Usage: tidk find [OPTIONS] [FASTA]

Arguments:
  [FASTA]  The input fasta file

Options:
  -w, --window [<WINDOW>]  Window size to calculate telomeric repeat counts in [default: 10000]
  -c, --clade <CLADE>      The clade of organism to identify telomeres in [possible values: Crassiclitellata, Hirudinida, Phyllodocida, Eucoccidiorida, Coleoptera, Hemiptera, Hymenoptera, Lepidoptera, Odonata, Orthoptera, Plecoptera, Symphypleona, Trichoptera, Cheilostomatida, Chlamydomonadales, Accipitriformes, Anura, Aplousobranchia, Caprimulgiformes, Carangiformes, Carcharhiniformes, Carnivora, Chiroptera, Cypriniformes, Labriformes, Perciformes, Phlebobranchia, Pleuronectiformes, Rodentia, Salmoniformes, Syngnathiformes, Actiniaria, Forcipulatida, Cardiida, Pectinida, Trochida, Venerida, Heteronemertea, Apiales, Asterales, Buxales, Caryophyllales, Fabales, Fagales, Hypnales, Lamiales, Malpighiales, Myrtales, Poales, Rosales, Sapindales, Solanales]
  -o, --output <OUTPUT>    Output filename for the TSVs (without extension)
  -d, --dir <DIR>          Output directory to write files to
  -p, --print              Print a table of clades, along with their telomeric sequences
      --log                Output a log file
  -h, --help               Print help
  -V, --version            Print version
```

### Search

`tidk search` will search the genome for an input string. If you know the telomeric repeat of your sequenced organism, this will find it and return counts of occurence in windows across the genome.

```
Search the input genome with a specific telomeric repeat search string.

Usage: tidk search [OPTIONS] --string <STRING> --output <OUTPUT> --dir <DIR> <FASTA>

Arguments:
  <FASTA>  The input fasta file

Options:
  -s, --string <STRING>          The DNA string to query the genome with
  -w, --window [<WINDOW>]        Window size to calculate telomeric repeat counts in [default: 10000]
  -o, --output <OUTPUT>          Output filename for the TSVs (without extension)
  -d, --dir <DIR>                Output directory to write files to
  -e, --extension [<EXTENSION>]  The extension, defining the output type of the file [default: tsv] [possible values: tsv, bedgraph]
      --log                      Output a log file
  -h, --help                     Print help
  -V, --version                  Print version
```

### Plot

`tidk plot` will plot the output of `tidk search`.

```
SVG plot of TSV generated from tidk search.

Usage: tidk plot [OPTIONS] --tsv <TSV>

Options:
  -t, --tsv <TSV>                     The input TSV file
      --height [<HEIGHT>]             The height of subplots (px). [default: 200]
  -w, --width [<WIDTH>]               The width of plot (px) [default: 1000]
  -o, --output [<OUTPUT>]             Output filename for the SVG (without extension) [default: tidk-plot]
      --fontsize [<FONT_SIZE>]        The font size of the axis labels in the plot [default: 12]
      --strokewidth [<STROKE_WIDTH>]  The stroke width of the line graph in the plot [default: 2]
  -h, --help                          Print help
  -V, --version                       Print version
```

As an example on the ol' Square Spot Rustic <i>Xestia xanthographa</i>:

```bash
tidk find -c lepidoptera -o Xes fastas/ilXesXant1_1.20201023.curated_primary.fa

tidk plot -t finder/Xes_telomeric_repeat_windows.tsv -o ilXes -h 120 -w 800
```

## Cite

If you use this software please cite:

Max R Brown, Pablo Gonzalez de La Rosa, Mark Blaxter, tidk: a toolkit to rapidly identify telomeric repeats from genomic datasets, Bioinformatics, 2025;, btaf049, https://doi.org/10.1093/bioinformatics/btaf049


## [Cited by](https://scholar.google.com/scholar?cites=2410970761109098762&as_sdt=2005&sciodt=0,5&hl=en):

- Schell, Tilman, Carola Greve, and Lars Podsiadlowski. "Establishing genome sequencing and assembly for non-model and emerging model organisms: a brief guide." Frontiers in Zoology 22.1 (2025): 7.
- Kon, Tetsuo, et al. "Chromosome-level genome assembly of the doctor fish (Garra rufa)." Scientific Data 12.1 (2025): 1-14.
- Ryan, Kara, et al. "New Genome Assemblies for Poeciliidae: A Foundation for Adaptation Studies." Genome Biology and Evolution 17.6 (2025): evaf111.
- de Mattos, Jacqueline S., et al. "The first chromosome-scale genome assembly of a microcyclic rust, Puccinia silphii." BMC genomics 26.1 (2025): 390.
- Heaven, Thomas, et al. "Chromosome-level Assemblies of Three Candidatus Liberibacter solanacearum Vectors: Dyspersa apicalis (Förster, 1848), Dyspersa pallida (Burckhardt, 1986), and Trioza urticae (Linnaeus, 1758)(Hemiptera: Psylloidea)." Genome Biology and Evolution 17.6 (2025): evaf116.
- Kurbessoian, Tania, et al. "In host evolution of Exophiala dermatitidis in cystic fibrosis lung micro-environment." **G3**  (2023) 13(8):jkad126. doi: [10.1093/g3journal/jkad126](https://doi.org/10.1093/g3journal/jkad126)
- Yin, Denghua, et al. "Gapless genome assembly of East Asian finless porpoise." **Scientific Data** 9.1 (2022): 765.
- Leonard, Guy, et al. "A genome sequence assembly of the phototactic and optogenetic model fungus Blastocladiella emersonii reveals a diversified nucleotide-cyclase repertoire." **Genome Biology and Evolution** 14.12 (2022): evac157.
- Edwards, Richard J., et al. "A phased chromosome-level genome and full mitochondrial sequence for the dikaryotic myrtle rust pathogen, Austropuccinia psidii." **BioRxiv** (2022): 2022-04.

