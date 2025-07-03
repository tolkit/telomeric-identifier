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

## Install

The easiest way to install is through conda:

```bash
conda install -c bioconda tidk
```

Otherwise...

As with other Rust projects, you will have to complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, `cd` into it, and then run:

`cargo install --path=.`

To install into `$PATH` as `tidk`.

## Usage

Below is some usage guidance. From 0.2.3 onwards there have been breaking changes to the CLI interface. They will be pointed out below, and in the release changelog.

### Build

Before using `tidk find`, you will need to fetch the data using `tidk build`. You can do this from version 0.2.6 onwards.

### Explore 

`tidk explore` will attempt to find the simple telomeric repeat unit in the genome provided. It will report this repeat in its canonical form (e.g. TTAGG -> AACCT). Unlike previous versions, only a simple TSV is printed to STDOUT. Use the `distance` parameter to search only in a proportion of the chromosome arms. The default is 1% of the length of the chromosome either side, but feel free to change this. In particular with raw reads (PacBio), I'd recommend setting the distance flag to 0.5 (`--distance 0.5` or `--distance=0.5`), to process the full length of each read.

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
  -v, --verbose                  Print verbose output.
      --log                      Output a log file.
  -h, --help                     Print help
  -V, --version                  Print version
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
  -c, --clade <CLADE>      The clade of organism to identify telomeres in [possible values: Accipitriformes, Actiniaria, Anura, Apiales, Aplousobranchia, Asterales, Buxales, Caprimulgiformes, Carangiformes, Carcharhiniformes, Cardiida, Carnivora, Caryophyllales, Cheilostomatida, Chiroptera, Chlamydomonadales, Coleoptera, Crassiclitellata, Cypriniformes, Eucoccidiorida, Fabales, Fagales, Forcipulatida, Hemiptera, Heteronemertea, Hirudinida, Hymenoptera, Hypnales, Labriformes, Lamiales, Lepidoptera, Malpighiales, Myrtales, Odonata, Orthoptera, Pectinida, Perciformes, Phlebobranchia, Phyllodocida, Plecoptera, Pleuronectiformes, Poales, Rodentia, Rosales, Salmoniformes, Sapindales, Solanales, Symphypleona, Syngnathiformes, Trichoptera, Trochida, Venerida]
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
  -t, --tsv <TSV>          The input TSV file
      --height [<HEIGHT>]  The height of subplots (px). [default: 200]
  -w, --width [<WIDTH>]    The width of plot (px) [default: 1000]
  -o, --output [<OUTPUT>]  Output filename for the SVG (without extension) [default: tidk-plot]
  -h, --help               Print help
  -V, --version            Print version
```

As an example on the ol' Square Spot Rustic <i>Xestia xanthographa</i>:

```bash
tidk find -c lepidoptera -o Xes fastas/ilXesXant1_1.20201023.curated_primary.fa

tidk plot -t finder/Xes_telomeric_repeat_windows.tsv -o ilXes -h 120 -w 800
```

## Cite

If you use this software please cite:

Max R Brown, Pablo Gonzalez de La Rosa, Mark Blaxter, tidk: a toolkit to rapidly identify telomeric repeats from genomic datasets, Bioinformatics, 2025;, btaf049, https://doi.org/10.1093/bioinformatics/btaf049


## Cited by:
- Kurbessoian, Tania, et al. "In host evolution of Exophiala dermatitidis in cystic fibrosis lung micro-environment." **G3**  (2023) 13(8):jkad126. doi: [10.1093/g3journal/jkad126](https://doi.org/10.1093/g3journal/jkad126)
- Yin, Denghua, et al. "Gapless genome assembly of East Asian finless porpoise." **Scientific Data** 9.1 (2022): 765.
- Leonard, Guy, et al. "A genome sequence assembly of the phototactic and optogenetic model fungus Blastocladiella emersonii reveals a diversified nucleotide-cyclase repertoire." **Genome Biology and Evolution** 14.12 (2022): evac157.
- Edwards, Richard J., et al. "A phased chromosome-level genome and full mitochondrial sequence for the dikaryotic myrtle rust pathogen, Austropuccinia psidii." **BioRxiv** (2022): 2022-04.

