# A Telomere Identification toolKit (`tidk`)

`tidk` is a toolkit to identify and visualise telomeric repeats for the Darwin Tree of Life genomes. `tidk` works especially well on chromosomal genomes, but can also work on PacBio HiFi reads well (see <a href="https://github.com/tolkit/a-telomeric-repeat-database">the telomeric repeat database</a> for many examples). There are a few modules in the tool, which may be useful to anyone investigating telomeric repeat sequences in a genome.

1. `explore` - tries to find the telomeric repeat unit in the genome.
2. `find` and `search` are essentially the same. They identify a repeat sequence in windows across the genome. `find` uses an in-built table of telomeric repeats, in `search` you supply your own.
3. `plot` does what is says on the tin, and plots the csv output of `find` or `search` as an SVG.

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

### Explore 

`tidk explore` will attempt to find the simple telomeric repeat unit in the genome provided. It will report this repeat in its canonical form (e.g. TTAGG -> AACCT). Unlike previous versions, a simple TSV is printed to STDOUT. Use the `distance` parameter to search only in a proportion of the chromosome arms.

For example:
`tidk explore --minimum 5 --maximum 12 fastas/iyBomHort1_1.20210303.curated_primary.fa` searches the genome for repeats from length 5 to length 12 sequentially on the <a href="https://www.ebi.ac.uk/ena/browser/view/PRJEB43539"><i>Bombus hortorum</i> genome</a>.

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
      --distance [<DISTANCE>]    The distance from the end of the chromosome as a proportion of chromosome length. [default: 0.1]
  -v, --verbose                  Print verbose output.
      --log                      Output a log file.
  -h, --help                     Print help
  -V, --version                  Print version
```

### Find

`tidk find` will take an input clade, and match the known telomeric repeat for that clade (or repeats plural) and search the genome. Uses the <a href="http://telomerase.asu.edu/sequences_telomere.html">telomeric repeat database</a>. As more telomeric repeats are found and added, the dictionary of sequences used will increase. We have a lot more clades of late, but do sanity check the repeats as the database is not yet curated. I'm actively working on a curated database.

```
Supply the name of a clade your organsim belongs to, and this submodule will find all telomeric repeat matches for that clade.

Usage: tidk find [OPTIONS] [FASTA]

Arguments:
  [FASTA]  The input fasta file

Options:
  -w, --window [<WINDOW>]  Window size to calculate telomeric repeat counts in [default: 10000]
  -c, --clade <CLADE>      The clade of organism to identify telomeres in [possible values: Accipitriformes, Actiniaria, Agaricales, Alismatales, Amphilepidida, Anura, Apiales, Aplousobranchia, Aquifoliales, Araneae, Artiodactyla, Asparagales, Asterales, Atheriniformes, Balanomorpha, Boraginales, Brassicales, Buxales, Camarodonta, Caprimulgiformes, Carcharhiniformes, Cardiida, Carnivora, Caryophyllales, Celastrales, Chaetocerotales, Cheilostomatida, Chiroptera, Chitonida, Chlamydomonadales, Coleoptera, Comatulida, Crassiclitellata, Cucurbitales, Cypriniformes, Decapoda, Dioctophymatida, Dipsacales, Ericales, Eucoccidiorida, Euglenales, Eulipotyphla, Fabales, Fagales, Forcipulatida, Fucales, Gentianales, Geophilomorpha, Geraniales, Gigartinales, Glomerida, Hemiptera, Heteronemertea, Hirudinida, Hymenoptera, Hypnales, Isochrysidales, Isopoda, Lamiales, Lepidoptera, Liliales, Lithobiomorpha, Littorinimorpha, Lunulariales, Lycopodiales, Malpighiales, Malvales, Megaloptera, Myrtales, Neuroptera, Nudibranchia, Odonata, Opiliones, Orthoptera, Ostreida, Palmariales, Pectinida, Pelecaniformes, Perciformes, Phlebobranchia, Phyllodocida, Plecoptera, Poales, Polytrichales, Primates, Procellariiformes, Pyrenomonadales, Ranunculales, Raphidioptera, Rhabditida, Rodentia, Rosales, Sabellida, Salmoniformes, Sapindales, Scombriformes, Scorpiones, Solanales, Sphagnales, Stolidobranchia, Symphypleona, Trichoptera, Trochida, Venerida]
  -o, --output <OUTPUT>    Output filename for the TSVs (without extension)
  -d, --dir <DIR>          Output directory to write files to
  -p, --print              Print a table of clades, along with their telomeric sequences
      --log                Output a log file
  -h, --help               Print help
  -V, --version            Print version
```

### Search

`tidk search` will search the genome for an input string. If you know the telomeric repeat of your sequenced organism, this will find it.

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

`tidk plot` will plot a CSV from the output of `tidk search`.

```
SVG plot of TSV generated from search or find.

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

## Cited by:

- Kurbessoian, Tania, et al. "In host evolution of Exophiala dermatitidis in cystic fibrosis lung micro-environment." **BioRxiv** (2022): 2022-09.
- Yin, Denghua, et al. "Gapless genome assembly of East Asian finless porpoise." **Scientific Data** 9.1 (2022): 765.
- Leonard, Guy, et al. "A genome sequence assembly of the phototactic and optogenetic model fungus Blastocladiella emersonii reveals a diversified nucleotide-cyclase repertoire." **Genome Biology and Evolution** 14.12 (2022): evac157.
- Edwards, Richard J., et al. "A phased chromosome-level genome and full mitochondrial sequence for the dikaryotic myrtle rust pathogen, Austropuccinia psidii." **BioRxiv** (2022): 2022-04.