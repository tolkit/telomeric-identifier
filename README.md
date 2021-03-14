# A Telomere Identification toolKit (tidk)

Simple and fast de-novo identification of telomeric repeats (`explore`), finding frequency of a telomeric sequence by clade in windows (`find`), or a user defined sequence in windows (`search`).

## Install

As with other Rust projects, you have to complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, and then run:

`cargo build --release`

Compiling takes ~2 minutes. The executable will be in `./target/release/tidk`.

## Usage

### Overall

```
TIDK 0.1.3
Max Brown <mb39@sanger.ac.uk>
A Telomere Identification Toolkit.

USAGE:
    tidk [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    explore    Use a search of all substrings of length k to query a genome for a telomere sequence.
    find       Supply the name of a clade your organsim belongs to, and this submodule will find all telomeric
               repeat matches for that clade.
    help       Prints this message or the help of the given subcommand(s)
    plot       SVG plot of CSV generated from search or find.
    search     Search the input genome with a specific telomeric repeat search string.
```

### Explore 

`tidk explore` will identify all sequences of length k, which repeat at least twice throughout a genome. Repeats of high number toward the beginning or end of sequences are likely candidates for telomeric repeats. The longest run of G's seems to be at the end of the canonical telomeric repeat unit.

For example:
`tidk explore --fasta fastas/iyBomHort1_1.20210303.curated_primary.fa --minimum 5 --maximum 12 -o test_dist -t 500` searches the genome for repeats from length 5 to length 12 sequentially (definite potential to be made concurrent) on the freshly minted <a href="https://www.ebi.ac.uk/ena/browser/view/PRJEB43539"><i>Bombus hortorum</i> genome</a>.

```
tidk-explore
Use a search of all substrings of length k to query a genome for a telomere sequence.

USAGE:
    tidk explore [OPTIONS] --fasta <fasta> --length <length> --maximum <maximum> --minimum <minimum>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -d, --distance <distance>      The distance in base pairs from the beginning or end of a chromosome, to report
                                   potential telomeric repeats in. [default: 150000]
    -f, --fasta <fasta>            The input fasta file.
    -l, --length <length>          Length of substring.
    -x, --maximum <maximum>        Maximum length of substring. [default: 12]
    -m, --minimum <minimum>        Minimum length of substring. [default: 5]
    -o, --output <output>          Output filename for the CSVs (without extension). [default: tidk-explore]
    -t, --threshold <threshold>    Positions of repeats are only reported if they occur sequentially in a greater number
                                   than the threshold. [default: 100]
```

### Find

`tidk find` will take an input clade, and match the known telomeric repeat for that clade (or repeats plural) and search the genome. Uses the <a href="http://telomerase.asu.edu/sequences_telomere.html">telomeric repeat database</a>. As more telomeric repeats are found and added, perhaps there is a more elegant way to parse the command line input.

```
tidk-find
Supply the name of a clade your organsim belongs to, and this submodule will find all telomeric repeat matches for that
clade.

USAGE:
    tidk find [FLAGS] --clade <clade> --fasta <fasta> --output <output> --window <window>

FLAGS:
    -h, --help       Prints help information
    -p, --print      Print a table of clades, along with their telomeric sequences.
    -V, --version    Prints version information

OPTIONS:
    -c, --clade <clade>      The clade of organism to identify telomeres in. [possible values: vertebrate, ascidian,
                             echinodermata, mollusca, coleoptera, hymenoptera, lepidoptera, megaloptera, trichoptera,
                             neuroptera, blattodea, orthoptera, nematoda, amoeba, plants, ciliates]
    -f, --fasta <fasta>      The input fasta file.
    -o, --output <output>    Output filename for the CSVs (without extension). [default: tidk-find]
    -w, --window <window>    Window size to calculate telomeric repeat counts in. [default: 10000]
```

### Search

`tidk search` will search the genome for an input string. If you know the telomeric repeat of your sequenced organism, this will hopefully find it.

```
tidk-search
Search the input genome with a specific telomeric repeat search string.

USAGE:
    tidk search [OPTIONS] --fasta <fasta> --string <string>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>      The input fasta file.
    -o, --output <output>    Output filename for the CSVs (without extension). [default: tidk-search]
    -s, --string <string>    Supply a DNA string to query the genome with.
    -w, --window <window>    Window size to calculate telomeric repeat counts in. [default: 10000]
```

### Plot

`tidk plot` will plot a CSV from the output of `tidk search`. Working on plotting for `tidk find` (i.e. extending to multiple telomeric repeat sequences in same CSV).

```
tidk-plot
SVG plot of CSV generated from search or find.

USAGE:
    tidk plot [OPTIONS] --csv <csv>

FLAGS:
        --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --csv <csv>                                The input CSV file.
    -h, --height <height>                          The height of subplots (px). [default: 200]
    -l, --length_chromosome <length_chromosome>
            Chromosomes shorter than this length will be excluded from the plot. Useful for unplaced scaffold exclusion.
            [default: 1000000]
    -o, --output <output>                          Output filename for the SVG (without extension). [default: tidk-plot]
    -w, --width <width>                            The width of plot (px). [default: 1000]
```

As an example on the ol' Square Spot Rustic <i>Xestia xanthographa</i>:

```bash
tidk find -f fastas/ilXesXant1_1.20201023.curated_primary.fa -c lepidopterac -o Xes

tidk plot -c finder/Xes_telomeric_repeat_windows.csv -o ilXes -h 120 -w 800
```
<img src="./ilXes.svg">

## TODO's

- Implement IUPAC in `tidk search`?
- Better summary of `tidk explore`?