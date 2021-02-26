# Identify telomeric repeats

## Aims

Simple and fast de-novo identification of telomeric repeats (explore), finding frequency of a telomeric sequence by clade in windows (find), or a user defined sequence in windows (search; not yet implemented).

## Install

As with other Rust projects, you have to complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, and then run:

`cargo build --release`

Compiling takes ~2 minutes.

## Usage

```
TIDK 0.1.2
Max Brown <mb39@sanger.ac.uk>
A Telomere Identification Toolkit.

USAGE:
    tidk [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    explore    Use a search of all substrings of length k to query a genome for a telomere sequence.
    find       Supply the name of a clade your organsim belongs to.
    help       Prints this message or the help of the given subcommand(s)
    search     Search the input genome with a specific telomeric repeat search string.
```

`telomeric-identifier explore` will identify all sequences of length k, which repeat at least twice throughout a genome. Repeats of high number toward the beginning or end of sequences are likely candidates for telomeric repeats. Prototype code done.

```
tidk-explore
Use a search of all substrings of length k to query a genome for a telomere sequence.

USAGE:
    tidk explore [OPTIONS] --fasta <fasta> --length <length>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>            The input fasta file.
    -l, --length <length>          Length of substring.
    -o, --output <output>          Output filename for the CSVs (without extension). [default: tidk-explore]
    -t, --threshold <threshold>    Positions of repeats are only reported if they occur sequentially in a greater number than the threshold. [default: 100]
```

`telomeric-identifier find` will take an input clade, and match the known telomeric repeat for that clade (or repeats plural) and search the genome. Not yet done. Uses the telomeric repeat database -> http://telomerase.asu.edu/sequences_telomere.html.

```
tidk-find
Supply the name of a clade your organsim belongs to.

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

`telomeric-identifier search` will search the genome for an input string. If you know the telomeric repeat of your sequenced organism, this will hopefully find it. Protoype code done but needs much tidying.

```
tidk-search
Search the input genome with a specific telomeric repeat search string.

USAGE:
    tidk search --fasta <fasta> --string <string>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>      The input fasta file.
    -s, --string <string>    Supply a DNA string to query the genome with.
```