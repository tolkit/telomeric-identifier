# Identify telomeric repeats de-novo, cladewise or using input string

## Aims

Motivated by this interesting bioinfomatic problem at the moment. In active development, so not a lot to see yet. It aims however, to be simple and fast. 

## To-do's

- What output format is most useful? CSV?
- Window approach for the find/search subcommands?

## Usage

There are three submodules. It takes an input fasta assembly, along with one other argument (so far) depending on the subcommand.

```
TIDK 0.1.0
Max Brown <mb39@sanger.ac.uk>
A Telomere Identification Toolkit.

USAGE:
    telomeric-identifier [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

SUBCOMMANDS:
    explore    Use a search of all substrings of length k to query a genome for a telomere sequence.
    find       Supply a DNA string or flag the clade your organsim belongs to.
    help       Prints this message or the help of the given subcommand(s)
    search     Search the input genome with a specific telomeric repeat search string.
```

`telomeric-identifier explore` will identify all sequences of length k, which repeat at least twice throughout a genome. Repeats of high number toward the beginning or end of sequences are likely candidates for telomeric repeats. Prototype code done.

```
telomeric-identifier-explore
Use a search of all substrings of length k to query a genome for a telomere sequence.

USAGE:
    telomeric-identifier explore --fasta <fasta> --length <length>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>      The input fasta file.
    -l, --length <length>    Length of substring.
```

`telomeric-identifier search` will search the genome for an input string. If you know the telomeric repeat of your sequenced organism, this will hopefully find it. Protoype code done but needs tidying.

```
telomeric-identifier-search
Search the input genome with a specific telomeric repeat search string.

USAGE:
    telomeric-identifier search --fasta <fasta> --string <string>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>      The input fasta file.
    -s, --string <string>    Supply a DNA string to query the genome with.
```

`telomeric-identifier find` will take an input clade, and match the known telomeric repeat for that clade (or repeats plural) and search the genome. Not yet done. Uses the telomeric repeat database -> http://telomerase.asu.edu/sequences_telomere.html.

```
telomeric-identifier-find
Supply a DNA string or flag the clade your organsim belongs to.

USAGE:
    telomeric-identifier find --clade <clade> --fasta <fasta>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --clade <clade>    The clade of organism to identify telomeres in. [values: vertebrate, ascidian, echinodermata,
                           mollusca, coleoptera, hymenoptera, lepidoptera, megaloptera, trichoptera, neuroptera,
                           blattodea, orthoptera, nematoda, amoeba, plants, ciliates]
    -f, --fasta <fasta>    The input fasta file.
```

