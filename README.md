# A Telomere Identification toolKit (tidk)

tidk is a toolkit to identify and visualise telomeric repeats for the Darwin Tree of Life genomes. tidk works especially well on chromosomal genomes, but can also work on PacBio HiFi reads well (see <a href="https://github.com/tolkit/a-telomeric-repeat-database">the telomeric repeat database</a> for many examples). There are a few modules in the tool, which may be useful to anyone investigating telomeric repeat sequences in a genome.

1. `explore` - tries to find the telomeric repeat unit in the genome.
2. `find` and `search` are essentially the same. They identify a repeat sequence in windows across the genome. `find` uses an in-built table of telomeric repeats, in `search` you supply your own.
3. `plot` does what is says on the tin, and plots the csv output of `find` or `search` as an SVG.
4. `min` returns the lexicographically minimal string of an input DNA string(s) or file of DNA sequences.
5. `trim` trims fasta sequences using a supplied base repeat string.

## Install

As with other Rust projects, you have to complile yourself. <a href="https://www.rust-lang.org/tools/install">Download rust</a>, clone this repo, and then run:

`cargo build --release`

Compiling takes anywhere from 1-6 minutes from fresh (tested on the farm). The executable will be at the location `./target/release/tidk`.

## Usage

### Overall

```
TIDK 0.1.5
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
    min        Emit the canonical lexicographically minimal DNA string.
    plot       SVG plot of CSV generated from search or find.
    search     Search the input genome with a specific telomeric repeat search string.
    trim       Trim a specific telomeric repeat from the input reads and yield reads oriented at the telomere start.
```

### Explore 

`tidk explore` will identify all sequences of length k, which repeat at least twice throughout a genome. Repeats of high number toward the beginning or end of sequences are likely candidates for telomeric repeats. The reported repeats are the lexicographically most minimal of all possible string rotations of the telomeric repeat, in both forward and reverse complement forms.

It outputs either a csv or bedgraph of potential telomeric repeats and their locations, in addition to a text file of the potential telomeric repeat sequences, and how often they are found in the genome.

For example:
`tidk explore --fasta fastas/iyBomHort1_1.20210303.curated_primary.fa --minimum 5 --maximum 12 -o test_dist -t 500` searches the genome for repeats from length 5 to length 12 sequentially (definite potential to be made concurrent) on the freshly minted <a href="https://www.ebi.ac.uk/ena/browser/view/PRJEB43539"><i>Bombus hortorum</i> genome</a>.

```
tidk-explore
Use a search of all substrings of length k to query a genome for a telomere sequence.

USAGE:
    tidk explore [OPTIONS] --extension <extension> --fasta <fasta> --length <length> --maximum <maximum> --minimum <minimum>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -d, --distance <distance>      The distance in base pairs from the beginning or end of a chromosome, to report
                                   potential telomeric repeats in. [default: 150000]
    -e, --extension <extension>    The extension, defining the output type of the file. [default: csv]  [possible
                                   values: csv, bedgraph]
    -f, --fasta <fasta>            The input fasta file.
    -l, --length <length>          Length of substring.
    -x, --maximum <maximum>        Maximum length of substring. [default: 12]
    -m, --minimum <minimum>        Minimum length of substring. [default: 5]
    -o, --output <output>          Output filename for the CSVs (without extension). [default: tidk-explore]
    -t, --threshold <threshold>    Positions of repeats are only reported if they occur sequentially in a greater number
                                   than the threshold. [default: 100]
```

### Find

`tidk find` will take an input clade, and match the known telomeric repeat for that clade (or repeats plural) and search the genome. Uses the <a href="http://telomerase.asu.edu/sequences_telomere.html">telomeric repeat database</a>. As more telomeric repeats are found and added, the dictionary of sequences used will increase (perhaps there is a more elegant way to parse the command line input?).

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
    tidk search [OPTIONS] --extension <extension> --fasta <fasta> --string <string>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -e, --extension <extension>    The extension, defining the output type of the file. [default: csv]  [possible
                                   values: csv, bedgraph]
    -f, --fasta <fasta>            The input fasta file.
    -o, --output <output>          Output filename for the CSVs (without extension). [default: tidk-search]
    -s, --string <string>          Supply a DNA string to query the genome with.
    -w, --window <window>          Window size to calculate telomeric repeat counts in. [default: 10000]
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
    -c, --csv <csv>          The input CSV file.
    -h, --height <height>    The height of subplots (px). [default: 200]
    -o, --output <output>    Output filename for the SVG (without extension). [default: tidk-plot]
    -w, --width <width>      The width of plot (px). [default: 1000]
```

As an example on the ol' Square Spot Rustic <i>Xestia xanthographa</i>:

```bash
tidk find -f fastas/ilXesXant1_1.20201023.curated_primary.fa -c lepidoptera -o Xes

tidk plot -c finder/Xes_telomeric_repeat_windows.csv -o ilXes -h 120 -w 800
```
<img src="./ilXes.svg">

### Min

`tidk min` returns the lexicographically minimal DNA string given an input. Useful for universally comparing repeats using this tool. Surprisingly, it works on even chromosomal length DNA strings quickly, not that you would want to do that...

Examples:

- `tidk min AATGCG` or process multiple `tidk min AATGCG AAGGTTC GGTTAAT`
- `tidk min -f input.fa` proccesses input fasta, and outputs fasta. Otherwise `tidk min -f input.txt` reads lines and outputs.
- `echo "AATTGC" | tidk min`, pipes work.
- `cat input.fasta | tidk min -x` outputs a fasta. Otherwise `cat input.txt | tidk min` reads lines and outputs. 

Note `cat` here is redundant and creates extra work, but it just shows the piping in action.

```
tidk-min 
Emit the canonical lexicographically minimal DNA string.

USAGE:
    tidk min [FLAGS] [OPTIONS] [DNA string]...

FLAGS:
    -x, --fasta      STDIN is in fasta format.
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --file <file>    The input file.

ARGS:
    <DNA string>...    Input DNA string. Multiple inputs allowed.
```

### Trim

`tidk trim` - a rust port of https://github.com/pgonzale60/telomeric-trim. Thanks Pablo!

```
tidk-trim 
Trim a specific telomeric repeat from the input reads and yield reads oriented at the telomere start.

USAGE:
    tidk trim [OPTIONS] --fasta <fasta> --string <string>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -f, --fasta <fasta>            The input fasta file.
    -l, --min_len <min_len>        Minimum length of trimmed reads. [default: 1000]
    -m, --min_occur <min_occur>    Number of contiguous occurrences of telomeric repeat to start trimming. [default: 3]
    -o, --output <output>          Output filename for the trimmed fasta output. [default: tidk-trim]
    -s, --string <string>          Supply a DNA string to trim the reads with.
```