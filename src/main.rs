// A telomere identification toolkit
// Max Brown 2021, Wellcome Sanger Institute
// mb39@sanger.ac.uk

use clap::{App, Arg};
use std::process;
use tidk::clades::clades::CLADES;
use tidk::explore::explore;
use tidk::finder::finder;
use tidk::min::min;
use tidk::plot::plot;
use tidk::search::search;
use tidk::trim::trim;

fn main() {
    // command line options
    let matches = App::new("TIDK")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("A Telomere Identification Toolkit.")
        .subcommand(
            clap::SubCommand::with_name("find")
                .about("Supply the name of a clade your organsim belongs to, and this submodule will find all telomeric repeat matches for that clade.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required_unless("print")
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::with_name("window")
                        .short("w")
                        .long("window")
                        .takes_value(true)
                        .required_unless("print")
                        .default_value("10000")
                        .help("Window size to calculate telomeric repeat counts in."),
                )
                .arg(
                    Arg::with_name("clade")
                        .short("c")
                        .long("clade")
                        .takes_value(true)
                        .required_unless("print")
                        .possible_values(CLADES)
                        .help("The clade of organism to identify telomeres in."),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .required_unless("print")
                        .default_value("tidk-find")
                        .help("Output filename for the CSVs (without extension)."),
                )
                .arg(
                    Arg::with_name("print")
                        .short("p")
                        .long("print")
                        .help("Print a table of clades, along with their telomeric sequences."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("explore")
                .about("Use a search of all substrings of length k to query a genome for a telomere sequence.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::with_name("length")
                        .short("l")
                        .long("length")
                        .takes_value(true)
                        .required_unless_all(&["minimum", "maximum"])
                        .default_value_if("minimum", None, "0")
                        .help("Length of substring."),
                )
                .arg(
                    Arg::with_name("minimum")
                        .short("m")
                        .long("minimum")
                        .takes_value(true)
                        .required_unless("length")
                        .default_value("5")
                        .help("Minimum length of substring."),
                )
                .arg(
                    Arg::with_name("maximum")
                        .short("x")
                        .long("maximum")
                        .takes_value(true)
                        .required_unless("length")
                        .default_value("12")
                        .help("Maximum length of substring."),
                )
                .arg(
                    Arg::with_name("threshold")
                        .short("t")
                        .long("threshold")
                        .takes_value(true)
                        .required(false)
                        .default_value("100")
                        .help("Positions of repeats are only reported if they occur sequentially in a greater number than the threshold."),
                )
                .arg(
                    Arg::with_name("distance")
                        .short("d")
                        .long("distance")
                        .takes_value(true)
                        .required(false)
                        .default_value("150000")
                        .help("The distance in base pairs from the beginning or end of a chromosome, to report potential telomeric repeats in."),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .required(false)
                        .default_value("tidk-explore")
                        .help("Output filename for the CSVs (without extension)."),
                )
                .arg(
                    Arg::with_name("extension")
                        .short("e")
                        .long("extension")
                        .takes_value(true)
                        .required_unless("print")
                        .default_value("csv")
                        .possible_values(&["csv", "bedgraph"])
                        .help("The extension, defining the output type of the file."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("search")
                .about("Search the input genome with a specific telomeric repeat search string.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::with_name("string")
                        .short("s")
                        .long("string")
                        .takes_value(true)
                        .required(true)
                        .help("Supply a DNA string to query the genome with."),
                )
                .arg(
                    Arg::with_name("window")
                        .short("w")
                        .long("window")
                        .takes_value(true)
                        .required(false)
                        .default_value("10000")
                        .help("Window size to calculate telomeric repeat counts in."),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .required(false)
                        .default_value("tidk-search")
                        .help("Output filename for the CSVs (without extension)."),
                )
                .arg(
                    Arg::with_name("extension")
                        .short("e")
                        .long("extension")
                        .takes_value(true)
                        .required_unless("print")
                        .default_value("csv")
                        .possible_values(&["csv", "bedgraph"])
                        .help("The extension, defining the output type of the file."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("trim")
                .about("Trim a specific telomeric repeat from the input reads and yield reads oriented at the telomere start.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::with_name("string")
                        .short("s")
                        .long("string")
                        .takes_value(true)
                        .required(true)
                        .help("Supply a DNA string to trim the reads with."),
                )
                .arg(
                    Arg::with_name("min_len")
                        .short("l")
                        .long("min_len")
                        .takes_value(true)
                        .required(false)
                        .default_value("1000")
                        .help("Minimum length of trimmed reads."),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .required(false)
                        .default_value("tidk-trim")
                        .help("Output filename for the trimmed fasta output."),
                )
                .arg(
                    Arg::with_name("min_occur")
                        .short("m")
                        .long("min_occur")
                        .takes_value(true)
                        .required(false)
                        .default_value("3")
                        .help("Number of contiguous occurrences of telomeric repeat to start trimming."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("plot")
                .about("SVG plot of CSV generated from search or find.")
                // output file name
                .arg(
                    Arg::with_name("csv")
                        .short("c")
                        .long("csv")
                        .takes_value(true)
                        .required(true)
                        .help("The input CSV file."),
                )
                .arg(
                    Arg::with_name("height")
                        .short("h")
                        .long("height")
                        .takes_value(true)
                        .required(false)
                        .default_value("200")
                        .help("The height of subplots (px)."),
                )
                .arg(
                    Arg::with_name("width")
                        .short("w")
                        .long("width")
                        .takes_value(true)
                        .required(false)
                        .default_value("1000")
                        .help("The width of plot (px)."),
                )
                .arg(
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .required(false)
                        .default_value("tidk-plot")
                        .help("Output filename for the SVG (without extension)."),
                )
        )
        .subcommand(
            clap::SubCommand::with_name("min")
                .about("Emit the canonical lexicographically minimal DNA string.")
                // output file name
                .arg(
                    Arg::with_name("DNA string")
                        // .required_unless("file")
                        .multiple(true)
                        .help("Input DNA string. Multiple inputs allowed."),
                )
                .arg(
                    Arg::with_name("fasta")
                        .short("x")
                        .long("fasta")
                        .help("STDIN is in fasta format."),
                )
                .arg(
                    Arg::with_name("file")
                        .short("f")
                        .long("file")
                        .takes_value(true)
                        // .required_unless("DNA string")
                        .help("The input file."),
                )
        )
        .get_matches();

    // feed command line options to each main function
    let subcommand = matches.subcommand();
    match subcommand.0 {
        "find" => {
            let matches = subcommand.1.unwrap();
            finder::finder(matches);
        }
        "explore" => {
            let matches = subcommand.1.unwrap();
            explore::explore(matches);
        }
        "search" => {
            let matches = subcommand.1.unwrap();
            search::search(matches);
        }
        "trim" => {
            let matches = subcommand.1.unwrap();
            trim::trim(matches);
        }
        "plot" => {
            let matches = subcommand.1.unwrap();
            plot::plot(matches).unwrap();
        }
        "min" => {
            let matches = subcommand.1.unwrap();
            min::min_dna_string(matches);
        }
        _ => {
            println!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            process::exit(1);
        }
    }
}
