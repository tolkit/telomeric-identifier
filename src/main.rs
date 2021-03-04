// how to build a telomere identification toolkit?

use clap::{App, Arg};
use std::process;
use tidk::clades::clades::CLADES;
use tidk::explore::explore;
use tidk::finder::finder;
use tidk::plot::plot;
use tidk::search::search;

fn main() {
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
                // maybe we can have default as 6-8?
                .arg(
                    Arg::with_name("length")
                        .short("l")
                        .long("length")
                        .takes_value(true)
                        .required(true)
                        .help("Length of substring."),
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
                    Arg::with_name("output")
                        .short("o")
                        .long("output")
                        .takes_value(true)
                        .required(false)
                        .default_value("tidk-explore")
                        .help("Output filename for the CSVs (without extension)."),
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
                        .required_unless("print")
                        .default_value("10000")
                        .help("Window size to calculate telomeric repeat counts in."),
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
                    Arg::with_name("length_chromosome")
                        .short("l")
                        .long("length_chromosome")
                        .takes_value(true)
                        .required(false)
                        .default_value("1000000") // 1 million bases
                        .help("Chromosomes shorter than this length will be excluded from the plot. Useful for unplaced scaffold exclusion."),
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
        "plot" => {
            let matches = subcommand.1.unwrap();
            // handle this potential error better please
            plot::plot(matches).unwrap();
        }
        _ => {
            println!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            process::exit(1);
        }
    }
}
