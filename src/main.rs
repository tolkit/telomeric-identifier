// how to build a telomere identification toolkit?

use clap::{App, Arg};
use std::process;
use telomeric_identifier::explore::explore;
use telomeric_identifier::finder::finder;
use telomeric_identifier::search::search;

fn main() {
    let matches = App::new("TIDK")
        .version(clap::crate_version!())
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("A Telomere Identification Toolkit.")
        .subcommand(
            clap::SubCommand::with_name("find")
                .about("Supply the name of a clade your organsim belongs to.")
                .arg(
                    Arg::with_name("fasta")
                        .short("f")
                        .long("fasta")
                        .takes_value(true)
                        .required_unless("print")
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::with_name("clade")
                        .short("c")
                        .long("clade")
                        .takes_value(true)
                        .required_unless("print")
                        .possible_values(&["vertebrate", "ascidian", "echinodermata", "mollusca", "coleoptera", "hymenoptera", "lepidoptera", "megaloptera", "trichoptera", "neuroptera", "blattodea", "orthoptera", "nematoda", "amoeba", "plants", "ciliates"])
                        .help("The clade of organism to identify telomeres in."),
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
        )
        .get_matches();

    // parse command line options
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
        _ => {
            println!("Subcommand invalid, run with '--help' for subcommand options. Exiting.");
            process::exit(1);
        }
    }
}
