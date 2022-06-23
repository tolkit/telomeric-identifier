use anyhow::Result;
use clap::{crate_version, Arg, Command};
use tidk::{clades::CLADES, explore, finder, min, plot, search, trim, SubCommand};

fn main() -> Result<()> {
    // command line options
    let matches = Command::new("tidk")
        .version(crate_version!())
        .propagate_version(true)
        .arg_required_else_help(true)
        .author("Max Brown <mb39@sanger.ac.uk>")
        .about("A Telomere Identification Toolkit.")
        .subcommand(
            Command::new("find")
                .about("Supply the name of a clade your organsim belongs to, and this submodule will find all telomeric repeat matches for that clade.")
                .arg(
                    Arg::new("fasta")
                        .short('f')
                        .long("fasta")
                        .takes_value(true)
                        .required_unless_present("print")
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::new("window")
                        .short('w')
                        .long("window")
                        .takes_value(true)
                        .required_unless_present("print")
                        .default_value("10000")
                        .help("Window size to calculate telomeric repeat counts in."),
                )
                .arg(
                    Arg::new("clade")
                        .short('c')
                        .long("clade")
                        .takes_value(true)
                        .required_unless_present("print")
                        .possible_values(CLADES)
                        .help("The clade of organism to identify telomeres in."),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .takes_value(true)
                        .required_unless_present("print")
                        .help("Output filename for the CSVs (without extension)."),
                )
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .takes_value(true)
                        .required_unless_present("print")
                        .help("Output directory to write files to."),
                )
                .arg(
                    Arg::new("print")
                        .short('p')
                        .long("print")
                        .help("Print a table of clades, along with their telomeric sequences."),
                )
                .arg(
                    Arg::new("log")
                        .long("log")
                        .help("Output a log file.")
                )
        )
        .subcommand(
            Command::new("explore")
                .about("Use a search of all substrings of length k to query a genome for a telomere sequence.")
                .arg(
                    Arg::new("fasta")
                        .short('f')
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::new("length")
                        .short('l')
                        .long("length")
                        .takes_value(true)
                        .required_unless_present_all(&["minimum", "maximum"])
                        .default_value_if("minimum", None, Some("0"))
                        .help("Length of substring."),
                )
                .arg(
                    Arg::new("minimum")
                        .short('m')
                        .long("minimum")
                        .takes_value(true)
                        .required_unless_present("length")
                        .default_value("5")
                        .help("Minimum length of substring."),
                )
                .arg(
                    Arg::new("maximum")
                        .short('x')
                        .long("maximum")
                        .takes_value(true)
                        .required_unless_present("length")
                        .default_value("12")
                        .help("Maximum length of substring."),
                )
                .arg(
                    Arg::new("threshold")
                        .short('t')
                        .long("threshold")
                        .takes_value(true)
                        .required(false)
                        .default_value("100")
                        .help("Positions of repeats are only reported if they occur sequentially in a greater number than the threshold."),
                )
                .arg(
                    Arg::new("distance")
                        .long("distance")
                        .takes_value(true)
                        .required(false)
                        .default_value("150000")
                        .help("The distance in base pairs from the beginning or end of a chromosome, to report potential telomeric repeats in."),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .takes_value(true)
                        .required(true)
                        .help("Output filename for the TSVs (without extension)."),
                )
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .takes_value(true)
                        .required(true)
                        .help("Output directory to write files to."),
                )
                .arg(
                    Arg::new("extension")
                        .short('e')
                        .long("extension")
                        .takes_value(true)
                        .required_unless_present("print")
                        .default_value("tsv")
                        .possible_values(&["tsv", "bedgraph"])
                        .help("The extension, defining the output type of the file."),
                )
                .arg(
                    Arg::new("verbose")
                        .short('v')
                        .long("verbose")
                        .help("Print verbose output."),
                )
                .arg(
                    Arg::new("log")
                        .long("log")
                        .help("Output a log file.")
                )
        )
        .subcommand(
            Command::new("search")
                .about("Search the input genome with a specific telomeric repeat search string.")
                .arg(
                    Arg::new("fasta")
                        .short('f')
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::new("string")
                        .short('s')
                        .long("string")
                        .takes_value(true)
                        .required(true)
                        .help("Supply a DNA string to query the genome with."),
                )
                .arg(
                    Arg::new("window")
                        .short('w')
                        .long("window")
                        .takes_value(true)
                        .required(false)
                        .default_value("10000")
                        .help("Window size to calculate telomeric repeat counts in."),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .takes_value(true)
                        .required(true)
                        .help("Output filename for the CSVs (without extension)."),
                )
                .arg(
                    Arg::new("dir")
                        .short('d')
                        .long("dir")
                        .takes_value(true)
                        .required(true)
                        .help("Output directory to write files to."),
                )
                .arg(
                    Arg::new("extension")
                        .short('e')
                        .long("extension")
                        .takes_value(true)
                        .required_unless_present("print")
                        .default_value("csv")
                        .possible_values(&["csv", "bedgraph"])
                        .help("The extension, defining the output type of the file."),
                )
                .arg(
                    Arg::new("log")
                        .long("log")
                        .help("Output a log file.")
                )
        )
        .subcommand(
            Command::new("trim")
                .about("Trim a specific telomeric repeat from the input reads and yield reads oriented at the telomere start.")
                .arg(
                    Arg::new("fasta")
                        .short('f')
                        .long("fasta")
                        .takes_value(true)
                        .required(true)
                        .help("The input fasta file."),
                )
                .arg(
                    Arg::new("string")
                        .short('s')
                        .long("string")
                        .takes_value(true)
                        .required(true)
                        .help("Supply a DNA string to trim the reads with."),
                )
                .arg(
                    Arg::new("min_len")
                        .short('l')
                        .long("min_len")
                        .takes_value(true)
                        .required(false)
                        .default_value("1000")
                        .help("Minimum length of trimmed reads."),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .takes_value(true)
                        .required(true)
                        .default_value("tidk-trim")
                        .help("Output filename for the trimmed fasta output."),
                )
                .arg(
                    Arg::new("min_occur")
                        .short('m')
                        .long("min_occur")
                        .takes_value(true)
                        .required(false)
                        .default_value("3")
                        .help("Number of contiguous occurrences of telomeric repeat to start trimming."),
                )
        )
        .subcommand(
            Command::new("plot")
                .about("SVG plot of CSV generated from search or find.")
                // output file name
                .arg(
                    Arg::new("csv")
                        .short('c')
                        .long("csv")
                        .takes_value(true)
                        .required(true)
                        .help("The input CSV file."),
                )
                .arg(
                    Arg::new("height")
                        .long("height")
                        .takes_value(true)
                        .required(false)
                        .default_value("200")
                        .help("The height of subplots (px)."),
                )
                .arg(
                    Arg::new("width")
                        .short('w')
                        .long("width")
                        .takes_value(true)
                        .required(false)
                        .default_value("1000")
                        .help("The width of plot (px)."),
                )
                .arg(
                    Arg::new("output")
                        .short('o')
                        .long("output")
                        .takes_value(true)
                        .required(true)
                        .default_value("tidk-plot")
                        .help("Output filename for the SVG (without extension)."),
                )
        )
        .subcommand(
            Command::new("min")
                .about("Emit the canonical lexicographically minimal DNA string.")
                // output file name
                .arg(
                    Arg::new("DNA string")
                        // .required_unless("file")
                        .multiple_values(true)
                        .help("Input DNA string. Multiple inputs allowed."),
                )
                .arg(
                    Arg::new("fasta")
                        .short('x')
                        .long("fasta")
                        .help("STDIN is in fasta format."),
                )
                .arg(
                    Arg::new("file")
                        .short('f')
                        .long("file")
                        .takes_value(true)
                        // .required_unless("DNA string")
                        .help("The input file."),
                )
        )
        .get_matches();

    // feed command line options to each main function
    match matches.subcommand() {
        Some(("find", matches)) => {
            finder::finder(matches, SubCommand::Find)?;
        }
        Some(("explore", matches)) => {
            explore::explore(matches, SubCommand::Explore)?;
        }
        Some(("search", matches)) => {
            search::search(matches, SubCommand::Search);
        }
        Some(("trim", matches)) => {
            trim::trim(matches);
        }
        Some(("plot", matches)) => {
            plot::plot(matches).unwrap();
        }
        Some(("min", matches)) => {
            min::min_dna_string(matches)?;
        }
        _ => {
            unreachable!()
        }
    }
    
    Ok(())
}
