use anyhow::Result;
use clap::{arg, builder::ArgPredicate, crate_version, value_parser, Arg, Command};
use std::path::PathBuf;
use tidk::{build, clades::get_clades, explore, finder, plot, search, SubCommand};

fn main() -> Result<()> {
    // command line options
    let matches = Command::new("tidk")
        .version(crate_version!())
        .propagate_version(true)
        .arg_required_else_help(true)
        .author("Max Brown <max.carter-brown@aru.ac.uk>")
        .about("A Telomere Identification Toolkit.")
        .subcommand(
            Command::new("build")
                .about("Build the reference database of telomeric repeat sequences. This is required for the 'find' subcommand.")
        )
        .subcommand(
            Command::new("find")
                .about("Supply the name of a clade your organsim belongs to, and this submodule will find all telomeric repeat matches for that clade.")
                .arg(
                    Arg::new("fasta")
                        .value_name("FASTA")
                        .value_parser(value_parser!(PathBuf))
                        .help("The input fasta file")
                        .required_unless_present("print")
                )
                .arg(
                    // no longer required.
                    arg!(-w --window [WINDOW] "Window size to calculate telomeric repeat counts in")
                        .value_parser(value_parser!(usize))
                        .default_value("10000")
                )
                .arg(
                    arg!(-c --clade <CLADE> "The clade of organism to identify telomeres in")
                        .required_unless_present("print")
                        .value_parser(get_clades().unwrap_or_else(|_| -> Vec<_> {
                            eprintln!("Warning! No clades found in the database. Run 'tidk build' to fetch the latest data.\n");
                            vec![]
                        }))
                )
                .arg(
                    arg!(-o --output <OUTPUT> "Output filename for the TSVs (without extension)")
                        .value_parser(value_parser!(PathBuf))
                        .required_unless_present("print")
                )
                .arg(
                    arg!(-d --dir <DIR> "Output directory to write files to")
                        .required_unless_present("print")
                        .value_parser(value_parser!(PathBuf))
                )
                .arg(
                    arg!(-p --print "Print a table of clades, along with their telomeric sequences")
                        .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    arg!(--log "Output a log file")
                        .action(clap::ArgAction::SetTrue)
                )
        )
        .subcommand(
            Command::new("explore")
                .about("Use a range of kmer sizes to find potential telomeric repeats.\nOne of either length, or minimum and maximum must be specified.")
                .arg(
                    Arg::new("fasta")
                        .value_name("FASTA")
                        .value_parser(value_parser!(PathBuf))
                        .required(true)
                        .help("The input fasta file")
                )
                .arg(
                    arg!(-l --length [LENGTH] "Length of substring")
                        .required_unless_present_all(["minimum", "maximum"])
                        .default_value_if("minimum", ArgPredicate::IsPresent, Some("0"))
                        .value_parser(value_parser!(usize))
                )
                .arg(
                    arg!(-m --minimum [MINIMUM] "Minimum length of substring")
                        .required_unless_present("length")
                        .value_parser(value_parser!(usize))
                        .default_value("5")
                )
                .arg(
                    arg!(-x --maximum [MAXIMUM] "Maximum length of substring")
                        .required_unless_present("length")
                        .value_parser(value_parser!(usize))
                        .default_value("12")
                )
                .arg(
                    arg!(-t --threshold [THRESHOLD] "Positions of repeats are only reported if they occur sequentially in a greater number than the threshold")
                        .value_parser(value_parser!(i32))
                        .default_value("100")
                )
                .arg(
                    arg!(--distance [DISTANCE] "The distance from the end of the chromosome as a proportion of chromosome length. Must range from 0-0.5.")
                        .value_parser(value_parser!(f64))
                        .default_value("0.01")
                )
                .arg(
                    arg!(-v --verbose "Print verbose output.")
                        .action(clap::ArgAction::SetTrue)
                )
                .arg(
                    arg!(--log "Output a log file.")
                        .action(clap::ArgAction::SetTrue)
                )
        )
        .subcommand(
            Command::new("search")
                .about("Search the input genome with a specific telomeric repeat search string.")
                .arg(
                    Arg::new("fasta")
                        .value_name("FASTA")
                        .value_parser(value_parser!(PathBuf))
                        .required(true)
                        .help("The input fasta file")
                )
                .arg(
                    arg!(-s --string <STRING> "The DNA string to query the genome with")
                        .required(true)
                )
                .arg(
                    arg!(-w --window [WINDOW] "Window size to calculate telomeric repeat counts in")
                        .value_parser(value_parser!(usize))
                        .default_value("10000")
                )
                .arg(
                    arg!(-o --output <OUTPUT> "Output filename for the TSVs (without extension)")
                        .required(true)
                )
                .arg(
                    arg!(-d --dir <DIR> "Output directory to write files to")
                        .required(true)
                        .value_parser(value_parser!(PathBuf))
                )
                .arg(
                    arg!(-e --extension [EXTENSION] "The extension, defining the output type of the file")
                        .default_value("tsv")
                        .value_parser(["tsv", "bedgraph"])
                )
                .arg(
                    arg!(--log "Output a log file")
                        .action(clap::ArgAction::SetTrue)
                )
        )
        .subcommand(
            Command::new("plot")
                .about("SVG plot of TSV generated from tidk search.")
                // output file name
                .arg(
                    arg!(-t --tsv <TSV> "The input TSV file")
                        .value_parser(value_parser!(PathBuf))
                        .required(true)
                )
                .arg(
                    arg!(--height [HEIGHT] "The height of subplots (px).")
                        .value_parser(value_parser!(i32))
                        .default_value("200")
                )
                .arg(
                    arg!(-w --width [WIDTH] "The width of plot (px)")
                        .value_parser(value_parser!(i32))
                        .default_value("1000")
                )
                .arg(
                    arg!(-o --output [OUTPUT] "Output filename for the SVG (without extension)")
                        .value_parser(value_parser!(PathBuf))
                        .default_value("tidk-plot")
                )
                .arg(
                    arg!(--fontsize [FONT_SIZE] "The font size of the axis labels in the plot")
                        .value_parser(value_parser!(i32))
                        .default_value("12")
                )
                .arg(
                    arg!(--strokewidth [STROKE_WIDTH] "The stroke width of the line graph in the plot")
                        .value_parser(value_parser!(i32))
                        .default_value("2")
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
            search::search(matches, SubCommand::Search)?;
        }
        Some(("plot", matches)) => {
            plot::plot(matches)?;
        }
        Some(("build", _)) => {
            build::fetch_and_save_data()?;
        }
        _ => {
            unreachable!()
        }
    }

    Ok(())
}
