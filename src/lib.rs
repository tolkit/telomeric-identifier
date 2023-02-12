use anyhow::Result;
use chrono::Local;
use clap::crate_version;
use std::{io::Write, path::PathBuf};

/// A module where the clades are defined, and their
/// respective telomeric repeats are enumerated.
pub mod clades;
/// The entry point for the `tidk explore` subcommand.
pub mod explore;
/// The entry point for the `tidk find` subcommand.
pub mod finder;
/// The entry point for the `tidk min` subcommand.
pub mod min;
/// Functions to plot output from `tidk search` and
/// `tidk find`.
pub mod plot;
/// The entry point for the `tidk search` subcommand.
pub mod search;
/// The entry point for `tidk trim` subcommand.
pub mod trim;
/// Module for utilities.
pub mod utils;

/// Three possible subcommands.
pub enum SubCommand {
    Find,
    Explore,
    Search,
}

/// A date format.
const DATE_FORMAT_STR: &str = "%Y-%m-%d: %H:%M:%S";

// this is not the optimal way to do this... but oh well.
// add optional log file directory
impl SubCommand {
    /// Make a log dependent on the subcommand that was run.
    pub fn log(&self, matches: &clap::ArgMatches) -> Result<()> {
        // only if log CLI arg is present
        if matches.get_flag("log") {
            match self {
                SubCommand::Find => {
                    let output = matches
                        .get_one::<PathBuf>("output")
                        .expect("errored by clap");
                    let outdir = matches.get_one::<PathBuf>("dir").expect("errored by clap");
                    let input_fasta = matches
                        .get_one::<PathBuf>("fasta")
                        .expect("errored by clap");
                    let clade = matches.get_one::<String>("clade").expect("errored by clap");
                    let clade_info = clades::return_telomere_sequence(&clade);
                    let window_size = *matches.get_one::<usize>("window").expect("errored by clap");

                    let file_name = format!(
                        "{}/{}{}",
                        outdir.display(),
                        output.display(),
                        "_telomeric_repeat_windows.csv"
                    );

                    let log_string = format!(
                        r#"tidk version: {}
Log information for output file: {}
Date: {}
`tidk find` was run with the following parameters:
    Input fasta: {}
    Window size: {}
    Clade chosen: {}
    Telomeric repeats queried: {}"#,
                        crate_version!(),
                        format!("{}", file_name),
                        Local::now().format(DATE_FORMAT_STR),
                        input_fasta.display(),
                        window_size,
                        clade,
                        clade_info.seq.0.join(", ")
                    );

                    // create file
                    let log_file_name =
                        format!("{}/{}{}", outdir.display(), output.display(), ".log");
                    let log_file = std::fs::File::create(&log_file_name)?;
                    let mut log_file = std::io::LineWriter::new(log_file);

                    writeln!(log_file, "{}", log_string)?;

                    Ok(eprintln!("[+]\tLog file written to: {}", log_file_name))
                }
                SubCommand::Explore => {
                    let outdir = matches.get_one::<PathBuf>("dir").expect("errored by clap");
                    let input_fasta = matches
                        .get_one::<PathBuf>("fasta")
                        .expect("errored by clap");
                    let length = matches.get_one::<usize>("length");

                    let minimum = matches.get_one::<usize>("minimum");
                    let maximum = matches.get_one::<usize>("maximum");

                    let threshold = matches.get_one::<i32>("threshold");
                    let output = matches
                        .get_one::<PathBuf>("output")
                        .expect("errored by clap");
                    let extension = matches
                        .get_one::<String>("extension")
                        .expect("errored by clap");

                    let dist_from_chromosome_end = matches.get_one::<usize>("distance");

                    // create file
                    let explore_file_name = format!(
                        "{}/{}{}{}",
                        outdir.display(),
                        output.display(),
                        "_telomeric_locations.",
                        extension
                    );

                    let putative_telomeric_file =
                        format!("./explore/{}{}", output.display(), ".txt");

                    let log_string = format!(
                        r#"tidk version: {}
Log information for output files: {}, {}
Date: {}
`tidk explore` was run with the following parameters:
    Input fasta: {}
    Explored telomeric repeat units of length: {}
    Or from length: {}
    To length: {}
    Threshold: {}
    Searching at {}bp distance from chromosome end"#,
                        crate_version!(),
                        explore_file_name,
                        putative_telomeric_file,
                        Local::now().format(DATE_FORMAT_STR),
                        input_fasta.display(),
                        {
                            if let Some(l) = length {
                                l.to_string()
                            } else {
                                "None".into()
                            }
                        },
                        {
                            if let Some(min) = minimum {
                                min.to_string()
                            } else {
                                "None".into()
                            }
                        },
                        {
                            if let Some(max) = maximum {
                                max.to_string()
                            } else {
                                "None".into()
                            }
                        },
                        // safely unwrap
                        threshold.unwrap(),
                        // safely unwrap
                        dist_from_chromosome_end.unwrap(),
                    );

                    // create file
                    let log_file_name =
                        format!("{}/{}{}", outdir.display(), output.display(), ".log");
                    let log_file = std::fs::File::create(&log_file_name)?;
                    let mut log_file = std::io::LineWriter::new(log_file);

                    writeln!(log_file, "{}", log_string)?;

                    Ok(eprintln!("[+]\tLog file written to: {}", log_file_name))
                }
                SubCommand::Search => {
                    let input_fasta = matches
                        .get_one::<PathBuf>("fasta")
                        .expect("errored by clap");
                    let telomeric_repeat = matches
                        .get_one::<String>("string")
                        .expect("errored by clap");
                    let extension = matches
                        .get_one::<String>("extension")
                        .expect("errored by clap");

                    let window_size = matches.get_one::<usize>("window").expect("errored by clap");
                    let outdir = matches.get_one::<PathBuf>("dir").expect("errored by clap");
                    let output = matches
                        .get_one::<PathBuf>("output")
                        .expect("errored by clap");

                    // create file
                    let file_name = format!(
                        "{}/{}{}{}",
                        outdir.display(),
                        output.display(),
                        "_telomeric_repeat_windows.",
                        extension
                    );

                    let log_string = format!(
                        r#"tidk version: {}
Log information for output file: {}
Date: {}
`tidk search` was run with the following parameters:
    Input fasta: {}
    Telomeric repeat search string: {}
    Window size: {}
                    "#,
                        crate_version!(),
                        file_name,
                        Local::now().format(DATE_FORMAT_STR),
                        input_fasta.display(),
                        telomeric_repeat,
                        window_size
                    );

                    // create file
                    let log_file_name =
                        format!("{}/{}{}", outdir.display(), output.display(), ".log");
                    let log_file = std::fs::File::create(&log_file_name)?;
                    let mut log_file = std::io::LineWriter::new(log_file);

                    writeln!(log_file, "{}", log_string)?;

                    Ok(eprintln!("[+]\tLog file written to: {}", log_file_name))
                }
            }
        } else {
            Ok(())
        }
    }
}
