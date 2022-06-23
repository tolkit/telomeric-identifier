use anyhow::{Context, Result};
use chrono::Local;
use clap::crate_version;
use std::io::Write;

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
const DATE_FORMAT_STR: &'static str = "%Y-%m-%d: %H:%M:%S";

// this is not the optimal way to do this... but oh well.
// add optional log file directory
impl SubCommand {
    /// Make a log dependent on the subcommand that was run.
    pub fn log(&self, matches: &clap::ArgMatches) -> Result<()> {
        // only if log CLI arg is present
        if matches.is_present("log") {
            match self {
                SubCommand::Find => {
                    let output = matches
                        .value_of("output")
                        .context("Could not get the value of `output`")?;

                    let outdir = matches
                        .value_of("dir")
                        .context("Could not get the value of `dir`")?;
                    let input_fasta = matches
                        .value_of("fasta")
                        .context("Could not get the value of `fasta`")?;
                    let clade: String = matches
                        .value_of_t("clade")
                        .context("Could not parse `clade` as a String.")?;
                    let clade_info = clades::return_telomere_sequence(&clade);
                    let window_size: usize = matches
                        .value_of_t("window")
                        .context("Could not parse `window` as usize.")?;

                    let file_name =
                        format!("{}/{}{}", outdir, output, "_telomeric_repeat_windows.csv");

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
                        input_fasta,
                        window_size,
                        clade,
                        clade_info.seq.0.join(", ")
                    );

                    // create file
                    let log_file_name = format!("{}/{}{}", outdir, output, ".log");
                    let log_file = std::fs::File::create(&log_file_name)?;
                    let mut log_file = std::io::LineWriter::new(log_file);

                    writeln!(log_file, "{}", log_string)?;

                    Ok(eprintln!("[+]\tLog file written to: {}", log_file_name))
                }
                SubCommand::Explore => {
                    let outdir = matches
                        .value_of("dir")
                        .context("Could not get the value of `dir`")?;
                    let input_fasta = matches
                        .value_of("fasta")
                        .context("Could not get the value of `fasta`")?;
                    let length: Option<String> = matches.value_of_t("length").ok();

                    let minimum: Option<String> = matches.value_of_t("minimum").ok();
                    let maximum: Option<String> = matches.value_of_t("maximum").ok();

                    let threshold: Option<String> = matches.value_of_t("threshold").ok();
                    let output = matches
                        .value_of("output")
                        .context("Could not get the value of `output`")?;
                    let extension: String = matches
                        .value_of_t("extension")
                        .context("Could not parse `extension` as a String.")?;

                    let dist_from_chromosome_end: Option<String> =
                        matches.value_of_t("distance").ok();

                    // create file
                    let explore_file_name = format!(
                        "{}/{}{}{}",
                        outdir, output, "_telomeric_locations.", extension
                    );

                    let putative_telomeric_file = format!("./explore/{}{}", output, ".txt");

                    let log_string = format!(
                        r#"tidk version: {}
Log information for output files: {}
Date: {}
`tidk explore` was run with the following parameters:
    Input fasta: {}
    Explored telomeric repeat units of length: {}
    Or from length: {}
    To length: {}
    Threshold: {}
    Searching at {}bp distance from chromosome end"#,
                        crate_version!(),
                        format!("{}, {}", explore_file_name, putative_telomeric_file),
                        Local::now().format(DATE_FORMAT_STR),
                        input_fasta,
                        length.unwrap_or("-".into()),
                        minimum.unwrap_or("-".into()),
                        maximum.unwrap_or("-".into()),
                        threshold.unwrap_or("No threshold specified".into()),
                        dist_from_chromosome_end.unwrap_or("No distance specified".into()),
                    );

                    // create file
                    let log_file_name = format!("{}/{}{}", outdir, output, ".log");
                    let log_file = std::fs::File::create(&log_file_name)?;
                    let mut log_file = std::io::LineWriter::new(log_file);

                    writeln!(log_file, "{}", log_string)?;

                    Ok(eprintln!("[+]\tLog file written to: {}", log_file_name))
                }
                SubCommand::Search => {
                    let input_fasta = matches
                        .value_of("fasta")
                        .context("Could not get the value of `fasta`")?;

                    let telomeric_repeat: String = matches
                        .value_of_t("string")
                        .context("Could not parse `string` as a String.")?;
                    let extension: String = matches
                        .value_of_t("extension")
                        .context("Could not parse `extension` as a String.")?;

                    let window_size: usize = matches
                        .value_of_t("window")
                        .context("Could not parse `window-size` as usize.")?;
                    let outdir = matches
                        .value_of("dir")
                        .context("Could not get the value of `dir`")?;
                    let output = matches
                        .value_of("output")
                        .context("Could not get the value of `output`")?;

                    // create file
                    let file_name = format!(
                        "{}/{}{}{}",
                        outdir, output, "_telomeric_repeat_windows.", extension
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
                        input_fasta,
                        telomeric_repeat,
                        window_size
                    );

                    // create file
                    let log_file_name = format!("{}/{}{}", outdir, output, ".log");
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
