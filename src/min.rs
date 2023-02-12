use crate::utils;
use anyhow::{bail, Result};
use bio::io::fasta;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::PathBuf;
use std::str;

// TODO:
// the CLI API is pretty bad currently...
// would be better to detect STDIN with atty.

/// The entry point for `tidk min`.
pub fn min_dna_string(matches: &clap::ArgMatches) -> Result<()> {
    let input_file = matches.get_one::<PathBuf>("file");
    let is_fasta = matches.get_flag("fasta");

    // if there is a file command line arg
    match input_file {
        Some(file) => {
            let fasta = file.ends_with(".fa") || file.ends_with(".fasta");
            if fasta {
                let reader = fasta::Reader::from_file(file)?;
                for record in reader.records() {
                    let record = record?;
                    let seq = str::from_utf8(record.seq())?;
                    let res = utils::lex_min(seq);
                    // write to fasta
                    let mut writer = fasta::Writer::new(io::stdout());
                    writer.write(record.id(), Some("tidk-min"), res.as_bytes())?;
                }
                Ok(())
            } else {
                let file = File::open(file)?;
                let reader = BufReader::new(file);
                for line in reader.lines() {
                    let line = line?;
                    let res = utils::lex_min(&line);
                    println!("{}", res);
                }
                Ok(())
            }
        }
        None => {
            let matches = matches.try_get_many("dna_string");

            if let Ok(Some(m)) = matches {
                let m: Vec<&String> = m.clone().collect();
                for el in m {
                    let res = utils::lex_min(el);
                    println!("{}", res);
                }
                Ok(())
            } else if is_fasta {
                let mut records = fasta::Reader::new(io::stdin()).records();
                while let Some(Ok(record)) = records.next() {
                    let seq = str::from_utf8(record.seq())?;
                    let res = utils::lex_min(seq);
                    // write to fasta
                    let mut writer = fasta::Writer::new(io::stdout());
                    writer.write(record.id(), Some("tidk-min"), res.as_bytes())?;
                }
                Ok(())
            } else {
                let stdin = io::stdin();
                for line in stdin.lock().lines() {
                    let line = line?;
                    if line.chars().next() == ">".chars().next() {
                        bail!("Use -x option, for fasta input.")
                    }
                    let res = utils::lex_min(&line);
                    println!("{}", res);
                }
                Ok(())
            }
        }
    }
}
