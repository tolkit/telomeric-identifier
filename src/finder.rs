use crate::clades::TelomereSeq;
use crate::{clades, utils, SubCommand};
use anyhow::{Context, Result};
use bio::io::fasta;
use std::fs::{create_dir_all, File};
use std::io::LineWriter;
use std::io::Write;
use std::path::PathBuf;
use std::process;
use std::str;

/// The entry point for `tidk find`.
///
/// Finder uses the clade specific telomere sequence and queries against the genome.
pub fn finder(matches: &clap::ArgMatches, sc: SubCommand) -> Result<()> {
    // print table of telomeric sequences
    if matches.get_flag("print") {
        clades::print_table();
        process::exit(1);
    }

    let input_fasta: PathBuf = matches
        .get_one::<PathBuf>("fasta")
        .expect("errored by clap")
        .clone();
    let reader = fasta::Reader::from_file(input_fasta)?;

    let clade = matches.get_one::<String>("clade").expect("errored by clap");
    let clade_info: TelomereSeq = clade.parse()?;
    match clade_info.seq.len() {
        0 => {}
        1 => eprintln!(
            "[+]\tSearching genome for a single telomeric repeat: {}",
            clade_info
                .seq
                .get(0)
                .context("Could not get the first element of `seq`.")?
        ),
        length => {
            eprintln!("[+]\tSearching genome for {} telomeric repeats:", length);
            for telomeric_repeat in 0..length {
                eprintln!(
                    "[+]\t\t{}",
                    clade_info.seq.get(telomeric_repeat).context(format!(
                        "Could not get the {} element of `seq`.",
                        telomeric_repeat
                    ))?
                );
            }
        }
    }

    let window_size: usize = *matches.get_one::<usize>("window").expect("errored by clap");
    let outdir = matches.get_one::<PathBuf>("dir").expect("errored by clap");
    let output = matches
        .get_one::<PathBuf>("output")
        .expect("errored by clap");

    // create directory for output
    create_dir_all(outdir)?;

    // create file
    let file_name = format!(
        "{}/{}{}",
        outdir.display(),
        output.display(),
        "_telomeric_repeat_windows.tsv"
    );
    let finder_file = File::create(file_name)?;
    let mut finder_file = LineWriter::new(finder_file);
    // add headers
    writeln!(
        finder_file,
        "id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat"
    )?;

    // extract the string from TelomereSeq struct
    // dereference here because of Box<T>
    let telomeric_repeat = clade_info.seq.inner;

    // iterate over the fasta records
    for result in reader.records() {
        let record = result?;
        let id = record.id().to_owned();

        // fn window counter
        write_window_counts(
            record,
            &mut finder_file,
            clade_info.clone(),
            telomeric_repeat,
            window_size,
            id.clone(),
        )?;

        eprintln!("[+]\tChromosome {} processed", id);
    }
    eprintln!("[+]\tFinished searching genome.");

    // optional log file
    sc.log(matches)?;

    Ok(())
}

/// Creates the window iterator and iterates over each iteration of the
/// fasta file, writing on the fly.
fn write_window_counts<T: std::io::Write>(
    sequence: bio::io::fasta::Record,
    file: &mut LineWriter<T>,
    clade_info: clades::TelomereSeq,
    telomeric_repeat: &[&str],
    window_size: usize,
    id: String,
) -> Result<()> {
    // needed as in some clades there is more than one telomeric repeat sequence
    let mut telomeric_repeat_index = 0;
    loop {
        // break this loop if we reach the end of &[&str] of telomeric repeats
        if clade_info.seq.len() == telomeric_repeat_index {
            break;
        }

        // get forward and reverse sequences, and length
        // to remove overlapping matches.
        let forward_telomeric_seq =
            *telomeric_repeat
                .get(telomeric_repeat_index)
                .context(format!(
                    "Could not get the telomeric repeat with index: {}.",
                    telomeric_repeat_index
                ))?;
        let reverse_telomeric_seq = utils::reverse_complement(forward_telomeric_seq);
        let current_telomeric_length = forward_telomeric_seq.len();

        // create the iterator in each loop iteration isnt costly is it?
        let windows = sequence.seq().chunks(window_size);
        // keep track of the window size

        let mut end = window_size;
        // iterate over windows
        for (i, window) in windows.enumerate() {
            // make window uppercase
            let windows_upper = str::from_utf8(window)?.to_uppercase();
            // for each window, find the motifs in this
            let forward_motif = utils::find_motifs(forward_telomeric_seq, &windows_upper);
            let reverse_motif = utils::find_motifs(&reverse_telomeric_seq, &windows_upper);

            // remove overlapping matches
            // not sure this is necessary, but thought it might be...
            let forward_motif_noverlap =
                utils::remove_overlapping_indexes(forward_motif, current_telomeric_length);
            let reverse_motif_noverlap =
                utils::remove_overlapping_indexes(reverse_motif, current_telomeric_length);

            // the number of matches for forward/reverse
            let forward_repeat_number = forward_motif_noverlap.len();
            let reverse_repeat_number = reverse_motif_noverlap.len();

            if i != 0 {
                end += window_size;
            }
            end = std::cmp::min(end, sequence.seq().len());
            // write to file
            writeln!(
                file,
                "{}\t{}\t{}\t{}\t{}",
                id, end, forward_repeat_number, reverse_repeat_number, forward_telomeric_seq
            )?;
        }
        // go to the next telomeric repeat (if there is one)
        telomeric_repeat_index += 1;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::io::{LineWriter, Read};

    use crate::clades::{Seq, TelomereSeq};

    use super::write_window_counts;

    // a wrapper for making a bio::io::fasta record
    fn make_record(id: &str, seq: &[u8]) -> bio::io::fasta::Record {
        bio::io::fasta::Record::with_attrs(id, None, seq)
    }

    // take a record, write to a vector (fake file), then read out of this the output.
    fn calc_windows(rec: bio::io::fasta::Record, ts: TelomereSeq, ws: usize) -> String {
        let file = Vec::new();
        let mut lw = LineWriter::new(file);
        let id = rec.id().to_owned();

        let telomeric_repeat = ts.seq.inner;
        write_window_counts(rec, &mut lw, ts, telomeric_repeat, ws, id).unwrap();

        // read file contents to new vec
        let mut out = Vec::new();
        // we created a line writer, so we need the underlying writer
        let c = lw.into_inner().unwrap();
        let mut d = c.as_slice();
        d.read_to_end(&mut out).unwrap();

        String::from_utf8(out).unwrap()
    }

    #[test]
    fn test_search_1() {
        let rec = make_record(
            "test1",
            b"AAACCCTAAACCCTAAACCCTTGAGAGAGGGGGTGTGGGGAGGGGTTGAGAAACCCT",
        );

        let apiales = TelomereSeq {
            clade: crate::clades::Clade::Apiales,
            seq: Seq {
                inner: &["AAACCCT"],
            },
        };

        let windows_calculation = calc_windows(rec, apiales, 20);

        let rows: Vec<&str> = windows_calculation.lines().collect();

        // three in first window
        assert_eq!(rows[0], "test1\t20\t2\t0\tAAACCCT");
        // none in second
        assert_eq!(rows[1], "test1\t40\t0\t0\tAAACCCT");
        // two in third
        assert_eq!(rows[2], "test1\t57\t1\t0\tAAACCCT");
    }
}
