// a reduced version of search
use crate::utils;
use bio::io::fasta;
use std::fs::{create_dir_all, File};
use std::io::LineWriter;
use std::io::Write;
use std::str;

// uses a user defined string to trim reads

pub fn trim(matches: &clap::ArgMatches) {
    let input_fasta = matches.value_of("fasta").unwrap();
    let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

    let min_occur = matches.value_of_t("min_occur").unwrap_or_else(|e| e.exit());
    let min_len = matches.value_of_t("min_len").unwrap_or_else(|e| e.exit());

    let telomeric_repeat: String = matches.value_of_t("string").unwrap_or_else(|e| e.exit());
    let reverse_telomeric_seq = utils::reverse_complement(&telomeric_repeat);
    // min_occur * telomeric repeats
    let multiple_telomeric_repeat = telomeric_repeat.repeat(min_occur);
    let multiple_reverse_repeat = reverse_telomeric_seq.repeat(min_occur);
    let telomeric_length = telomeric_repeat.len();

    eprintln!(
        "[+]\tSearching genome for telomeric repeat: {}",
        telomeric_repeat
    );

    // create directory for output
    if let Err(e) = create_dir_all("./trim/") {
        eprintln!("[-]\tCreate directory error: {}", e.to_string());
    }

    let output = matches.value_of("output").unwrap();

    // create file
    let file_name = format!("./trim/{}{}", output, "_trimmed.fasta");

    let search_file = File::create(&file_name).unwrap();
    let mut search_file = LineWriter::new(search_file);
    let mut num_trimmed_reads = 0;

    // iterate over the fasta records
    for result in reader.records() {
        let record = result.expect("[-]\tError during fasta record parsing.");
        let id = record.id().to_owned();
        // check if the start of the read matches the reverse complemented telomeric repeat
        let matches_start = str::from_utf8(&record.seq()[..telomeric_length * 3])
            .unwrap()
            .find(&reverse_telomeric_seq);
        // check if the end of the read matches the telomeric repeat
        let matches_end =
            str::from_utf8(&record.seq()[record.seq().len() - telomeric_length * 3..])
                .unwrap()
                .find(&telomeric_repeat);

        if !matches_end.is_none() {
            let telo_pos = str::from_utf8(record.seq())
                .unwrap()
                .find(&multiple_telomeric_repeat);

            // catch the none
            if telo_pos.is_none() {
                eprintln!(
                    "[-]\tAt `matches_end` for sequence ID {}: no multiple telomeric repeat found.",
                    id
                );
                continue;
            }
            if telo_pos.unwrap() < min_len {
                continue;
            }
            let trimmed_seq = utils::reverse_complement(
                &str::from_utf8(&record.seq()[..telo_pos.unwrap()]).unwrap(),
            );

            if let Err(e) = writeln!(search_file, ">{}\n{}", id, trimmed_seq) {
                eprintln!("[-]\tCouldn't write trimmed reads: {}", e.to_string());
            }
            num_trimmed_reads += 1;
        }

        if !matches_start.is_none() {
            let telo_pos = str::from_utf8(record.seq())
                .unwrap()
                .rfind(&multiple_reverse_repeat);
            // catch the none
            if telo_pos.is_none() {
                eprintln!("[-]\tAt `matches_start` for sequence ID {}: no multiple telomeric repeat found.", id);
                continue;
            }

            if record.seq().len() - telo_pos.unwrap() < min_len {
                continue;
            }
            let trimmed_seq =
                str::from_utf8(&record.seq()[telo_pos.unwrap() + telomeric_length * min_occur..])
                    .unwrap();

            if let Err(e) = writeln!(search_file, ">{}\n{}", id, trimmed_seq) {
                eprintln!("[-]\tCouldn't write trimmed reads: {}", e.to_string());
            }
            num_trimmed_reads += 1;
        }
    }
    eprintln!(
        "[+]\tWrote {} reads longer than {} nucleotides after trimming {} repeat.",
        num_trimmed_reads, min_len, telomeric_repeat
    );
}