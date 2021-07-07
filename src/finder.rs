pub mod finder {
    use crate::clades::clades;
    use crate::utils::utils;
    use bio::io::fasta;
    use clap::value_t;
    use std::fs::{create_dir_all, File};
    use std::io::LineWriter;
    use std::io::Write;
    use std::process;
    use std::str;

    // finder uses the clade specific telomere sequence and queries against the genome.

    pub fn finder(matches: &clap::ArgMatches) {
        // print table of telomeric sequences
        if matches.is_present("print") {
            clades::print_table();
            process::exit(1);
        }

        let input_fasta = matches.value_of("fasta").unwrap();
        let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

        let clade = value_t!(matches.value_of("clade"), String).unwrap_or_else(|e| e.exit());
        let clade_info = clades::return_telomere_sequence(&clade);

        if clade_info.length == 1 {
            println!(
                "[+]\tSearching genome for a single telomeric repeat: {}",
                clade_info.seq.get(0).unwrap().to_owned()
            );
        } else if clade_info.length > 1 {
            println!(
                "[+]\tSearching genome for {} telomeric repeats:",
                clade_info.length
            );
            for telomeric_repeat in 0..clade_info.length {
                println!("[+]\t\t{}", clade_info.seq.get(telomeric_repeat).unwrap());
            }
        }

        let window_size = value_t!(matches.value_of("window"), usize).unwrap_or_else(|e| e.exit());
        let output = matches.value_of("output").unwrap();

        // create directory for output
        if let Err(e) = create_dir_all("./finder/") {
            println!("[-]\tCreate directory error: {}", e.to_string());
        }
        // create file
        let file_name = format!("./finder/{}{}", output, "_telomeric_repeat_windows.csv");
        let finder_file = File::create(&file_name).unwrap();
        let mut finder_file = LineWriter::new(finder_file);
        // add headers
        writeln!(
            finder_file,
            "id,window,forward_repeat_number,reverse_repeat_number,telomeric_repeat"
        )
        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));

        // extract the string from TelomereSeq struct
        // dereference here because of Box<T>
        let telomeric_repeat = *clade_info.seq;

        // iterate over the fasta records
        for result in reader.records() {
            let record = result.expect("[-]\tError during fasta record parsing.");
            let id = record.id().to_owned();

            // fn window counter
            write_window_counts(
                record,
                &mut finder_file,
                clade_info.clone(),
                telomeric_repeat,
                window_size,
                id.clone(),
            )
            .expect("Could not write to file.");

            println!("[+]\tChromosome {} processed", id);
        }
        println!("[+]\tFinished searching genome.");
    }

    // creates the window iterator and iterates over each iteration of the
    // fasta file, writing on the fly.

    fn write_window_counts<T: std::io::Write>(
        sequence: bio::io::fasta::Record,
        file: &mut LineWriter<T>,
        clade_info: clades::TelomereSeq,
        telomeric_repeat: &[&str],
        window_size: usize,
        id: String,
    ) -> std::io::Result<()> {
        // needed as in some clades there is more than one telomeric repeat sequence
        let mut telomeric_repeat_index = 0;
        loop {
            // break this loop if we reach the end of &[&str] of telomeric repeats
            if clade_info.length == telomeric_repeat_index {
                break;
            }

            // get forward and reverse sequences, and length
            // to remove overlapping matches.
            let forward_telomeric_seq = *telomeric_repeat.get(telomeric_repeat_index).unwrap();
            let reverse_telomeric_seq = utils::reverse_complement(forward_telomeric_seq);
            let current_telomeric_length = forward_telomeric_seq.len();

            // create the iterator in each loop iteration isnt costly is it?
            let windows = sequence.seq().chunks(window_size);
            // keep track of the window size
            let mut window_index = window_size;
            // iterate over windows
            for window in windows {
                // make window uppercase
                let windows_upper = str::from_utf8(window).unwrap().to_uppercase();
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
                // write to file
                writeln!(
                    file,
                    "{},{},{},{},{}",
                    id,
                    window_index,
                    forward_repeat_number,
                    reverse_repeat_number,
                    forward_telomeric_seq
                )
                .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
                // increment window
                window_index += window_size;
            }
            // go to the next telomeric repeat (if there is one)
            telomeric_repeat_index += 1;
        }
        Ok(())
    }
}
