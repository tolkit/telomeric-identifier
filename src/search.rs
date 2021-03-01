pub mod search {

    // essentially the same as finder
    use crate::utils::utils;
    use bio::io::fasta;
    use clap::value_t;
    use std::fs::{create_dir_all, File};
    use std::io::LineWriter;
    use std::io::Write;
    use std::str;

    // essentially same as finder

    pub fn search(matches: &clap::ArgMatches) {
        let input_fasta = matches.value_of("fasta").unwrap();
        let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

        let telomeric_repeat =
            value_t!(matches.value_of("string"), String).unwrap_or_else(|e| e.exit());

        let window_size = value_t!(matches.value_of("window"), usize).unwrap_or_else(|e| e.exit());
        let output = matches.value_of("output").unwrap();

        // create directory for output
        if let Err(e) = create_dir_all("./search/") {
            println!("[-]\tCreate directory error: {}", e.to_string());
        }
        // create file
        let file_name = format!("./search/{}{}", output, "_telomeric_repeat_windows.csv");
        let search_file = File::create(&file_name).unwrap();
        let mut search_file = LineWriter::new(search_file);
        // add headers
        writeln!(
            search_file,
            "ID,window,forward_repeat_number,reverse_repeat_number,telomeric_repeat"
        )
        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));

        // iterate over the fasta records
        for result in reader.records() {
            let record = result.expect("[-]\tError during fasta record parsing.");
            let id = record.id().to_owned();

            // fn window counter
            write_window_counts(record, &mut search_file, &telomeric_repeat, window_size, id);
        }
    }

    // write windows to file
    fn write_window_counts<T: Write>(
        sequence: bio::io::fasta::Record,
        file: &mut LineWriter<T>,
        telomeric_repeat: &str,
        window_size: usize,
        id: String,
    ) {
        // get forward and reverse sequences, and length
        // to remove overlapping matches.
        let forward_telomeric_seq = telomeric_repeat;
        let reverse_telomeric_seq = utils::reverse_complement(forward_telomeric_seq);
        let telomeric_length = forward_telomeric_seq.len();

        // create the iterator in each loop iteration isnt costly is it?
        let windows = sequence.seq().chunks(window_size);
        // keep track of the window size
        let mut window_index = window_size;
        // iterate over windows
        for window in windows {
            // for each window, find the motifs in this
            let forward_motif =
                utils::find_motifs(forward_telomeric_seq, str::from_utf8(window).unwrap());
            let reverse_motif =
                utils::find_motifs(&reverse_telomeric_seq, str::from_utf8(window).unwrap());

            // remove overlapping matches
            // not sure this is necessary, but thought it might be...
            let forward_motif_noverlap =
                utils::remove_overlapping_indexes(forward_motif, telomeric_length);
            let reverse_motif_noverlap =
                utils::remove_overlapping_indexes(reverse_motif, telomeric_length);

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
    }
}
