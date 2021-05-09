pub mod trim {

    // a reduced version of search
    use crate::utils::utils;
    use bio::io::fasta;
    use clap::value_t;
    use std::fs::File;
    use std::io::LineWriter;
    use std::io::Write;
    use std::str;

    // uses a user defined string to trim reads

    pub fn trim(matches: &clap::ArgMatches) {
        let input_fasta = matches.value_of("fasta").unwrap();
        let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

        let min_occur = value_t!(matches.value_of("min_occur"), usize).unwrap_or_else(|e| e.exit());
        let min_len = value_t!(matches.value_of("min_len"), usize).unwrap_or_else(|e| e.exit());

        let telomeric_repeat =
            value_t!(matches.value_of("string"), String).unwrap_or_else(|e| e.exit());
        let reverse_telomeric_seq = utils::reverse_complement(&telomeric_repeat);
        // min_occur * telomeric repeats
        let multiple_telomeric_repeat = telomeric_repeat.repeat(min_occur);
        let multiple_reverse_repeat = reverse_telomeric_seq.repeat(min_occur);
        let telomeric_length = telomeric_repeat.len();


        println!(
            "[+]\tSearching genome for telomeric repeat: {}",
            telomeric_repeat
        );

        let file_name = matches.value_of("output").unwrap();

        let search_file = File::create(&file_name).unwrap();
        let mut search_file = LineWriter::new(search_file);
        
        let mut num_trimmed_reads = 0;

        

        // iterate over the fasta records
        for result in reader.records() {
            let record = result.expect("[-]\tError during fasta record parsing.");
            let id = record.id().to_owned();
            // check if the start of the read matches the reverse complemented telomeric repeat 
            let matches_start = str::from_utf8(&record.seq()[.. telomeric_length * 3]).unwrap().find(&reverse_telomeric_seq);
            // check if the end of the read matches the telomeric repeat 
            let matches_end = str::from_utf8(&record.seq()[record.seq().len() - telomeric_length * 3 ..]).unwrap().find(&telomeric_repeat);

            if !matches_end.is_none() {
                let telo_pos = str::from_utf8(record.seq()).unwrap().find(&multiple_telomeric_repeat).unwrap();
                if telo_pos < min_len {
                    continue;
                }
                let trimmed_seq = utils::reverse_complement(&str::from_utf8(&record.seq()[.. telo_pos]).unwrap());
                // println!("[+]\tLocation of telomere in read {} is {}", id, telo_pos);
                // println!("[+]\tSubread is {}", trimmed_seq);
                if let Err(e) = writeln!(search_file, ">{}\n{}", id, trimmed_seq) {
                    println!("Couldn't write trimmed reads: {}", e.to_string());   
                }
                num_trimmed_reads+=1;
            }

            if !matches_start.is_none() {
                let telo_pos = str::from_utf8(record.seq()).unwrap().rfind(&multiple_reverse_repeat).unwrap();
                if record.seq().len() - telo_pos < min_len {
                    continue;
                }
                let trimmed_seq = str::from_utf8(&record.seq()[telo_pos + telomeric_length * min_occur ..]).unwrap();
                // println!("[+]\tLocation of telomere in read {} is {}", id, telo_pos);
                // println!("[+]\tSubread is {}", trimmed_seq);
                if let Err(e) = writeln!(search_file, ">{}\n{}", id, trimmed_seq) {
                    println!("Couldn't write trimmed reads: {}", e.to_string());   
                }
                num_trimmed_reads+=1;
            }

            


        }
        println!("[+]\t Wrote {} reads longer than {} nts after trimming {} repeat.", num_trimmed_reads, min_len, telomeric_repeat);
    }
}
