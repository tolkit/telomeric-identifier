pub mod search {
    // search the genome with a string
    // needs forward and reverse complement of the string
    // what's the output? Where should it look?

    // for simplicity it can look everywhere
    // in windows of defined size?

    use bio::io::fasta;
    use clap::value_t;
    pub fn search(matches: &clap::ArgMatches) {
        let input_fasta = matches.value_of("fasta").unwrap();
        let search_string =
            value_t!(matches.value_of("string"), String).unwrap_or_else(|e| e.exit());
        println!(
            "Search string is {}. Nothing else implemented.",
            search_string
        );
        let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

        // for record in reader.records() {}
    }
}
