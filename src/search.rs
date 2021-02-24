pub mod search {
    // search the genome with a string
    // needs forward and reverse complement of the string
    // what's the output? Where should it look?

    // for simplicity it can look everywhere
    // in windows of defined size?

    use crate::utils::utils;
    use bio::io::fasta;
    use clap::value_t;
    use std::str;

    pub fn search(matches: &clap::ArgMatches) {
        let input_fasta = matches.value_of("fasta").unwrap();
        let search_string =
            value_t!(matches.value_of("string"), String).unwrap_or_else(|e| e.exit());

        let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

        for result in reader.records() {
            let record = result.expect("[-]\tError during fasta record parsing.");
            let id = record.id().to_owned();

            let x = utils::find_motifs(&search_string, str::from_utf8(record.seq()).unwrap());
            let test = utils::longest_repeat(x);
            println!("{:?}", test);
        }
    }
}
