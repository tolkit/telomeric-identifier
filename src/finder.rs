pub mod finder {

    use crate::patterns::patterns;
    use bio::io::fasta;
    use clap::value_t;

    // finder uses the clade specific telomere sequence and queries against the genome.
    pub fn finder(matches: &clap::ArgMatches) {
        let input_fasta = matches.value_of("fasta").unwrap();
        let clade = value_t!(matches.value_of("clade"), String).unwrap_or_else(|e| e.exit());
        let tostr = clade.to_owned();
        let s_slice: &str = &tostr[..];
        let capture = patterns::return_telomere_sequence(s_slice);

        let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");
        // do stuff...
        println!("Nothing implemented yet {:?}", capture.seq);
    }
}
