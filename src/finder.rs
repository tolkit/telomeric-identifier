pub mod finder {

    use crate::patterns::patterns;
    use bio::io::fasta;
    use clap::value_t;

    pub fn finder(matches: &clap::ArgMatches) {
        let clade = value_t!(matches.value_of("clade"), String).unwrap_or_else(|e| e.exit());
        let tostr = clade.to_owned();
        let s_slice: &str = &tostr[..];
        let capture = patterns::return_telomere_sequence(s_slice);

        println!("{:?}", capture.seq);
    }
}
