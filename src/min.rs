pub mod min {
    use crate::utils::utils;
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    pub fn min_dna_string(matches: &clap::ArgMatches) {
        let input_file = matches.value_of("file");

        // if there is a file command line arg
        match input_file {
            Some(file) => {
                let file = File::open(file).expect("Could not open file.");
                let reader = BufReader::new(file);
                for line in reader.lines() {
                    let line = line.expect("Could not read line.");
                    let res = utils::lex_min(&line);
                    println!("{}", res);
                }
            }
            None => {
                if let Some(i) = matches.values_of("DNA string") {
                    for el in i {
                        let res = utils::lex_min(el);
                        println!("{}", res);
                    }
                }
            }
        }
    }
}
