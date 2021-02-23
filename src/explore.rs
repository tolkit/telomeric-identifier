pub mod explore {
    // I think this is a stupid idea.
    use crate::patterns::patterns;
    use crate::utils::utils;
    use bio::io::fasta;
    use clap::value_t;
    use std::collections::HashMap;
    use std::str;

    // pub fn explore(matches: &clap::ArgMatches) {
    //     let input_fasta = matches.value_of("fasta").unwrap();
    //     let length = value_t!(matches.value_of("length"), usize).unwrap_or_else(|e| e.exit());

    //     // read in the fasta from file
    //     let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

    //     let test = patterns::enumerate_mers(length);

    //     for result in reader.records() {
    //         let record = result.expect("[-]\tError during fasta record parsing.");
    //         let seq = str::from_utf8(record.seq()).unwrap();

    //         let mut kmer_value: HashMap<&str, usize> = HashMap::new();
    //         for i in test.iter() {
    //             let telomeric_repeats = utils::find_motifs(i, seq);
    //             let mut m: HashMap<i32, usize> = HashMap::new();
    //             let repeats = utils::longest_repeat(telomeric_repeats);
    //             for x in repeats {
    //                 *m.entry(x as i32).or_default() += 1;
    //             }
    //             let max = m.get(&(length as i32)).cloned();

    //             kmer_value.insert(&i, max.unwrap_or(0));
    //         }
    //         let map = kmer_value
    //             .iter()
    //             .max_by(|a, b| a.1.cmp(&b.1))
    //             .map(|(k, v)| (k, v));

    //         println!("{:?}", map);
    //     }
    // }

    // maybe a better idea to split the fasta sequence into chunks
    // while the ith and i + 1 chunk are the same, push tuple to a vector
    // if this count is less than x discard

    // (index, sequence)

    // then consecutive sequences can be counted at their index position
    // do we want a threshold? 5?
    // 200 - 235: 5 x repeat sequence AAATTAA

    pub fn explore(matches: &clap::ArgMatches) {
        let input_fasta = matches.value_of("fasta").unwrap();
        let length = value_t!(matches.value_of("length"), usize).unwrap_or_else(|e| e.exit());

        let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

        for result in reader.records() {
            let record = result.expect("[-]\tError during fasta record parsing.");
            let chunks = record.seq().chunks(length);
            let chunks_plus_one = record.seq()[length..record.seq().len()].chunks(length);

            // now we have consecutive records a and b
            let mut indexes = Vec::new();
            for (i, (a, b)) in chunks.zip(chunks_plus_one).enumerate() {
                // if chunk contains N, skip.
                if a.contains(&78) || b.contains(&78) {
                    continue;
                } else if a == b {
                    indexes.push((i, str::from_utf8(a).unwrap()));
                }
            }
            let mut adjacent_indexes = Vec::new();
            for i in 1..indexes.len() {
                if i + 1 < indexes.len() {
                    adjacent_indexes.push((
                        indexes[i].0,
                        indexes[i + 1].0 - indexes[i].0,
                        indexes[i].1,
                    ))
                }
            }
            for i in adjacent_indexes.into_iter() {
                println!("{:?}", i);
            }
        }
    }
}
