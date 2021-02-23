pub mod explore {
    use bio::io::fasta;
    use clap::value_t;
    use std::str;

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

            // store the index positions and adjacent equivalent sequences.
            let mut indexes = Vec::new();
            // need this otherwise we lose the position in the sequence.
            // need to check this is actually correct.
            let mut pos = 1;

            // this is the heavy lifting.
            for (i, (a, b)) in chunks.zip(chunks_plus_one).enumerate() {
                // if chunk contains N, skip.
                if a.contains(&78) || b.contains(&78) {
                    pos += length;
                    continue;
                } else if a == b {
                    indexes.push((pos, str::from_utf8(a).unwrap()));
                }
                pos += length;
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
                pos += length;
            }

            // collect all of this vector into a nice format to print for the moment.
            // iteration
            let mut it = 1;
            let mut count = 0;
            // threshold
            let threshold = 10;
            let mut start_pos = 1;
            let mut start;
            let mut end;

            loop {
                start = adjacent_indexes[start_pos].0;
                if it == adjacent_indexes.len() - 1 {
                    break;
                }
                if adjacent_indexes[it].1 == length {
                    count += 1;
                    it += 1;
                } else {
                    end = adjacent_indexes[it].0;
                    if count > threshold {
                        println!(
                            "{:?} - {:?}: {:?} x repeat sequence: {:?}",
                            start, end, count, adjacent_indexes[it].2
                        );
                    }

                    start_pos = it;
                    count = 1;
                    it += 1;
                }
            }
        }
    }
}
