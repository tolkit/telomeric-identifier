pub mod explore {
    use crate::utils::utils;
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
            let id = record.id().to_owned();

            let indexes = chunk_fasta(record, length);
            let adjacents = calculate_indexes(indexes);
            let formatted = generate_explore_data(adjacents, length);

            for result in formatted {
                println!(
                    "Chromosome {}: {:?} - {:?}: {:?} x repeat sequence: {:?}",
                    id, result.start, result.end, result.count, result.sequence
                );
            }
        }
    }

    // have to return a String. Otherwise Rust whines.
    fn chunk_fasta(sequence: bio::io::fasta::Record, chunk_length: usize) -> Vec<(usize, String)> {
        let chunks = sequence.seq().chunks(chunk_length);
        let chunks_plus_one =
            sequence.seq()[chunk_length..sequence.seq().len()].chunks(chunk_length);

        // store the index positions and adjacent equivalent sequences.
        let mut indexes = Vec::new();
        // need this otherwise we lose the position in the sequence.
        // need to check this is actually correct.
        let mut pos = 1;

        // this is the heavy lifting.
        for (_i, (a, b)) in chunks.zip(chunks_plus_one).enumerate() {
            // if chunk contains N, skip.
            if a.contains(&78) || b.contains(&78) {
                pos += chunk_length;
                continue;
            } else if a == b {
                indexes.push((pos, str::from_utf8(a).unwrap().to_owned()));
            }
            pos += chunk_length;
        }
        indexes
    }

    fn calculate_indexes(indexes: Vec<(usize, String)>) -> Vec<(usize, usize, String)> {
        let mut adjacent_indexes = Vec::new();
        for i in 1..indexes.len() {
            if i + 1 < indexes.len() {
                adjacent_indexes.push((
                    indexes[i].0,
                    indexes[i + 1].0 - indexes[i].0,
                    indexes[i].1.clone(), //otherwise Rust complains
                ))
            }
        }
        adjacent_indexes
    }

    pub struct TelomericRepeatExplore {
        pub start: usize,
        pub end: usize,
        pub count: i32,
        pub sequence: String,
    }

    fn generate_explore_data(
        adjacent_indexes: Vec<(usize, usize, String)>,
        chunk_length: usize,
    ) -> Vec<TelomericRepeatExplore> {
        // collect all of this vector into a nice format to print for the moment.
        // iteration
        let mut it = 1;
        let mut count = 0;
        // threshold could become an input parameter
        let threshold = 10;
        let mut start_pos = 1;
        let mut start;
        let mut end;
        let mut potential_telomeric_repeats = Vec::new();

        loop {
            start = adjacent_indexes[start_pos].0;
            if it == adjacent_indexes.len() - 1 {
                break;
            }
            if adjacent_indexes[it].1 == chunk_length
            //&& utils::string_rotation(adjacent_indexes[it].2, adjacent_indexes[it - 1].2)
            {
                count += 1;
                it += 1;
            } else {
                end = adjacent_indexes[it].0;
                if count > threshold {
                    potential_telomeric_repeats.push(TelomericRepeatExplore {
                        start: start,
                        end: end,
                        count: count,
                        sequence: adjacent_indexes[it].2.clone(),
                    });
                }

                start_pos = it;
                count = 1;
                it += 1;
            }
        }
        potential_telomeric_repeats
    }

    // bit experimental, not sure it totally works...
    // fn merge_rotated_repeats(data: Vec<(usize, usize, i32, std::string::String)>) {
    //     let mut it = 1;
    //     let mut start_pos = 1;
    //     let mut start;
    //     let mut count = 0;
    //     let mut end;

    //     loop {
    //         start = data[start_pos].0;
    //         if it == data.len() - 2 {
    //             break;
    //         }
    //         // if adjacent sequences are rotations
    //         if utils::string_rotation(&data[it].3, &data[it + 1].3) {
    //             count += data[it].2;
    //             it += 1;
    //         } else {
    //             end = data[it].0;
    //             if count > 0 {
    //                 println!("{:?}, {:?}, {:?}, {:?}", start, end, count, &data[it].3);
    //             }
    //             count = 0;
    //         }
    //         start_pos = it;
    //         it += 1;
    //     }
    // }
}
