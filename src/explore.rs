pub mod explore {
    use crate::utils::utils;
    use bio::io::fasta;
    use clap::value_t;
    use std::fs::{create_dir_all, File};
    use std::io::prelude::*;
    use std::io::LineWriter;
    use std::str;

    // split the fasta sequence into chunks
    // while the ith and i + 1 chunk are the same, push tuple to a vector
    // if this count is less than x discard

    // the main useful data structure is:
    // (start, end, count, sequence)

    // the entry function called from main
    // iterates over the fasta file and calls the functions below.
    // I think it could be in general condensed, there seems to be some redundancy.

    // this should identify the telomeric repeat as reverse complement at either
    // end of a chromosome
    // I am not sure this will be as effective a tool on raw reads, but this remains
    // to be tested.

    pub fn explore(matches: &clap::ArgMatches) {
        // parse arguments from main
        let input_fasta = matches.value_of("fasta").unwrap();
        let length = value_t!(matches.value_of("length"), usize).unwrap_or_else(|e| e.exit());

        // if length is not set, these are the lengths (and length itself is set to zero)
        let minimum = value_t!(matches.value_of("minimum"), usize).unwrap_or_else(|e| e.exit());
        let maximum = value_t!(matches.value_of("maximum"), usize).unwrap_or_else(|e| e.exit());

        let threshold = value_t!(matches.value_of("threshold"), i32).unwrap_or_else(|e| e.exit());
        let output = matches.value_of("output").unwrap();

        let dist_from_chromosome_end =
            value_t!(matches.value_of("distance"), usize).unwrap_or_else(|e| e.exit());

        // create directory for output
        if let Err(e) = create_dir_all("./explore/") {
            println!("[-]\tCreate directory error: {}", e.to_string());
        }

        // create file
        let file_name = format!("./explore/{}{}", output, "_telomeric_locations.csv");
        let explore_file = File::create(&file_name).unwrap();
        let mut explore_file = LineWriter::new(explore_file);

        // add headers
        writeln!(
            explore_file,
            "id,start_pos,end_pos,repeat_number,repeat_sequence,sequence_length"
        )
        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));

        // to report the telomeres...
        let mut output_vec: Vec<FormatTelomericRepeat> = Vec::new();
        // i.e. if you chose a length, as opposed to a minmum/maximum
        if length > 0 {
            println!(
                "[+]\tExploring genome for potential telomeric repeats of length: {}",
                length
            );
            let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

            for result in reader.records() {
                let record = result.expect("[-]\tError during fasta record parsing.");
                let id = record.id().to_owned();
                let seq_len = record.seq().len();

                let indexes = chunk_fasta(record, length);
                let adjacents = calculate_indexes(indexes);
                let formatted = generate_explore_data(adjacents, id.clone(), length);
                output_vec.append(
                    &mut merge_rotated_repeats(
                        formatted.unwrap_or(vec![]),
                        length,
                        id.clone(),
                        threshold,
                        &mut explore_file,
                        seq_len,
                        dist_from_chromosome_end,
                    )
                    .unwrap_or(vec![]),
                );

                println!("[+]\tChromosome {} processed", id);
            }
        } else {
            // how can I parallelise this??
            println!(
                "[+]\tExploring genome for potential telomeric repeats between lengths {} and {}.",
                minimum, maximum
            );
            for length in minimum..maximum + 1 {
                println!("[+]\t\tFinding telomeric repeat length: {}", length);

                // have to call reader in the loop, as otherwise `reader` doesn't live long enough.
                // I expect it's not an expensive call anyway.
                let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

                for result in reader.records() {
                    let record = result.expect("[-]\tError during fasta record parsing.");
                    let id = record.id().to_owned();
                    let seq_len = record.seq().len();

                    let indexes = chunk_fasta(record, length);
                    let adjacents = calculate_indexes(indexes);
                    let formatted = generate_explore_data(adjacents, id.clone(), length);
                    output_vec.append(
                        &mut merge_rotated_repeats(
                            formatted.unwrap_or(vec![]),
                            length,
                            id.clone(),
                            threshold,
                            &mut explore_file,
                            seq_len,
                            dist_from_chromosome_end,
                        )
                        .unwrap_or(vec![]),
                    );

                    println!("[+]\tChromosome {} processed", id);
                }
            }
        }
        println!("[+]\tFinished searching genome");
        // print likely telomeric repeat
        get_single_telomeric_repeat_estimate(&mut output_vec);
    }

    // split the fasta into chunks of size k, where k is the potential telomeric repeat length
    // compare consecutive iterations of these chunks and decide if they are equal
    // record the position of the repeat and the sequence.

    pub struct ChunkedFasta {
        pub position: usize,
        pub sequence: String,
    }

    fn chunk_fasta(sequence: bio::io::fasta::Record, chunk_length: usize) -> Vec<ChunkedFasta> {
        let chunks = sequence.seq().chunks(chunk_length);
        let chunks_plus_one =
            sequence.seq()[chunk_length..sequence.seq().len()].chunks(chunk_length);

        // store the index positions and adjacent equivalent sequences.
        let mut indexes = Vec::new();
        // need this otherwise we lose the position in the sequence.
        // need to check this is actually correct.
        let mut pos = 0;

        // this is the heavy lifting.
        for (_i, (a, b)) in chunks.zip(chunks_plus_one).enumerate() {
            // if chunk contains N, skip.
            if a.contains(&78) || b.contains(&78) {
                pos += chunk_length;
                continue;
            } else if a == b {
                indexes.push(ChunkedFasta {
                    position: pos,
                    sequence: str::from_utf8(a).unwrap().to_owned(),
                });
            }
            pos += chunk_length;
        }
        indexes
    }

    // takes the positions from chunk_fasta
    // take the current iteration position away from next iteration position

    pub struct RepeatRuns {
        pub position: usize,
        pub subtracted_position: usize,
        pub sequence: String,
    }

    fn calculate_indexes(indexes: Vec<ChunkedFasta>) -> Vec<RepeatRuns> {
        let mut adjacent_indexes = Vec::new();
        for i in 1..indexes.len() {
            if i + 1 < indexes.len() {
                adjacent_indexes.push(RepeatRuns {
                    position: indexes[i].position,
                    subtracted_position: indexes[i + 1].position - indexes[i].position,
                    sequence: indexes[i].sequence.clone(),
                })
            }
        }
        adjacent_indexes
    }

    #[derive(Debug, Clone)]
    pub struct TelomericRepeatExplore {
        pub start: usize,
        pub end: usize,
        pub count: i32,
        pub sequence: String,
        pub sequence_len: usize,
    }

    // inital data generation, using the indexes above
    // if there are a run of indexes == chunk_size,
    // these are the repeats we are looking for.

    fn generate_explore_data(
        adjacent_indexes: Vec<RepeatRuns>,
        id: String,
        chunk_length: usize,
    ) -> Option<Vec<TelomericRepeatExplore>> {
        if adjacent_indexes.is_empty() {
            println!(
                "[-]\tChromosome {}: No consecutive repeats of length {} were identified.",
                id, chunk_length
            );
            return None;
        }
        // collect all of this vector into a nice format to print for the moment.
        // iteration
        let mut it = 0;
        let mut count = 0;
        // threshold could become an input parameter
        // this is the first level of filtering which is useful
        // as most matches only occur twice (once repeated)
        // this local threshold parameter is quite interesting, too low and you get too much
        // output, but too high and you miss information.
        let local_threshold = 0;
        let mut start_pos = 0;
        let mut start;
        let mut end;
        let mut potential_telomeric_repeats = Vec::new();

        loop {
            start = adjacent_indexes[start_pos].position;
            if it == adjacent_indexes.len() - 1 {
                break;
            }
            if adjacent_indexes[it].subtracted_position == chunk_length {
                count += 1;
                it += 1;
            } else {
                end = adjacent_indexes[it].position;
                if count > local_threshold {
                    potential_telomeric_repeats.push(TelomericRepeatExplore {
                        start: start,
                        end: end,
                        count: count,
                        sequence: adjacent_indexes[it].sequence.clone(),
                        sequence_len: chunk_length,
                    });
                }

                start_pos = it;
                count = 1;
                it += 1;
            }
        }
        Some(potential_telomeric_repeats)
    }

    // basically pretty prints generate_explore_data
    // but crucially aggregates runs of records which are
    // string rotations of one another, yielding better summaries.

    // TODO: can the sequences be summarised? I.e. ID reverse complement sets.
    // -> Option<std::io::Result<()>> {

    #[derive(Debug)]
    pub struct FormatTelomericRepeat {
        sequence: String,
        count: i32,
        sequence_len: usize,
    }

    fn merge_rotated_repeats<T: std::io::Write>(
        data: Vec<TelomericRepeatExplore>,
        chunk_length: usize,
        id: String,
        threshold: i32,
        file: &mut LineWriter<T>,
        seq_len: usize,
        dist_from_chromosome_end: usize,
    ) -> Option<Vec<FormatTelomericRepeat>> {
        let mut output_vec: Vec<FormatTelomericRepeat> = Vec::new();

        if data.is_empty() {
            println!(
                "[-]\tChromosome {}: No consecutive repeats of length {} were identified.",
                id, chunk_length
            );
            return None;
        }

        let mut it = 0;
        let mut count = data[0].count;
        let mut start_index = 0;
        let mut start;
        let mut end;

        loop {
            // the starting value for the first result
            start = data[start_index].start;
            end = data[it].end;
            // explicit break in the loop
            if it == data.len() - 1 {
                // if all telomere repeat to the end, this is not printed.
                if count > threshold {
                    if start < dist_from_chromosome_end || end > seq_len - dist_from_chromosome_end
                    {
                        writeln!(
                            file,
                            "{},{},{},{},{},{}",
                            id, start, end, count, data[it].sequence, chunk_length
                        )
                        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
                        output_vec.push(FormatTelomericRepeat {
                            sequence: data[it].sequence.clone(),
                            count: count,
                            sequence_len: chunk_length,
                        });
                    }
                }
                // this seems weird, fix this?
                break Some(output_vec);
            }
            // if consecutive sequences are rotations
            if utils::string_rotation(&data[it].sequence, &data[it + 1].sequence) {
                // increment count by number of counts in next iteration
                count += data[it + 1].count;
                it += 1;
            } else {
                end = data[it].end;
                if count > threshold {
                    if start < dist_from_chromosome_end || end > seq_len - dist_from_chromosome_end
                    {
                        writeln!(
                            file,
                            "{},{},{},{},{},{}",
                            id, start, end, count, data[it].sequence, chunk_length
                        )
                        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
                        output_vec.push(FormatTelomericRepeat {
                            sequence: data[it].sequence.clone(),
                            count: count,
                            sequence_len: chunk_length,
                        });
                    }
                }
                it += 1;
                start_index = it;
                count = data[it].count;
            }
        }
    }

    // takes the final aggregation of potential telomeric repeats across chromosomes
    // and also potentially across different lengths, and tries to find the most likely
    // telomeric repeat. See utils::format_telomeric_repeat() for the explanation of the
    // formatting.
    // the algorithm here sorts the vector of formatted telomeric repeats by the sequence and their
    // length, then loops over this structure.

    fn get_single_telomeric_repeat_estimate(telomeric_repeats: &mut Vec<FormatTelomericRepeat>) {
        if telomeric_repeats.is_empty() {
            println!("[-]\tNo potential telomeric repeats found.");
            return;
        }
        // sort the sequence, then by sequence length
        // these sort in place
        // make from here downwards into a function
        telomeric_repeats.sort_by(|d1, d2| d1.sequence.cmp(&d2.sequence));
        telomeric_repeats.sort_by(|d1, d2| d2.sequence_len.cmp(&d1.sequence_len));

        // keep track of iterations (only incremented when a string rotation NOT matched...)
        let mut it = 0;
        // keep track of the lengths of the telomeric_repeats vec over its lifetime.
        let mut len_vec = Vec::new();

        // now loop
        loop {
            // create a vector of lengths of the telomeric repeats
            // as if all goes well, elements are removed.
            len_vec.push(telomeric_repeats.len());
            // we only increment the iteration if a string rotation is found
            // if the iteration reaches the length of the reduced telomeric repeats
            if it == telomeric_repeats.len() - 1 {
                // i.e. sort and go back to start...
                it = 0;
                telomeric_repeats.sort_by(|d1, d2| {
                    d2.sequence_len
                        .cmp(&d1.sequence_len)
                        .then(d2.count.cmp(&d1.count))
                        .then(d1.sequence.cmp(&d2.sequence))
                });
                // if we are left with only one element, this is the best case scenario
                // so break!
                if telomeric_repeats.len() == 1 {
                    break;
                }
                // if we get to the situation where in the length vector, the last two elements
                // have the same length, there is not going to be a more optimal solution.
                // so break!
                // should the 2 be 1 here..?
                if len_vec.len() > 2 {
                    if len_vec[len_vec.len() - 1] == len_vec[len_vec.len() - 2] {
                        break;
                    }
                }
            }

            if utils::string_rotation(
                &telomeric_repeats[it].sequence,
                &telomeric_repeats[it + 1].sequence,
            ) || utils::string_rotation(
                &utils::reverse_complement(&telomeric_repeats[it].sequence),
                &telomeric_repeats[it + 1].sequence,
            ) {
                telomeric_repeats[it].count =
                    telomeric_repeats[it].count + telomeric_repeats[it + 1].count;
                telomeric_repeats.remove(it + 1);
            } else {
                it += 1;
            }
        }
        // take the max count and then report the sequence
        let max = telomeric_repeats.iter().max_by_key(|i| i.count);
        // is this unwrap safe?
        println!(
            "[+]\tThe likely telomeric repeat is: {}",
            utils::format_telomeric_repeat(max.unwrap().sequence.clone())
        )
    }
}
