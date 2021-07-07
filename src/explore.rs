pub mod explore {
    use crate::utils::utils;
    use bio::io::fasta;
    use clap::value_t;
    use itertools::Itertools;
    use rayon::prelude::*;
    use std::collections::HashMap;
    use std::fs::{create_dir_all, File};
    use std::io::prelude::*;
    use std::io::LineWriter;
    use std::str;
    use std::sync::mpsc::channel;

    // function called from main.rs
    // TODO: sort out output file formats.

    pub fn explore(matches: &clap::ArgMatches) {
        // parse arguments from main
        let input_fasta = matches.value_of("fasta").unwrap();
        let length = value_t!(matches.value_of("length"), usize).unwrap_or_else(|e| e.exit());

        // if length is not set, these are the lengths (and length itself is set to zero)
        let minimum = value_t!(matches.value_of("minimum"), usize).unwrap_or_else(|e| e.exit());
        let maximum = value_t!(matches.value_of("maximum"), usize).unwrap_or_else(|e| e.exit());

        let threshold = value_t!(matches.value_of("threshold"), i32).unwrap_or_else(|e| e.exit());
        let output = matches.value_of("output").unwrap();
        let extension =
            value_t!(matches.value_of("extension"), String).unwrap_or_else(|e| e.exit());

        let dist_from_chromosome_end =
            value_t!(matches.value_of("distance"), usize).unwrap_or_else(|e| e.exit());

        // create directory for output
        if let Err(e) = create_dir_all("./explore/") {
            println!("[-]\tCreate directory error: {}", e.to_string());
        }

        // create file
        let file_name = format!(
            "./explore/{}{}{}",
            output, "_telomeric_locations.", extension
        );
        let explore_file = File::create(&file_name).unwrap();
        let mut explore_file = LineWriter::new(explore_file);

        let putative_telomeric_file = format!("./explore/{}{}", output, ".txt");
        let putative_telomeric_file_txt = File::create(&putative_telomeric_file).unwrap();
        let mut putative_telomeric_file_txt = LineWriter::new(putative_telomeric_file_txt);

        // add header to txt file
        writeln!(
            putative_telomeric_file_txt,
            "telomeric_repeat\treverse_complement\tfrequency"
        )
        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));

        // to report the telomeres...
        let mut output_vec: Vec<FormatTelomericRepeat> = Vec::new();
        let mut output_vec_bed: Vec<TsvTelomericRepeat> = Vec::new();
        // i.e. if you chose a length, as opposed to a minmum/maximum
        if length > 0 {
            println!(
                "[+]\tExploring genome for potential telomeric repeats of length: {}",
                length
            );
            let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

            // try parallelising
            let (sender, receiver) = channel();

            reader
                .records()
                .par_bridge()
                .for_each_with(sender, |s, record| {
                    let record = record.expect("[-]\tError during fasta record parsing.");
                    let id = record.id().to_owned();
                    let seq_len = record.seq().len();

                    let indexes = chunk_fasta(record, length);
                    let adjacents = calculate_indexes(indexes);
                    let formatted =
                        generate_explore_data(adjacents, id.clone(), length).unwrap_or(vec![]);

                    s.send(
                        merge_rotated_repeats(
                            formatted,
                            length,
                            &id,
                            threshold,
                            seq_len,
                            dist_from_chromosome_end,
                        )
                        .unwrap_or(Output {
                            telomeric_repeats: vec![],
                            bed_file: vec![],
                        }),
                    )
                    .expect("[-]\tDid not send.");
                });

            // this bit is a little chaotic
            // collect output into a vector
            let output: Vec<Output> = receiver.iter().collect();

            // so it can be cloned here
            // to extract repeats
            let mut telomeric_repeats = output
                .clone()
                .iter()
                .map(|a| a.telomeric_repeats.clone())
                .flatten()
                .collect();
            // and appended to the output vec
            output_vec.append(&mut telomeric_repeats);

            // and also cloned here
            let mut bed_file: Vec<TsvTelomericRepeat> = output
                .clone()
                .iter()
                .map(|a| a.bed_file.clone())
                .flatten()
                .collect();
            // to get a bed file of all potential repeat locations.
            output_vec_bed.append(&mut bed_file);
        } else {
            // if a range was chosen.
            println!(
                "[+]\tExploring genome for potential telomeric repeats between lengths {} and {}.",
                minimum, maximum
            );
            for length in minimum..maximum + 1 {
                println!("[+]\t\tFinding telomeric repeat length: {}", length);

                // have to call reader in the loop, as otherwise `reader` doesn't live long enough.
                // I expect it's not an expensive call anyway.
                let reader = fasta::Reader::from_file(input_fasta).expect("[-]\tPath invalid.");

                // try parallelising
                let (sender, receiver) = channel();
                reader
                    .records()
                    .par_bridge()
                    .for_each_with(sender, |s, record| {
                        let record = record.expect("[-]\tError during fasta record parsing.");
                        let id = record.id().to_owned();
                        let seq_len = record.seq().len();

                        let indexes = chunk_fasta(record, length);
                        let adjacents = calculate_indexes(indexes);
                        let formatted =
                            generate_explore_data(adjacents, id.clone(), length).unwrap_or(vec![]);

                        s.send(
                            merge_rotated_repeats(
                                formatted,
                                length,
                                &id,
                                threshold,
                                seq_len,
                                dist_from_chromosome_end,
                            )
                            .unwrap_or(Output {
                                telomeric_repeats: vec![],
                                bed_file: vec![],
                            }),
                        )
                        .expect("[-]\tDid not send.");
                    });
                // this bit is a little chaotic
                // collect output into a vector
                let output: Vec<Output> = receiver.iter().collect();

                // so it can be cloned here
                // to extract repeats
                let mut telomeric_repeats = output
                    .clone()
                    .iter()
                    .map(|a| a.telomeric_repeats.clone())
                    .flatten()
                    .collect();
                // and appended to the output vec
                output_vec.append(&mut telomeric_repeats);

                // and also cloned here
                let mut bed_file: Vec<TsvTelomericRepeat> = output
                    .clone()
                    .iter()
                    .map(|a| a.bed_file.clone())
                    .flatten()
                    .collect();
                // to get a bed file of all potential repeat locations.
                output_vec_bed.append(&mut bed_file);
            }
        }
        println!("[+]\tFinished searching genome");
        println!("[+]\tGenerating output");
        // print likely telomeric repeat
        // costly calculation if threshold is too low.
        get_telomeric_repeat_estimates(&mut output_vec, &mut putative_telomeric_file_txt);
        // write the bed file
        writeln!(
            explore_file,
            "id\tstart_pos\tend_pos\trepeat_number\trepeat_sequence\tsequence_length"
        )
        .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
        for line in output_vec_bed {
            writeln!(
                explore_file,
                "{}\t{}\t{}\t{}\t{}\t{}",
                line.id, line.start, line.end, line.count, line.sequence, line.seq_len
            )
            .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
        }
    }

    // split the fasta into chunks of size k, where k is the potential telomeric repeat length
    // compare consecutive iterations of these chunks and decide if they are equal
    // record the position of the repeat and the sequence.

    pub struct ChunkedFasta {
        pub position: usize,
        pub sequence: String,
    }

    // is it possible to filter this iterator based on distance from chromosome end..?

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
                    sequence: str::from_utf8(a).unwrap().to_uppercase(),
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
            eprintln!(
                "[-]\t\tChromosome {}: No consecutive repeats of length {} were identified.",
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

    #[derive(Debug, Clone, PartialEq, Eq, Hash)]
    pub struct FormatTelomericRepeat {
        sequence: String,
        count: i32,
        sequence_len: usize,
    }

    #[derive(Debug, Clone)]
    pub struct TsvTelomericRepeat {
        id: String,
        start: usize,
        end: usize,
        count: i32,
        sequence: String,
        seq_len: usize,
    }

    #[derive(Debug, Clone)]
    pub struct Output {
        telomeric_repeats: Vec<FormatTelomericRepeat>,
        bed_file: Vec<TsvTelomericRepeat>,
    }

    fn merge_rotated_repeats<'a>(
        data: Vec<TelomericRepeatExplore>,
        chunk_length: usize,
        id: &'a str,
        threshold: i32,
        seq_len: usize,
        dist_from_chromosome_end: usize,
    ) -> Option<Output> {
        let mut output_vec: Vec<FormatTelomericRepeat> = Vec::new();
        let mut output_vec_tsv: Vec<TsvTelomericRepeat> = Vec::new();

        // check for absence of data in the vector
        if data.is_empty() {
            eprintln!(
                "[-]\t\tChromosome {}: No consecutive repeats of length {} were identified.",
                id.clone(),
                chunk_length
            );
            return None;
        }

        // keep track of iterations
        let mut it = 0;
        // initiate count as the count of the first element
        let mut count = data[0].count;
        // increment this only when no more string rotations are found
        let mut start_index = 0;
        // start and end chromosome positions of aggregated repeats.
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
                        // to output later
                        output_vec_tsv.push(TsvTelomericRepeat {
                            id: id.to_owned(),
                            start: start,
                            end: end,
                            count: count,
                            sequence: data[it].sequence.clone(),
                            seq_len: chunk_length,
                        });
                        // and collect for guessing telomeric repeat
                        output_vec.push(FormatTelomericRepeat {
                            sequence: data[it].sequence.clone(),
                            count: count,
                            sequence_len: chunk_length,
                        });
                    }
                }
                break Some(Output {
                    telomeric_repeats: output_vec,
                    bed_file: output_vec_tsv,
                });
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
                        // to output later
                        output_vec_tsv.push(TsvTelomericRepeat {
                            id: id.to_owned(),
                            start: start,
                            end: end,
                            count: count,
                            sequence: data[it].sequence.clone(),
                            seq_len: chunk_length,
                        });
                        // and collect for guessing telomeric repeat
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

    fn get_telomeric_repeat_estimates<T: std::io::Write>(
        telomeric_repeats: &mut Vec<FormatTelomericRepeat>,
        putative_telomeric_file: &mut LineWriter<T>,
    ) {
        if telomeric_repeats.is_empty() {
            eprintln!("[-]\tNo potential telomeric repeats found.");
            return;
        }

        // we need to compare all elements against all others
        let mut map: HashMap<String, i32> = HashMap::new();
        // so we don't compare the same thing twice.
        let mut tracker: Vec<usize> = Vec::new();
        // create all combinations of indices
        let it = (0..telomeric_repeats.len()).combinations(2);

        // iterate over combinations
        for comb in it {
            // if the combination is a string rotation (or its reverse complement)
            // then combine
            if utils::string_rotation(
                &telomeric_repeats[comb[0]].sequence,
                &telomeric_repeats[comb[1]].sequence,
            ) || utils::string_rotation(
                &utils::reverse_complement(&telomeric_repeats[comb[0]].sequence),
                &telomeric_repeats[comb[1]].sequence,
            ) || utils::string_rotation(
                &utils::reverse_complement(&telomeric_repeats[comb[1]].sequence),
                &telomeric_repeats[comb[0]].sequence,
            ) {
                // if comb[0] || comb[1] not in tracker...
                // as we already added the contents of the tracked telomeric repeats
                // we do not want to count them again.
                if !tracker.contains(&comb[0]) && !tracker.contains(&comb[1]) {
                    let count = map
                        // relies on the telomeric repeat string resolving to a 'canonical'
                        // or unique form of the string, see utils::lms()
                        .entry(utils::lms(
                            &telomeric_repeats[comb[0]].sequence,
                            &telomeric_repeats[comb[1]].sequence,
                        ))
                        .or_insert(
                            telomeric_repeats[comb[0]].count + telomeric_repeats[comb[1]].count,
                        );
                    *count += telomeric_repeats[comb[0]].count + telomeric_repeats[comb[1]].count;

                    tracker.push(comb[0]);
                    tracker.push(comb[1]);
                }
            }
        }

        let mut count_vec: Vec<_> = map.iter().collect();
        count_vec.sort_by(|a, b| b.1.cmp(a.1));

        let mut it = 0;
        for (seq, count) in count_vec {
            if it == 0 {
                println!(
                    "[+]\tThe likely telomeric repeat is: {}, found {} times.",
                    seq, count
                );
            }
            writeln!(
                putative_telomeric_file,
                "{}\t{}\t{}",
                seq,
                utils::reverse_complement(seq),
                count
            )
            .unwrap_or_else(|_| println!("[-]\tError in writing to file."));
            it += 1;
        }
    }
}
