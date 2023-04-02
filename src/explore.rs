use crate::{utils, SubCommand};
use anyhow::Result;
use bio::io::fasta;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::path::PathBuf;
use std::str;
use std::sync::mpsc::channel;

static REPEAT_PERIOD_THRESHOLD: usize = 3;

/// The function called from `tidk explore`. It takes the [`clap::Argmatches`]
/// from the user and also a [`SubCommand`].
pub fn explore(matches: &clap::ArgMatches, sc: SubCommand) -> Result<()> {
    // parse arguments from main
    let input_fasta = matches
        .get_one::<PathBuf>("fasta")
        .expect("errored by clap");
    let length = *matches.get_one::<usize>("length").expect("errored by clap");

    // if length is not set, these are the lengths (and length itself is set to zero)
    let minimum = *matches
        .get_one::<usize>("minimum")
        .expect("errored by clap");
    let maximum = *matches
        .get_one::<usize>("maximum")
        .expect("errored by clap");

    let threshold = *matches
        .get_one::<i32>("threshold")
        .expect("errored by clap");

    let dist_from_chromosome_end = *matches.get_one::<f64>("distance").expect("errored by clap");

    let verbose = matches.get_flag("verbose");

    // to report the telomeres...
    let mut output_vec: Vec<RepeatPositions> = Vec::new();
    // i.e. if you chose a length, as opposed to a minmum/maximum
    if length > 0 {
        eprintln!(
            "[+]\tExploring genome for potential telomeric repeats of length: {}",
            length
        );
        let reader = fasta::Reader::from_file(input_fasta)?;

        // try parallelising
        let (sender, receiver) = channel();

        reader
            .records()
            .par_bridge()
            .for_each_with(sender, |s, record| {
                let record = record.expect("[-]\tError during fasta record parsing.");
                let id = record.id().to_owned();
                let seq_len = record.seq().len();

                let indexes =
                    chunk_fasta(record, length, dist_from_chromosome_end, seq_len, verbose);

                if let Some(r) = calculate_indexes(indexes, length, verbose, id, threshold as usize)
                {
                    s.send(r).expect("Did not send!");
                }
            });

        // this bit is a little chaotic
        // collect output into a vector
        let mut output = receiver.into_iter().collect();

        output_vec.append(&mut output);
    } else {
        // if a range was chosen.
        eprintln!(
            "[+]\tExploring genome for potential telomeric repeats between lengths {} and {}.",
            minimum, maximum
        );
        for length in minimum..maximum + 1 {
            eprintln!("[+]\t\tFinding telomeric repeat length: {}", length);

            // have to call reader in the loop, as otherwise `reader` doesn't live long enough.
            // I expect it's not an expensive call anyway.
            let reader = fasta::Reader::from_file(input_fasta)?;

            // try parallelising
            let (sender, receiver) = channel();
            reader
                .records()
                .par_bridge()
                .for_each_with(sender, |s, record| {
                    let record = record.expect("[-]\tError during fasta record parsing.");
                    let id = record.id().to_owned();
                    let seq_len = record.seq().len();

                    let indexes =
                        chunk_fasta(record, length, dist_from_chromosome_end, seq_len, verbose);

                    if let Some(r) =
                        calculate_indexes(indexes, length, verbose, id, threshold as usize)
                    {
                        s.send(r).expect("Did not send!");
                    }
                });
            let mut output = receiver.iter().collect();

            output_vec.append(&mut output);
        }
    }
    eprintln!("[+]\tFinished searching genome");
    eprintln!("[+]\tGenerating output");

    let mut repeat_postitions = RepeatPositions::new();
    for mut el in output_vec {
        repeat_postitions.add(&mut el.0);
    }

    // print likely telomeric repeat
    // costly calculation if threshold is too low.
    let est = get_telomeric_repeat_estimates(&mut repeat_postitions)?;

    println!("canonical_repeat_unit\tcount");
    for (cru, count) in est {
        println!("{}\t{}", cru, count);
    }

    // optional log file
    sc.log(matches)?;

    Ok(())
}

/// A chunked fasta segment with a position and a sequence.
/// We split the fasta into chunks of size k, where k is the
/// potential telomeric repeat length. Consecutive iterations
/// of these chunks are compared for equality.
#[derive(Debug, PartialEq, Eq)]
pub struct ChunkedFasta {
    /// Position of the sequence in the
    /// fasta file.
    pub position: usize,
    /// The sequence itself.
    pub sequence: String,
}

/// Chunk a fasta into a [`Vec<ChunkedFasta>`], i.e. split a fasta into chunks
/// and compare adjacent chunks for equality. Store the positions and sequences
/// if they are equivalent.
fn chunk_fasta(
    sequence: bio::io::fasta::Record,
    chunk_length: usize,
    dist_from_chromosome_end: f64,
    seq_len: usize,
    verbose: bool,
) -> Vec<ChunkedFasta> {
    // seems a bit slower now for pacbio data?
    // probably need more error handling in this function
    let dist = (seq_len as f64 * dist_from_chromosome_end).floor() as usize;

    let filtered_sequence1 = &sequence.seq()[0..dist];
    let filtered_sequence2 = &sequence.seq()[(seq_len - dist)..];

    let filtered_sequence = [filtered_sequence1, filtered_sequence2].concat();

    let filtered_sequence_len = filtered_sequence.len();

    let chunks = filtered_sequence.chunks(chunk_length);
    // catch edge cases where chunk length greater than sequence length.
    if filtered_sequence_len <= chunk_length {
        if verbose {
            eprintln!(
                "[-]\tChunk length ({}) greater than filtered sequence length ({}) for {}
[-]\tConsider increasing proportion of chromosome length covered. Skipping.",
                chunk_length,
                filtered_sequence_len,
                sequence.id()
            );
        }
        return vec![];
    }

    let chunks_plus_one =
        filtered_sequence[chunk_length..filtered_sequence_len].chunks(chunk_length);

    // store the index positions and adjacent equivalent sequences.
    let mut indexes = Vec::new();
    // need this otherwise we lose the position in the sequence.
    // need to check this is actually correct.
    let mut pos = 0;

    // this is the heavy lifting.
    // can use the enumerate to check whether the position is < dist from start or > dist from end.
    for (a, b) in chunks.zip(chunks_plus_one) {
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

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct RepeatPosition {
    id: String,
    pub start: usize,
    pub end: usize,
    pub sequence: String,
}

impl RepeatPosition {
    fn get_count(&self) -> usize {
        (self.end - self.start) / self.sequence.len()
    }
    // true if it's not a simple repeat
    fn is_simple_repeat(&self) -> bool {
        check_telomeric_repeat(&self.sequence)
    }
}

#[derive(Debug)]
pub struct RepeatPositions(Vec<RepeatPosition>);

impl RepeatPositions {
    fn new() -> Self {
        Self(Vec::new())
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn add(&mut self, elem: &mut Vec<RepeatPosition>) {
        self.0.append(elem);
    }

    fn filter_by_frequency(&mut self, frequency: usize) -> Self {
        let inner: &Vec<RepeatPosition> = &self
            .0
            .clone() // can we remove this?
            .into_iter()
            .filter(|e| e.get_count() > frequency && !e.is_simple_repeat())
            .collect();

        Self(inner.to_vec())
    }
}

// logic messed up here - the start/end don't exclusively include telomeric repeats.
// it's merging two consective runs, even if they are separated by non (canonical)-telomeric sequence.
fn calculate_indexes(
    indexes: Vec<ChunkedFasta>,
    chunk_length: usize,
    verbose: bool,
    id: String,
    frequency: usize,
) -> Option<RepeatPositions> {
    let mut start = 0usize;

    let mut collection: Vec<RepeatPosition> = Vec::new();

    let mut iter = indexes.iter().zip(indexes.iter().skip(1)).peekable();

    while let Some((
        ChunkedFasta {
            position: position1,
            sequence: sequence1,
        },
        ChunkedFasta {
            position: position2,
            sequence: sequence2,
        },
    )) = iter.next()
    {
        if iter.peek().is_none() {
            // this is techinically incorrect - as this will almost always
            // overshoot the last index of the genome.
            collection.push(RepeatPosition {
                id: id.clone(),
                start,
                end: *position2 + chunk_length,
                sequence: sequence1.to_string(),
            });
        } else if sequence1 == sequence2 {
            eprintln!("The telomeric repeat for {} is: {}", id.clone(), sequence1);
            continue;
        } else if sequence1 != sequence2 {
            collection.push(RepeatPosition {
                id: id.clone(),
                start,
                end: *position1 + chunk_length,
                sequence: sequence1.to_string(),
            });
            start = *position2;
        }
    }
    if collection.is_empty() {
        if verbose {
            eprintln!(
                "[-]\t\tChromosome {}: No consecutive repeats of length {} were identified.",
                id, chunk_length
            );
        }
        None
    } else {
        let filtered_repeat_positions = RepeatPositions(collection).filter_by_frequency(frequency);
        Some(filtered_repeat_positions)
    }
}

/// a quick method to check if
/// a sequence looks like it is not
/// a telomeric repeat
fn check_telomeric_repeat(sequence: &str) -> bool {
    let repeat_period = check_repeats(sequence);
    repeat_period < REPEAT_PERIOD_THRESHOLD
}

/// Takes the final aggregation of potential telomeric repeats across
/// chromosomes and also potentially across different lengths and tries
/// to find the most likely telomeric repeat. See [`utils::format_telomeric_repeat()`]
/// for the explanation of the formatting.
fn get_telomeric_repeat_estimates(
    telomeric_repeats: &mut RepeatPositions,
) -> Result<Vec<(String, i32)>> {
    // can't get all the combinations if we only have 1 element.
    eprintln!("{:?}", telomeric_repeats);
    if telomeric_repeats.len() == 1 {
        let count = telomeric_repeats.0[0].get_count();
        let seq = &telomeric_repeats.0[0].sequence;
        return Ok(vec![(seq.clone(), count as i32)]);
    }
    // we need to compare all elements against all others
    let mut map: HashMap<String, i32> = HashMap::new();
    // so we don't compare the same thing twice.
    let mut tracker: Vec<usize> = Vec::new();
    // create all combinations of indices
    let it = (0..telomeric_repeats.len()).combinations(2);

    // iterate over combinations
    for comb in it {
        let first = &telomeric_repeats.0[comb[0]];
        let second = &telomeric_repeats.0[comb[1]];
        // if the combination is a string rotation (or its reverse complement)
        // then combine
        if utils::string_rotation(&first.sequence, &second.sequence)
            || utils::string_rotation(
                &utils::reverse_complement(&first.sequence),
                &second.sequence,
            )
            || utils::string_rotation(
                &first.sequence,
                &utils::reverse_complement(&second.sequence),
            )
        {
            // if comb[0] || comb[1] not in tracker...
            // as we already added the contents of the tracked telomeric repeats
            // we do not want to count them again.
            if !tracker.contains(&comb[0]) && !tracker.contains(&comb[1]) {
                let count = map
                    // relies on the telomeric repeat string resolving to a 'canonical'
                    // or unique form of the string, see utils::lms()
                    .entry(utils::lms(&first.sequence, &second.sequence))
                    .or_insert(first.get_count() as i32 + second.get_count() as i32);
                *count += first.get_count() as i32 + second.get_count() as i32;

                tracker.push(comb[0]);
                tracker.push(comb[1]);
            }
        }
    }

    let mut count_vec: Vec<_> = map.into_iter().collect();
    count_vec.sort_by(|a, b| b.1.cmp(&a.1));
    filter_count_vec(&mut count_vec)?;

    Ok(count_vec)
}

/// Returns the shortest period of repetition in s.
/// If s does not repeat, returns the number of characters in s.
///
/// See https://users.rust-lang.org/t/checking-simple-repeats-in-strings/79729
/// for a small discussion.
fn check_repeats(s: &str) -> usize {
    let mut delays: BTreeMap<_, std::str::Chars> = BTreeMap::new();
    for (i, c) in s.chars().enumerate() {
        delays.retain(|_, iter| iter.next() == Some(c));
        delays.insert(i + 1, s.chars());
    }
    delays.into_keys().next().unwrap()
}

/// A function to filter the final count vec of certain kinds of
/// short repeat which are probably not telomeric repeats. Should
/// clean up output dramatically.
///
/// These are:
/// - Monomeric
/// - Dimeric
/// - Trimeric
fn filter_count_vec(v: &mut Vec<(String, i32)>) -> Result<()> {
    // monomers
    // not sure I need this.
    v.retain(|(s, _)| {
        let repeat_period = check_repeats(s);
        repeat_period > REPEAT_PERIOD_THRESHOLD
    });

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const R1: &str = "AAAAAAAA";
    const R2: &str = "ATATATAT";
    const R3: &str = "AATAATAAT";
    // A tiny genome with telomeric repeats at the end
    const GENOME: &str =
        "AACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTTATAGGGACATTGACCGAAGGGGGCACATAGACAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCTAACCT";
    const CHUNK_LENGTH: usize = 5;
    const DIST_FROM_CHROM_END: f64 = 0.2;

    fn generate_chunks() -> Vec<ChunkedFasta> {
        let record = bio::io::fasta::Record::with_attrs("id1", None, GENOME.as_bytes());
        chunk_fasta(
            record,
            CHUNK_LENGTH,
            DIST_FROM_CHROM_END,
            GENOME.len(),
            false,
        )
    }

    fn generate_indexes() -> RepeatPositions {
        let chunks = generate_chunks();
        calculate_indexes(chunks, CHUNK_LENGTH, false, "test".into(), 0).unwrap()
    }

    #[test]
    fn check_repeat_period1() {
        let p = check_repeats(R1);
        assert_eq!(p, 1)
    }
    #[test]
    fn check_repeat_period2() {
        let p = check_repeats(R2);
        assert_eq!(p, 2)
    }
    #[test]
    fn check_repeat_period3() {
        let p = check_repeats(R3);
        assert_eq!(p, 3)
    }
    #[test]
    fn test_chunk_fasta() {
        let chunks = generate_chunks();
        assert_eq!(
            vec![
                ChunkedFasta {
                    position: 0,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 5,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 10,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 15,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 20,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 25,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 30,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 35,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 40,
                    sequence: "AACCT".into()
                }
            ],
            chunks
        )
    }
    #[test]
    fn test_index_calculation1() {
        let indices = generate_indexes();
        assert_eq!(
            indices.0,
            // i.e. there are 9 matches for AACCT
            vec![RepeatPosition {
                start: 0,
                end: 45,
                sequence: "AACCT".into()
            }]
        )
    }
    #[test]
    fn test_get_telomeric_repeat_estimates() {
        let mut indices = generate_indexes();
        let estimates = get_telomeric_repeat_estimates(&mut indices).unwrap();

        assert_eq!(estimates, vec![("AACCT".into(), 9)])
    }
}
