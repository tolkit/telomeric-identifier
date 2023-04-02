use crate::{utils, SubCommand};
use anyhow::bail;
use anyhow::Result;
use bio::io::fasta;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::path::PathBuf;
use std::str;
use std::sync::mpsc::channel;

// when distance == 1, we get lower estimate of telomeric repeat number
// than if we use distance == 0.1
// it's somehow not splitting RepeatPositions correctly.

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

    if dist_from_chromosome_end > 0.5 {
        bail!("Distance from chromosome end as a proportion can't be more than 0.5.")
    }

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

                let sequences = split_seq_by_distance(record, dist_from_chromosome_end, seq_len);

                for sequence in sequences {
                    let indexes = chunk_fasta(sequence, length, verbose, id.clone());

                    if let Some(r) =
                        calculate_indexes(indexes, length, verbose, id.clone(), threshold as usize)
                    {
                        s.send(r).expect("Did not send!");
                    }
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

                    let sequences =
                        split_seq_by_distance(record, dist_from_chromosome_end, seq_len);

                    for sequence in sequences {
                        let indexes = chunk_fasta(sequence, length, verbose, id.clone());

                        if let Some(r) = calculate_indexes(
                            indexes,
                            length,
                            verbose,
                            id.clone(),
                            threshold as usize,
                        ) {
                            s.send(r).expect("Did not send!");
                        }
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

pub fn split_seq_by_distance(
    sequence: bio::io::fasta::Record,
    dist_from_chromosome_end: f64,
    seq_len: usize,
) -> [Vec<u8>; 2] {
    let dist = (seq_len as f64 * dist_from_chromosome_end).ceil() as usize;
    let filtered_sequence1 = sequence.seq()[0..dist].to_vec();
    let filtered_sequence2 = sequence.seq()[(seq_len - dist)..].to_vec();
    [filtered_sequence1, filtered_sequence2]
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
    sequence: Vec<u8>,
    chunk_length: usize,
    verbose: bool,
    id: String,
) -> Vec<ChunkedFasta> {
    let sequence_len = sequence.len();
    let chunks = sequence.chunks(chunk_length);
    // catch edge cases where chunk length greater than sequence length.
    if sequence_len <= chunk_length {
        if verbose {
            eprintln!(
                "[-]\tChunk length ({}) greater than filtered sequence length ({}) for {}
[-]\tConsider increasing proportion of chromosome length covered. Skipping.",
                chunk_length, sequence_len, id
            );
        }
        return vec![];
    }

    let chunks_plus_one = sequence.chunks(chunk_length).skip(1);

    // store the index positions and adjacent equivalent sequences.
    let mut indexes = Vec::new();
    // need this otherwise we lose the position in the sequence.
    let mut pos = 0;
    let mut is_first_consecutive = true;

    // this is the heavy lifting.
    // can use the enumerate to check whether the position is < dist from start or > dist from end.
    // FIXME: how do I check if this is the first in a run?
    for (a, b) in chunks.zip(chunks_plus_one) {
        // if chunk contains N, skip.
        if a.contains(&78) || b.contains(&78) {
            pos += chunk_length;
            continue;
        } else if a == b {
            if is_first_consecutive {
                indexes.push(ChunkedFasta {
                    position: pos,
                    sequence: str::from_utf8(a).unwrap().to_uppercase(),
                });
                pos += chunk_length;
                indexes.push(ChunkedFasta {
                    position: pos,
                    sequence: str::from_utf8(a).unwrap().to_uppercase(),
                });
                is_first_consecutive = false;
            } else {
                pos += chunk_length;
                indexes.push(ChunkedFasta {
                    position: pos,
                    sequence: str::from_utf8(a).unwrap().to_uppercase(),
                });
            }
        } else if a != b {
            pos += chunk_length;
            is_first_consecutive = true;
        }
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
    // group into HashMap<usize, Vec<RepeatPosition>>
    // where usize is the length of the telomeric repeat
    fn make_length_groups(&self) -> HashMap<usize, Vec<RepeatPosition>> {
        let mut groups = HashMap::new();
        for el in &self.0 {
            let length = el.sequence.len();
            groups.entry(length).or_insert(Vec::new()).push(el.clone());
        }
        groups
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
    // eprintln!("INDEXES: {:#?}", indexes);
    let mut start = 0usize;
    // let mut end = 0usize;

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
            // end = *position2 + chunk_length;
            // eprintln!("FINAL:: START: {}\tEND: {}", start, end);
            // FIXME: ITS HERE? ALSO NEED TO CHECK IF START/END IS DIFFERENT?
        } else if sequence1 == sequence2 && (position2 - position1) == sequence1.len() {
            // start = *position1;
            // end = *position2;
            // eprintln!("EQUAL:: START: {}\tEND: {}", start, end);
            continue;
        } else if !(sequence1 == sequence2 && (position2 - position1) == sequence1.len()) {
            // eprintln!("UNEQUAL:: START: {}\tEND: {}", start, end);
            collection.push(RepeatPosition {
                id: id.clone(),
                start,
                end: *position1 + chunk_length,
                sequence: sequence1.to_string(),
            });
            // end = *position1 + chunk_length;
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

/// check if a sequence looks like it is not
/// a telomeric repeat
fn check_telomeric_repeat(sequence: &str) -> bool {
    let repeat_period = check_repeats(sequence);
    repeat_period < REPEAT_PERIOD_THRESHOLD
}

/// take two repeat sequences of the same length and
/// return bool, if they represent the same canonical repeat
/// and the canonical repeat
fn test_repeats(repeat1: &RepeatPosition, repeat2: &RepeatPosition) -> (bool, Option<String>) {
    let r1_seq = &repeat1.sequence;
    let r2_seq = &repeat2.sequence;
    let is_equal = utils::string_rotation(r1_seq, r2_seq)
        || utils::string_rotation(&utils::reverse_complement(r1_seq), r2_seq)
        || utils::string_rotation(r1_seq, &utils::reverse_complement(r2_seq));

    if is_equal {
        (is_equal, Some(utils::lms(r1_seq, r2_seq)))
    } else {
        (is_equal, None)
    }
}

/// Takes the final aggregation of potential telomeric repeats across
/// chromosomes and also potentially across different lengths and tries
/// to find the most likely telomeric repeat. See [`utils::format_telomeric_repeat()`]
/// for the explanation of the formatting.
/// FIXME: we compare telomeric repeats of different lengths here! sort that out.
fn get_telomeric_repeat_estimates(
    telomeric_repeats: &mut RepeatPositions,
) -> Result<Vec<(String, i32)>> {
    let groups = telomeric_repeats.make_length_groups();

    // we need to compare all elements against all others
    let mut map: HashMap<String, i32> = HashMap::new();

    for (_, telomeric_repeats_i) in groups {
        // deal with this separately, as if there's only one repeat
        // estimated, we can't go further
        if telomeric_repeats_i.len() == 1 {
            let count = telomeric_repeats_i[0].get_count();
            let seq = &telomeric_repeats_i[0].sequence;
            map.insert(seq.clone(), count as i32);
            continue;
        }
        // so we don't compare the same thing twice.
        let mut tracker: Vec<usize> = Vec::new();
        // create all combinations of indices
        let it = (0..telomeric_repeats_i.len()).combinations(2);

        // iterate over combinations
        for comb in it {
            let first = &telomeric_repeats_i[comb[0]];
            let second = &telomeric_repeats_i[comb[1]];

            let (is_rotation, potential_sequence) = test_repeats(first, second);
            // eprintln!("first comb: {:?}\tsecond comb: {:?}", first, second);
            // if the combination is a string rotation (or its reverse complement)
            // then combine
            if is_rotation {
                let sequence = potential_sequence.unwrap();
                // if comb[0] || comb[1] not in tracker...
                // as we already added the contents of the tracked telomeric repeats
                // we do not want to count them again.
                let count = map
                    // relies on the telomeric repeat string resolving to a 'canonical'
                    // or unique form of the string, see utils::lms()
                    .entry(sequence)
                    .or_insert(first.get_count() as i32 + second.get_count() as i32);
                if !tracker.contains(&comb[0]) && !tracker.contains(&comb[1]) {
                    *count += first.get_count() as i32 + second.get_count() as i32;
                }
            }
            tracker.push(comb[0]);
            tracker.push(comb[1]);
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

    // three kinds of short repeat sequences
    // that aren't telomeric repeats.
    const R1: &str = "AAAAAAAA";
    const R2: &str = "ATATATAT";
    const R3: &str = "AATAATAAT";

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

    // A tiny genome with telomeric repeats at the end
    const GENOME: &str = "AACCTAACCTAACATATCGTAACCTAACCTAACCTAACCTAACATATCGTAACCTAACCT";
    //                                                  ^ repeats here
    const GENOME_2: &str = "AACCTAACCTTAAATTAAATAACCTAACCTAACCTAACCTTAAATTAAATAACCTAACCT";
    //                                                    ^ repeats here
    // we are looking at 5-mers
    const CHUNK_LENGTH: usize = 5;
    // include the whole sequence.
    const DIST_FROM_CHROM_END: f64 = 0.5;

    fn split_by_dist(genome: &str) -> [Vec<u8>; 2] {
        let record = bio::io::fasta::Record::with_attrs("id1", None, genome.as_bytes());
        split_seq_by_distance(record, DIST_FROM_CHROM_END, genome.len())
    }

    // GENOME/GENOME_2 are just two meta-repeats, so this should just be in half
    #[test]
    fn test_split_left() {
        let seq = &split_by_dist(GENOME)[0];
        let left = std::str::from_utf8(seq).unwrap();
        assert_eq!(left, "AACCTAACCTAACATATCGTAACCTAACCT")
    }

    #[test]
    fn test_split_right() {
        let seq = &split_by_dist(GENOME)[1];
        let left = std::str::from_utf8(seq).unwrap();
        assert_eq!(left, "AACCTAACCTAACATATCGTAACCTAACCT")
    }

    fn generate_chunks_left(genome: &str) -> Vec<ChunkedFasta> {
        let left = &split_by_dist(genome)[0];
        chunk_fasta(left.clone(), CHUNK_LENGTH, false, "".into())
    }

    #[test]
    fn test_chunks_left() {
        let chunks = generate_chunks_left(GENOME);
        // in the left
        assert_eq!(
            chunks,
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
                    position: 20,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 25,
                    sequence: "AACCT".into()
                }
            ]
        )
    }

    fn generate_chunks_right() -> Vec<ChunkedFasta> {
        let left = &split_by_dist(GENOME)[1];
        chunk_fasta(left.clone(), CHUNK_LENGTH, false, "".into())
    }

    #[test]
    fn test_chunks_right() {
        let chunks = generate_chunks_right();
        // in the right
        assert_eq!(
            chunks,
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
                    position: 20,
                    sequence: "AACCT".into()
                },
                ChunkedFasta {
                    position: 25,
                    sequence: "AACCT".into()
                }
            ]
        )
    }

    fn generate_indexes_left(genome: &str) -> RepeatPositions {
        let chunks = generate_chunks_left(genome);
        calculate_indexes(chunks, CHUNK_LENGTH, false, "test".into(), 0).unwrap()
    }

    #[test]
    fn test_index_left() {
        let indices = generate_indexes_left(GENOME);
        assert_eq!(
            indices.0,
            // i.e. there are 9 matches for AACCT
            vec![
                RepeatPosition {
                    id: "test".into(),
                    start: 0,
                    end: 10,
                    sequence: "AACCT".into()
                },
                RepeatPosition {
                    id: "test".into(),
                    start: 20,
                    end: 30,
                    sequence: "AACCT".into()
                }
            ]
        )
    }
    #[test]
    fn test_get_length_groups() {
        let indices = generate_indexes_left(GENOME_2);
        // we have AACCT 0-10, TAAAT 10-20, AACCT 20-30
        let map_len = indices.make_length_groups().get(&5).unwrap().len();
        assert_eq!(map_len, 3);
    }

    #[test]
    fn test_get_telomeric_repeat_estimates() {
        let mut indices = generate_indexes_left(GENOME_2);
        let res = get_telomeric_repeat_estimates(&mut indices).unwrap();
        assert_eq!(res, vec![("AACCT".to_string(), 4)]);
    }
}
