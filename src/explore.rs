use crate::{open_fasta_reader, utils, SubCommand};
use anyhow::bail;
use anyhow::Result;
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

// Scoring for error tolerant runs. Each base of a chunk that matches the
// seed scores +1, and each mismatch scores -MISMATCH_PENALTY, so runs extend
// through error rates below 1 / (1 + MISMATCH_PENALTY) = 25%.
const MISMATCH_PENALTY: i64 = 3;
// Cost of matching a different rotation of the seed than the previous chunk,
// i.e. an indel shifted the frame. A wrong kmer length (e.g. 5 on a 6bp
// telomere) changes rotation every chunk, so pays this every time and the run
// dies, rather than reporting a substring of the true unit.
const ROTATION_SWITCH_PENALTY: i64 = 8;

// Warn when one of the top units has less than this proportion of its
// copies on one strand. A telomeric repeat should read both ways round
// (both chromosome ends of an assembly, both strands across reads), so a
// strong bias in reads usually means strand-specific basecalling errors.
const STRAND_BIAS_WARNING: f64 = 0.1;
// how many of the top units to check for strand bias
const STRAND_BIAS_TOP_UNITS: usize = 5;

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
    let error_tolerant = !matches.get_flag("exact");
    if !error_tolerant {
        eprintln!("[+]\tOnly counting runs of exactly identical repeats");
    }

    // to report the telomeres...
    let mut output_vec: Vec<RepeatPositions> = Vec::new();
    // i.e. if you chose a length, as opposed to a minmum/maximum
    if length > 0 {
        eprintln!("[+]\tExploring genome for potential telomeric repeats of length: {length}");
        output_vec.append(&mut explore_length(
            input_fasta,
            length,
            dist_from_chromosome_end,
            verbose,
            threshold as usize,
            error_tolerant,
        )?);
    } else {
        // if a range was chosen.
        eprintln!(
            "[+]\tExploring genome for potential telomeric repeats between lengths {minimum} and {maximum}."
        );
        for length in minimum..maximum + 1 {
            eprintln!("[+]\t\tFinding telomeric repeat length: {length}");
            output_vec.append(&mut explore_length(
                input_fasta,
                length,
                dist_from_chromosome_end,
                verbose,
                threshold as usize,
                error_tolerant,
            )?);
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

    warn_strand_bias(&est);

    // copies of the unit in runs longer than the threshold
    println!("canonical_repeat_unit\tcopies\tcount_as_unit\tcount_as_revcomp");
    for e in est {
        let fmt = |c: Option<usize>| c.map_or("NA".to_string(), |c| c.to_string());
        println!(
            "{}\t{}\t{}\t{}",
            e.unit,
            e.count,
            fmt(e.as_unit),
            fmt(e.as_revcomp)
        );
    }

    // optional log file
    sc.log(matches)?;

    Ok(())
}

/// Scan every record in the fasta for tandem runs of a single
/// kmer length. Run coordinates are relative to the whole record,
/// so runs found at different lengths can be merged later.
fn explore_length(
    input_fasta: &PathBuf,
    length: usize,
    dist_from_chromosome_end: f64,
    verbose: bool,
    threshold: usize,
    error_tolerant: bool,
) -> Result<Vec<RepeatPositions>> {
    let reader = open_fasta_reader(input_fasta)?;

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
            // the right hand arm starts this far into the record
            let offsets = [0, seq_len - sequences[1].len()];

            for (sequence, offset) in sequences.into_iter().zip(offsets) {
                let positions = if error_tolerant {
                    tolerant_positions(&sequence, length, verbose, &id, threshold)
                } else {
                    let indexes = chunk_fasta(sequence, length, verbose, id.clone());
                    calculate_indexes(indexes, length, verbose, id.clone(), threshold)
                };

                if let Some(mut r) = positions {
                    r.shift(offset);
                    s.send(r).expect("Did not send!");
                }
            }
        });

    Ok(receiver.into_iter().collect())
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
                "[-]\tChunk length ({chunk_length}) greater than filtered sequence length ({sequence_len}) for {id}
[-]\tConsider increasing proportion of chromosome length covered. Skipping."
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
    for (a, b) in chunks.zip(chunks_plus_one) {
        if a == b {
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
    // move all positions along by `offset`
    fn shift(&mut self, offset: usize) {
        for el in &mut self.0 {
            el.start += offset;
            el.end += offset;
        }
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
    // the first run starts at the first index, not the start of the sequence
    let mut start = indexes.first().map_or(0, |c| c.position);
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
        } else if sequence1 == sequence2 && (position2 - position1) == sequence1.len() {
            continue;
        } else if !(sequence1 == sequence2 && (position2 - position1) == sequence1.len()) {
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
                "[-]\t\tChromosome {id}: No consecutive repeats of length {chunk_length} were identified."
            );
        }
        None
    } else {
        let filtered_repeat_positions = RepeatPositions(collection).filter_by_frequency(frequency);
        Some(filtered_repeat_positions)
    }
}

/// The error tolerant counterpart of [`chunk_fasta`] followed by
/// [`calculate_indexes`].
fn tolerant_positions(
    sequence: &[u8],
    chunk_length: usize,
    verbose: bool,
    id: &str,
    frequency: usize,
) -> Option<RepeatPositions> {
    let collection = tolerant_runs(sequence, chunk_length, id);
    if collection.is_empty() {
        if verbose {
            eprintln!(
                "[-]\t\tChromosome {id}: No consecutive repeats of length {chunk_length} were identified."
            );
        }
        None
    } else {
        Some(RepeatPositions(collection).filter_by_frequency(frequency))
    }
}

/// Find runs of a repeat which may contain sequencing errors.
///
/// A run is seeded where two adjacent chunks are identical, as in the exact
/// algorithm. It is then extended chunk by chunk, scoring each chunk against
/// the closest rotation of the seed (see [`MISMATCH_PENALTY`] and
/// [`ROTATION_SWITCH_PENALTY`]). Extension stops once the score drops more
/// than a fixed amount below its best, and the run is trimmed back to the
/// chunk where the score was highest.
fn tolerant_runs(sequence: &[u8], chunk_length: usize, id: &str) -> Vec<RepeatPosition> {
    let sequence = sequence.to_ascii_uppercase();
    let chunks: Vec<&[u8]> = sequence.chunks_exact(chunk_length).collect();
    // enough to survive one indel: a chunk spanning it, then the rotation switch
    let x_drop = 2 * chunk_length as i64 + ROTATION_SWITCH_PENALTY;

    let mut runs = Vec::new();
    let mut i = 0;
    while i + 1 < chunks.len() {
        if chunks[i] != chunks[i + 1] {
            i += 1;
            continue;
        }
        let seed = chunks[i];
        let rotations: Vec<Vec<u8>> = (0..chunk_length)
            .map(|r| [&seed[r..], &seed[..r]].concat())
            .collect();

        // the rotation each chunk matched, starting with the two seed chunks
        let mut matched = vec![0, 0];
        let mut score = 2 * chunk_length as i64;
        let mut best = score;
        let mut best_end = i + 2;
        for (j, chunk) in chunks.iter().enumerate().skip(i + 2) {
            let (r, chunk_score) = best_rotation(chunk, &rotations, *matched.last().unwrap());
            matched.push(r);
            score += chunk_score;
            if score > best {
                best = score;
                best_end = j + 1;
            }
            if best - score > x_drop {
                break;
            }
        }
        matched.truncate(best_end - i);

        runs.push(RepeatPosition {
            id: id.to_string(),
            start: i * chunk_length,
            end: best_end * chunk_length,
            sequence: consensus(&chunks[i..best_end], &matched),
        });
        i = best_end;
    }
    runs
}

/// The per position majority base across the chunks of a run, after undoing
/// the rotation each chunk matched. A seed can itself carry an error (both
/// seed chunks sharing it), so this, rather than the seed, labels the run.
/// Ties go to the seed's base, as the first chunk always has rotation 0.
fn consensus(chunks: &[&[u8]], rotations: &[usize]) -> String {
    let k = chunks[0].len();
    let mut counts: Vec<HashMap<u8, usize>> = vec![HashMap::new(); k];
    for (chunk, &r) in chunks.iter().zip(rotations) {
        // a chunk on rotation r has seed position (r + p) % k at position p
        for (p, &base) in chunk.iter().enumerate() {
            *counts[(r + p) % k].entry(base).or_insert(0) += 1;
        }
    }
    counts
        .iter()
        .enumerate()
        .map(|(p, c)| {
            let seed_base = chunks[0][p];
            let seed_count = c[&seed_base];
            let (&base, &count) = c.iter().max_by_key(|(_, &n)| n).unwrap();
            (if count > seed_count { base } else { seed_base }) as char
        })
        .collect()
}

/// Score a chunk against one rotation of the seed.
fn chunk_score(chunk: &[u8], rotation: &[u8]) -> i64 {
    let mismatches = chunk.iter().zip(rotation).filter(|(a, b)| a != b).count() as i64;
    chunk.len() as i64 - (1 + MISMATCH_PENALTY) * mismatches
}

/// The best scoring rotation of the seed for this chunk, preferring to stay
/// on the current rotation unless switching pays for itself.
fn best_rotation(chunk: &[u8], rotations: &[Vec<u8>], current: usize) -> (usize, i64) {
    let mut best = (current, chunk_score(chunk, &rotations[current]));
    for (r, rotation) in rotations.iter().enumerate() {
        if r == current {
            continue;
        }
        let score = chunk_score(chunk, rotation) - ROTATION_SWITCH_PENALTY;
        if score > best.1 {
            best = (r, score);
        }
    }
    best
}

/// check if a sequence looks like it is not
/// a telomeric repeat
fn check_telomeric_repeat(sequence: &str) -> bool {
    let repeat_period = check_repeats(sequence);
    repeat_period < REPEAT_PERIOD_THRESHOLD
}

/// Which way round a run reads, relative to its canonical unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
enum Orientation {
    /// A rotation of the unit itself, e.g. CCCTAA for AACCCT.
    Unit,
    /// A rotation of the unit's reverse complement, e.g. TTAGGG for AACCCT.
    Revcomp,
    /// The unit is a rotation of its own reverse complement (e.g. AATT), so
    /// there is no way to tell.
    Either,
}

/// The orientation of a run's sequence relative to its primitive canonical
/// unit (see [`utils::primitive_telomere_unit()`]).
fn orientation(sequence: &str, unit: &str) -> Orientation {
    if utils::string_rotation(unit, &utils::reverse_complement(unit)) {
        return Orientation::Either;
    }
    // the sequence may be a multimer of the unit, but any unit length
    // window of it is a rotation of the unit or its reverse complement
    if utils::string_rotation(&sequence[..unit.len()], unit) {
        Orientation::Unit
    } else {
        Orientation::Revcomp
    }
}

/// Runs are grouped by unit, sequence id and orientation.
type RunKey<'a> = (String, &'a str, Orientation);

/// A candidate telomeric repeat unit and how many copies of it were found.
#[derive(Debug, PartialEq, Eq)]
struct UnitEstimate {
    unit: String,
    count: usize,
    /// Copies reading as the unit, None if orientation can't be told.
    as_unit: Option<usize>,
    /// Copies reading as its reverse complement, None if orientation can't be told.
    as_revcomp: Option<usize>,
}

impl UnitEstimate {
    /// The proportion of copies on the less common strand.
    fn minor_strand_proportion(&self) -> Option<f64> {
        let (u, r) = (self.as_unit?, self.as_revcomp?);
        if u + r == 0 {
            return None;
        }
        Some(u.min(r) as f64 / (u + r) as f64)
    }
}

/// Takes the final aggregation of potential telomeric repeats across
/// chromosomes and also potentially across different lengths and tries
/// to find the most likely telomeric repeat.
///
/// Each run is reduced to its primitive canonical unit (see
/// [`utils::primitive_telomere_unit()`]), so rotations, reverse complements
/// and exact multimers (e.g. AACCTAACCT) of a repeat all count towards the
/// same unit. Runs of the same unit on the same sequence are then merged where
/// they overlap, so a telomere found at several kmer lengths is only counted
/// once. The count is the number of copies of the unit in the merged runs,
/// also split by which way round the runs read (see [`Orientation`]).
fn get_telomeric_repeat_estimates(
    telomeric_repeats: &mut RepeatPositions,
) -> Result<Vec<UnitEstimate>> {
    // (unit, sequence id, orientation) -> run intervals. A stretch of sequence
    // only reads one way round, so splitting on orientation never splits a
    // telomere found at several kmer lengths.
    let mut runs: HashMap<RunKey, Vec<(usize, usize)>> = HashMap::new();
    for el in &telomeric_repeats.0 {
        let unit = utils::primitive_telomere_unit(&el.sequence);
        let orientation = orientation(&el.sequence, &unit);
        runs.entry((unit, el.id.as_str(), orientation))
            .or_default()
            .push((el.start, el.end));
    }

    let mut map: HashMap<String, (usize, usize, usize)> = HashMap::new();
    for ((unit, _, orientation), mut intervals) in runs {
        let copies = merged_length(&mut intervals) / unit.len();
        let counts = map.entry(unit).or_insert((0, 0, 0));
        match orientation {
            Orientation::Unit => counts.0 += copies,
            Orientation::Revcomp => counts.1 += copies,
            Orientation::Either => counts.2 += copies,
        }
    }

    let mut estimates: Vec<_> = map
        .into_iter()
        .map(|(unit, (as_unit, as_revcomp, either))| {
            let palindromic = either > 0;
            UnitEstimate {
                unit,
                count: as_unit + as_revcomp + either,
                as_unit: (!palindromic).then_some(as_unit),
                as_revcomp: (!palindromic).then_some(as_revcomp),
            }
        })
        .collect();
    // ties broken on the unit so output is deterministic
    estimates.sort_by(|a, b| b.count.cmp(&a.count).then_with(|| a.unit.cmp(&b.unit)));
    filter_count_vec(&mut estimates)?;

    Ok(estimates)
}

/// Print a warning for any of the top units which mostly read one way round.
fn warn_strand_bias(estimates: &[UnitEstimate]) {
    for e in estimates.iter().take(STRAND_BIAS_TOP_UNITS) {
        if let Some(p) = e.minor_strand_proportion() {
            if p < STRAND_BIAS_WARNING {
                eprintln!(
                    "[!]\t{}: only {:.1}% of copies are on the minor strand. In reads this usually means strand-specific basecalling errors: the unit may be an artefact, or a real repeat basecalled badly on one strand.",
                    e.unit,
                    p * 100.0
                );
            }
        }
    }
}

/// Total length covered by a set of half open intervals,
/// counting overlapping stretches once.
fn merged_length(intervals: &mut [(usize, usize)]) -> usize {
    intervals.sort_unstable();
    let mut total = 0;
    let mut current: Option<(usize, usize)> = None;
    for &(start, end) in intervals.iter() {
        match current {
            Some((cs, ce)) if start <= ce => current = Some((cs, ce.max(end))),
            _ => {
                if let Some((cs, ce)) = current {
                    total += ce - cs;
                }
                current = Some((start, end));
            }
        }
    }
    if let Some((cs, ce)) = current {
        total += ce - cs;
    }
    total
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
fn filter_count_vec(v: &mut Vec<UnitEstimate>) -> Result<()> {
    // monomers
    // not sure I need this.
    v.retain(|e| {
        let repeat_period = check_repeats(&e.unit);
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

    /// (unit, count) pairs from the estimates
    fn unit_counts(estimates: Vec<UnitEstimate>) -> Vec<(String, usize)> {
        estimates.into_iter().map(|e| (e.unit, e.count)).collect()
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
    fn test_get_telomeric_repeat_estimates() {
        let mut indices = generate_indexes_left(GENOME_2);
        let res = unit_counts(get_telomeric_repeat_estimates(&mut indices).unwrap());
        // AACCT 0-10 and 20-30, TAAAT (canonical AAATT) 10-20
        assert_eq!(
            res,
            vec![("AACCT".to_string(), 4), ("AAATT".to_string(), 2)]
        );
    }

    fn run(id: &str, start: usize, end: usize, sequence: &str) -> RepeatPosition {
        RepeatPosition {
            id: id.into(),
            start,
            end,
            sequence: sequence.into(),
        }
    }

    #[test]
    fn test_estimates_collapse_multimers() {
        // a 7-mer run, and a separate run seen as its 14-mer (rotated) multimer
        let mut positions = RepeatPositions(vec![
            run("chr1", 0, 70, "TTTAGGG"),
            run("chr1", 1000, 1140, "AGGGTTTAGGGTTT"),
        ]);
        let res = unit_counts(get_telomeric_repeat_estimates(&mut positions).unwrap());
        assert_eq!(res, vec![("AAACCCT".to_string(), 30)]);
    }

    #[test]
    fn test_estimates_do_not_double_count_across_lengths() {
        // the same telomere found at k = 5 and k = 10
        let mut positions = RepeatPositions(vec![
            run("chr1", 0, 500, "TTAGG"),
            run("chr1", 0, 500, "TTAGGTTAGG"),
            // same coordinates on another sequence are a different telomere
            run("chr2", 0, 500, "TTAGG"),
        ]);
        let res = unit_counts(get_telomeric_repeat_estimates(&mut positions).unwrap());
        assert_eq!(res, vec![("AACCT".to_string(), 200)]);
    }

    #[test]
    fn test_merged_length() {
        let mut intervals = vec![(20, 30), (0, 10), (5, 15), (30, 40)];
        assert_eq!(merged_length(&mut intervals), 35);
    }

    #[test]
    fn test_right_arm_offset() {
        let mut positions = generate_indexes_left(GENOME);
        positions.shift(30);
        assert_eq!(positions.0[0].start, 30);
        assert_eq!(positions.0[1].end, 60);
    }

    #[test]
    fn test_first_run_starts_at_first_index() {
        // junk before the first run must not be included in it
        let seq = format!("{}{}", "GATTC", "AACCT".repeat(4));
        let chunks = chunk_fasta(seq.into_bytes(), CHUNK_LENGTH, false, "".into());
        let positions = calculate_indexes(chunks, CHUNK_LENGTH, false, "test".into(), 0).unwrap();
        assert_eq!(positions.0, vec![run("test", 5, 25, "AACCT")]);
    }

    /// a small deterministic pseudo random sequence
    fn random_seq(len: usize, mut state: u64) -> String {
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 62) as usize] as char
            })
            .collect()
    }

    fn tolerant(seq: &str, k: usize) -> Vec<RepeatPosition> {
        tolerant_runs(seq.as_bytes(), k, "test")
    }

    #[test]
    fn test_tolerant_matches_exact_on_clean_repeat() {
        let seq = "TTAGGG".repeat(50);
        assert_eq!(tolerant(&seq, 6), vec![run("test", 0, 300, "TTAGGG")]);
    }

    #[test]
    fn test_tolerant_spans_substitutions() {
        // a substitution every 20 copies breaks exact runs into short pieces
        let mut seq = "TTAGGG".repeat(100).into_bytes();
        for i in (60..600).step_by(120) {
            seq[i + 3] = b'C';
        }
        let seq = String::from_utf8(seq).unwrap();
        assert_eq!(tolerant(&seq, 6), vec![run("test", 0, 600, "TTAGGG")]);
    }

    #[test]
    fn test_tolerant_spans_indels() {
        // a deletion and, later, an insertion
        let seq = format!(
            "{}{}{}{}{}",
            "TTAGGG".repeat(30),
            "TAGGG",
            "TTAGGG".repeat(30),
            "TTAAGGG",
            "TTAGGG".repeat(30)
        );
        let runs = tolerant(&seq, 6);
        assert_eq!(runs.len(), 1);
        assert_eq!(runs[0].start, 0);
        assert!(seq.len() - runs[0].end < 6);
    }

    #[test]
    fn test_tolerant_trims_trailing_junk() {
        let seq = format!("{}{}", "TTAGGG".repeat(50), random_seq(300, 7));
        assert_eq!(tolerant(&seq, 6), vec![run("test", 0, 300, "TTAGGG")]);
    }

    #[test]
    fn test_tolerant_rejects_wrong_length() {
        // 5-mers of a 6bp telomere must not extend into long runs
        let seq = "TTAGGG".repeat(200);
        assert!(tolerant(&seq, 5).iter().all(|r| r.get_count() < 10));
    }

    #[test]
    fn test_tolerant_random_sequence_has_no_long_runs() {
        let seq = random_seq(100_000, 42);
        for k in 5..=12 {
            assert!(tolerant(&seq, k).iter().all(|r| r.get_count() < 10));
        }
    }

    #[test]
    fn test_tolerant_labels_run_with_consensus() {
        // both seed chunks carry the same error, the rest of the run does not
        let seq = format!("{}{}", "TTCGGG".repeat(2), "TTAGGG".repeat(40));
        assert_eq!(tolerant(&seq, 6), vec![run("test", 0, 252, "TTAGGG")]);
    }

    #[test]
    fn test_consensus_undoes_rotation() {
        let chunks: Vec<&[u8]> = vec![b"TTAGGG", b"TTAGGG", b"AGGGTT", b"AGGGTT"];
        assert_eq!(consensus(&chunks, &[0, 0, 2, 2]), "TTAGGG");
    }

    #[test]
    fn test_orientation() {
        assert_eq!(orientation("CCCTAA", "AACCCT"), Orientation::Unit);
        assert_eq!(orientation("TTAGGG", "AACCCT"), Orientation::Revcomp);
        // multimers, rotated
        assert_eq!(orientation("CTAACCCTAACC", "AACCCT"), Orientation::Unit);
        assert_eq!(orientation("GGGTTAGGGTTA", "AACCCT"), Orientation::Revcomp);
        // AATT is its own reverse complement
        assert_eq!(orientation("ATTA", "AATT"), Orientation::Either);
    }

    #[test]
    fn test_estimates_split_by_strand() {
        let mut positions = RepeatPositions(vec![
            // start of chr1, C-rich strand
            run("chr1", 0, 600, "CCCTAA"),
            // end of chr1 at two kmer lengths, G-rich strand
            run("chr1", 10_000, 10_300, "TTAGGG"),
            run("chr1", 10_000, 10_300, "GGGTTAGGGTTA"),
        ]);
        let res = get_telomeric_repeat_estimates(&mut positions).unwrap();
        assert_eq!(
            res,
            vec![UnitEstimate {
                unit: "AACCCT".into(),
                count: 150,
                as_unit: Some(100),
                as_revcomp: Some(50),
            }]
        );
        assert_eq!(res[0].minor_strand_proportion(), Some(50.0 / 150.0));
    }

    #[test]
    fn test_estimates_palindromic_unit_has_no_strand() {
        let mut positions = RepeatPositions(vec![run("chr1", 0, 400, "AATTAATT")]);
        let res = get_telomeric_repeat_estimates(&mut positions).unwrap();
        assert_eq!(res[0].as_unit, None);
        assert_eq!(res[0].minor_strand_proportion(), None);
    }
}
