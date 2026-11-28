use crate::{open_fasta, utils, SubCommand};
use anyhow::bail;
use anyhow::Result;
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::path::PathBuf;
use std::str;
use std::sync::mpsc::channel;
use std::sync::Mutex;

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

    // i.e. if you chose a length, as opposed to a minmum/maximum
    let lengths: Vec<usize> = if length > 0 {
        eprintln!("[+]\tExploring genome for potential telomeric repeats of length: {length}");
        vec![length]
    } else {
        // if a range was chosen.
        eprintln!(
            "[+]\tExploring genome for potential telomeric repeats between lengths {minimum} and {maximum}."
        );
        (minimum..maximum + 1).collect()
    };
    // All lengths are explored in one pass over the input. Progress and
    // verbose messages come out as if each length had its own pass: those of
    // the first length as they happen, the others after the pass, each
    // after its own "Finding" line.
    let range = length == 0;
    let finding = |l: usize| {
        if range {
            eprintln!("[+]\t\tFinding telomeric repeat length: {l}");
        }
    };
    let mut output_vec: Vec<RepeatPositions> = Vec::new();
    if let Some(&first) = lengths.first() {
        finding(first);
        let (mut positions, messages) = explore_lengths(
            input_fasta,
            &lengths,
            dist_from_chromosome_end,
            verbose,
            threshold as usize,
            error_tolerant,
        )?;
        output_vec.append(&mut positions);
        for (&l, messages) in lengths.iter().zip(messages).skip(1) {
            finding(l);
            for m in messages {
                eprintln!("{m}");
            }
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
    println!("canonical_repeat_unit\tcopies\tcopies_as_unit\tcopies_as_revcomp");
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

/// Where verbose messages go: printed, or kept to print later.
type Log<'a> = Option<&'a (dyn Fn(String) + Sync)>;

/// Scan every record in the fasta for tandem runs of each kmer length, in
/// one pass. Run coordinates are relative to the whole record, so runs found
/// at different lengths can be merged later. Verbose messages for the first
/// length are printed as they happen; those for the others are returned, per
/// length, to print later.
fn explore_lengths(
    input_fasta: &PathBuf,
    lengths: &[usize],
    dist_from_chromosome_end: f64,
    verbose: bool,
    threshold: usize,
    error_tolerant: bool,
) -> Result<(Vec<RepeatPositions>, Vec<Vec<String>>)> {
    let reader = open_fasta(input_fasta)?;
    let kept: Vec<Mutex<Vec<String>>> = lengths.iter().map(|_| Mutex::new(Vec::new())).collect();
    let print = |m: String| eprintln!("{m}");
    let keep: Vec<Box<dyn Fn(String) + Sync + '_>> = kept
        .iter()
        .map(|k| Box::new(move |m: String| k.lock().unwrap().push(m)) as Box<dyn Fn(String) + Sync>)
        .collect();
    let log = |i: usize| -> Log {
        if !verbose {
            None
        } else if i == 0 {
            Some(&print)
        } else {
            Some(keep[i].as_ref())
        }
    };

    let (sender, receiver) = channel();
    reader
        .records()
        .par_bridge()
        .for_each_with(sender, |s, record| {
            let record = record.expect("[-]\tError during fasta record parsing.");
            let id = record.id();
            let seq_len = record.seq().len();

            let sequences = split_seq_by_distance(&record, dist_from_chromosome_end, seq_len);
            // the right hand arm starts this far into the record
            let offsets = [0, seq_len - sequences[1].len()];

            // each length, then each arm, as separate passes would have
            let tasks: Vec<(usize, usize)> = (0..lengths.len())
                .flat_map(|i| (0..2).map(move |a| (i, a)))
                .collect();
            let found: Vec<RepeatPositions> = tasks
                .into_par_iter()
                .filter_map(|(i, a)| {
                    let (sequence, length) = (&sequences[a], lengths[i]);
                    let positions = if error_tolerant {
                        tolerant_positions(sequence, length, log(i), id, threshold)
                    } else {
                        exact_positions(sequence, length, log(i), id, threshold)
                    };
                    positions.map(|mut r| {
                        r.shift(offsets[a]);
                        r
                    })
                })
                .collect();
            for r in found {
                s.send(r).expect("Did not send!");
            }
        });

    drop(keep);
    let messages = kept.into_iter().map(|k| k.into_inner().unwrap()).collect();
    Ok((receiver.into_iter().collect(), messages))
}

pub fn split_seq_by_distance(
    sequence: &crate::fasta::Record,
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
    sequence: &[u8],
    chunk_length: usize,
    log: Log,
    id: String,
) -> Vec<ChunkedFasta> {
    let sequence_len = sequence.len();
    let chunks = sequence.chunks(chunk_length);
    // catch edge cases where chunk length greater than sequence length.
    if sequence_len <= chunk_length {
        if let Some(log) = log {
            log(format!(
                "[-]\tChunk length ({chunk_length}) greater than filtered sequence length ({sequence_len}) for {id}
[-]\tConsider increasing proportion of chromosome length covered. Skipping."
            ));
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
        Self(
            std::mem::take(&mut self.0)
                .into_iter()
                .filter(|e| e.get_count() > frequency && !e.is_simple_repeat())
                .collect(),
        )
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
    log: Log,
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
        if let Some(log) = log {
            log(format!(
                "[-]\t\tChromosome {id}: No consecutive repeats of length {chunk_length} were identified."
            ));
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
    log: Log,
    id: &str,
    frequency: usize,
) -> Option<RepeatPositions> {
    let (any, kept) = filtered_runs(sequence, chunk_length, id, true, frequency);
    if !any {
        if let Some(log) = log {
            log(format!(
                "[-]\t\tChromosome {id}: No consecutive repeats of length {chunk_length} were identified."
            ));
        }
        None
    } else {
        Some(RepeatPositions(kept))
    }
}

/// [`chunk_fasta`] followed by [`calculate_indexes`], with the same messages.
fn exact_positions(
    sequence: &[u8],
    chunk_length: usize,
    log: Log,
    id: &str,
    frequency: usize,
) -> Option<RepeatPositions> {
    if chunk_length == 0 || !sequence.is_ascii() {
        let indexes = chunk_fasta(sequence, chunk_length, log, id.to_string());
        return calculate_indexes(indexes, chunk_length, log, id.to_string(), frequency);
    }
    let sequence_len = sequence.len();
    if sequence_len <= chunk_length {
        if let Some(log) = log {
            log(format!(
                "[-]\tChunk length ({chunk_length}) greater than filtered sequence length ({sequence_len}) for {id}
[-]\tConsider increasing proportion of chromosome length covered. Skipping."
            ));
        }
    }
    let (any, kept) = filtered_runs(sequence, chunk_length, id, false, frequency);
    if !any {
        if let Some(log) = log {
            log(format!(
                "[-]\t\tChromosome {id}: No consecutive repeats of length {chunk_length} were identified."
            ));
        }
        None
    } else {
        Some(RepeatPositions(kept))
    }
}

/// Whether there are any runs at all, and the runs passing
/// [`RepeatPositions::filter_by_frequency`], allocating only for those.
fn filtered_runs(
    sequence: &[u8],
    chunk_length: usize,
    id: &str,
    error_tolerant: bool,
    frequency: usize,
) -> (bool, Vec<RepeatPosition>) {
    let k = chunk_length;
    if !error_tolerant && (k == 0 || !sequence.is_ascii()) {
        // the reference, whose exact runs are already without simple repeats
        let runs = find_runs_reference(sequence, k, id, false);
        let any = !runs.is_empty();
        return (any, RepeatPositions(runs).filter_by_frequency(frequency).0);
    }
    let mut any = false;
    let mut kept = Vec::new();
    // as RepeatPosition::get_count and is_simple_repeat
    let mut take = |start: usize, end: usize, unit: &[u8]| {
        any = true;
        if (end - start) / unit.len() > frequency && smallest_period(unit) >= REPEAT_PERIOD_THRESHOLD
        {
            kept.push(RepeatPosition {
                id: id.to_string(),
                start,
                end,
                sequence: String::from_utf8_lossy(unit).into_owned(),
            });
        }
    };
    if error_tolerant {
        scan_runs(sequence, k, true, &mut take);
    } else {
        // every run, simple repeats too, for whether there are any
        exact_scan(sequence, k, |start, end, unit| take(start * k, end * k, unit));
    }
    (any, kept)
}

/// Find runs of a repeat which may contain sequencing errors.
///
/// A run is seeded where two adjacent chunks are identical, as in the exact
/// algorithm. It is then extended chunk by chunk, scoring each chunk against
/// the closest rotation of the seed (see [`MISMATCH_PENALTY`] and
/// [`ROTATION_SWITCH_PENALTY`]). Extension stops once the score drops more
/// than a fixed amount below its best, and the run is trimmed back to the
/// chunk where the score was highest.
#[cfg(test)]
fn tolerant_runs(sequence: &[u8], chunk_length: usize, id: &str) -> Vec<RepeatPosition> {
    let mut runs = Vec::new();
    scan_runs(sequence, chunk_length, true, |start, end, unit| {
        runs.push(RepeatPosition {
            id: id.to_string(),
            start,
            end,
            sequence: String::from_utf8_lossy(unit).into_owned(),
        })
    });
    runs
}

/// [`tolerant_runs`] as first written: the reference for the fast scanner,
/// and used for non-ASCII sequence.
fn tolerant_runs_reference(sequence: &[u8], chunk_length: usize, id: &str) -> Vec<RepeatPosition> {
    let sequence = sequence.to_ascii_uppercase();
    let chunks: Vec<&[u8]> = sequence.chunks_exact(chunk_length).collect();
    // enough to survive one indel: a chunk spanning it, then the rotation switch
    let x_drop = 2 * chunk_length as i64 + ROTATION_SWITCH_PENALTY;

    let mut runs = Vec::new();
    // chunks before this belong to the previous run
    let mut previous_end = 0;
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

        // extend right from the seed, then left back to the previous run, so
        // a run is found the same whichever way round the sequence reads
        let (forward, end) = extend(&chunks, (i + 2..chunks.len()).collect(), &rotations, x_drop);
        let (backward, start) = extend(
            &chunks,
            (previous_end..i).rev().collect(),
            &rotations,
            x_drop,
        );
        let start = start.unwrap_or(i);
        let end = end.map_or(i + 2, |e| e + 1);

        // the rotation each chunk matched, in order
        let matched: Vec<usize> = backward
            .into_iter()
            .rev()
            .chain([0, 0])
            .chain(forward)
            .collect();

        runs.push(RepeatPosition {
            id: id.to_string(),
            start: start * chunk_length,
            end: end * chunk_length,
            sequence: consensus(&chunks[start..end], &matched, seed),
        });
        previous_end = end;
        i = end;
    }
    runs
}

/// Extend a run from its seed over the chunks at `indices`, in order, until
/// the score drops more than `x_drop` below its best. Returns the rotation
/// each chunk up to the best point matched, and the index of the last chunk
/// in the run (None if the run doesn't extend at all).
fn extend(
    chunks: &[&[u8]],
    indices: Vec<usize>,
    rotations: &[Vec<u8>],
    x_drop: i64,
) -> (Vec<usize>, Option<usize>) {
    let mut matched = Vec::new();
    let mut rotation = 0;
    let mut score = 0;
    let mut best = 0;
    let mut best_len = 0;
    for &j in &indices {
        let (r, chunk_score) = best_rotation(chunks[j], rotations, rotation);
        rotation = r;
        matched.push(r);
        score += chunk_score;
        if score > best {
            best = score;
            best_len = matched.len();
        }
        if best - score > x_drop {
            break;
        }
    }
    matched.truncate(best_len);
    let last = best_len.checked_sub(1).map(|n| indices[n]);
    (matched, last)
}

/// The per position majority base across the chunks of a run, after undoing
/// the rotation each chunk matched. A seed can itself carry an error (both
/// seed chunks sharing it), so this, rather than the seed, labels the run.
/// Ties go to the seed's base.
fn consensus(chunks: &[&[u8]], rotations: &[usize], seed: &[u8]) -> String {
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
            let seed_base = seed[p];
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

/// All runs of repeats of one kmer length in a sequence, without filtering
/// on run length or repeat type. Positions are relative to the sequence.
#[cfg(test)]
pub(crate) fn find_runs(
    sequence: &[u8],
    chunk_length: usize,
    id: &str,
    error_tolerant: bool,
) -> Vec<RepeatPosition> {
    let mut runs = Vec::new();
    scan_runs(sequence, chunk_length, error_tolerant, |start, end, unit| {
        runs.push(RepeatPosition {
            id: id.to_string(),
            start,
            end,
            sequence: String::from_utf8_lossy(unit).into_owned(),
        })
    });
    runs
}

/// [`find_runs`] as first written, the reference for [`scan_runs`].
fn find_runs_reference(
    sequence: &[u8],
    chunk_length: usize,
    id: &str,
    error_tolerant: bool,
) -> Vec<RepeatPosition> {
    if error_tolerant {
        tolerant_runs_reference(sequence, chunk_length, id)
    } else {
        let indexes = chunk_fasta(sequence, chunk_length, None, id.to_string());
        calculate_indexes(indexes, chunk_length, None, id.to_string(), 0)
            .map(|r| r.0)
            .unwrap_or_default()
    }
}

/// The runs [`find_runs`] finds, in order, passed to `emit` as (start, end,
/// unit as the run reads) without allocating for each. ASCII sequence is
/// scanned by [`tolerant_scan`] or [`exact_scan`]; other sequence (and a
/// kmer length of 0, which panics as before) by the reference code.
pub(crate) fn scan_runs(
    sequence: &[u8],
    chunk_length: usize,
    error_tolerant: bool,
    mut emit: impl FnMut(usize, usize, &[u8]),
) {
    if chunk_length == 0 || !sequence.is_ascii() {
        for r in find_runs_reference(sequence, chunk_length, "", error_tolerant) {
            emit(r.start, r.end, r.sequence.as_bytes());
        }
        return;
    }
    let k = chunk_length;
    if error_tolerant {
        let upper = ascii_uppercase(sequence);
        tolerant_scan(&upper, k, |start, end, unit| emit(start * k, end * k, unit));
    } else {
        // calculate_indexes drops simple repeats, even at frequency 0
        exact_scan(sequence, k, |start, end, unit| {
            if smallest_period(unit) >= REPEAT_PERIOD_THRESHOLD {
                emit(start * k, end * k, unit)
            }
        });
    }
}

/// An uppercase view of ASCII sequence, copied only if it has lowercase.
fn ascii_uppercase(sequence: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    if sequence.iter().any(u8::is_ascii_lowercase) {
        std::borrow::Cow::Owned(sequence.to_ascii_uppercase())
    } else {
        std::borrow::Cow::Borrowed(sequence)
    }
}

/// The smallest period of a string, as [`check_repeats`] finds it for ASCII:
/// the smallest d with s[i] == s[i - d] for all i >= d (s.len() if none).
fn smallest_period(s: &[u8]) -> usize {
    if s.is_empty() || !s.is_ascii() {
        return check_repeats(&String::from_utf8_lossy(s));
    }
    (1..s.len())
        .find(|&d| s[d..] == s[..s.len() - d])
        .unwrap_or(s.len())
}

/// The exact runs of [`chunk_fasta`] and [`calculate_indexes`] in ASCII
/// sequence, in chunks: two or more identical consecutive chunks of `k` bp
/// (compared as in the sequence, case and all), with consecutive runs
/// merged when they are the same once uppercased. `emit` gets (first chunk,
/// chunk after the last, uppercase unit).
fn exact_scan(seq: &[u8], k: usize, mut emit: impl FnMut(usize, usize, &[u8])) {
    let n = seq.len() / k;
    let chunk = |i: usize| &seq[i * k..(i + 1) * k];
    let same = |i: usize| {
        let (a, b) = (chunk(i), chunk(i + 1));
        a[0] == b[0] && a == b
    };
    let mut current: Option<(usize, usize)> = None;
    let mut unit: Vec<u8> = Vec::with_capacity(k);
    let mut next: Vec<u8> = Vec::with_capacity(k);
    let mut i = 0;
    while i + 1 < n {
        if !same(i) {
            i += 1;
            continue;
        }
        let mut j = i + 1;
        while j + 1 < n && same(j) {
            j += 1;
        }
        // a run of identical chunks i..=j
        next.clear();
        next.extend(chunk(i).iter().map(u8::to_ascii_uppercase));
        current = match current {
            Some((start, end)) if end == i && next == unit => Some((start, j + 1)),
            other => {
                if let Some((start, end)) = other {
                    emit(start, end, &unit);
                }
                std::mem::swap(&mut unit, &mut next);
                Some((i, j + 1))
            }
        };
        i = j + 1;
    }
    if let Some((start, end)) = current {
        emit(start, end, &unit);
    }
}

/// k-mers of up to this length are compared as packed integers.
const PACKED_MAX: usize = 16;

/// Up to 16 bytes as an integer, the rest zero.
fn pack(bytes: &[u8]) -> u128 {
    let mut b = [0u8; 16];
    b[..bytes.len()].copy_from_slice(bytes);
    u128::from_le_bytes(b)
}

/// Bytes that differ between two packed k-mers.
fn packed_mismatches(a: u128, b: u128) -> i64 {
    const LOW: u128 = u128::from_le_bytes([0x7f; 16]);
    const HIGH: u128 = u128::from_le_bytes([0x80; 16]);
    let x = a ^ b;
    ((((x & LOW) + LOW) | x) & HIGH).count_ones() as i64
}

/// [`tolerant_runs_reference`] on uppercase ASCII sequence, in chunks,
/// without allocating per seed. `emit` gets (first chunk, chunk after the
/// last, consensus unit).
fn tolerant_scan(seq: &[u8], k: usize, mut emit: impl FnMut(usize, usize, &[u8])) {
    let n = seq.len() / k;
    let chunk = |i: usize| &seq[i * k..(i + 1) * k];
    let packed = k <= PACKED_MAX;
    // packed chunks, for the comparisons with rotations of the seed
    let pchunk = |i: usize| -> u128 {
        if i * k + 16 <= seq.len() {
            let mask = if k == 16 { u128::MAX } else { (1u128 << (8 * k)) - 1 };
            u128::from_le_bytes(seq[i * k..i * k + 16].try_into().unwrap()) & mask
        } else {
            pack(chunk(i))
        }
    };
    let x_drop = 2 * k as i64 + ROTATION_SWITCH_PENALTY;
    // rotation r of the seed at [r * k, (r + 1) * k)
    let mut rotations = vec![0u8; k * k];
    let mut packed_rotations = vec![0u128; k];
    let (mut forward, mut backward) = (Vec::new(), Vec::new());
    let mut unit = vec![0u8; k];

    // as `best_rotation`: stay on the current rotation unless another,
    // paying the switch penalty, scores higher; the first such on ties
    let best_rotation = |j: usize, rotations: &[u8], packed_rotations: &[u128], current: usize| {
        let score = |r: usize| -> i64 {
            let mismatches = if packed {
                packed_mismatches(pchunk(j), packed_rotations[r])
            } else {
                chunk(j)
                    .iter()
                    .zip(&rotations[r * k..(r + 1) * k])
                    .filter(|(a, b)| a != b)
                    .count() as i64
            };
            k as i64 - (1 + MISMATCH_PENALTY) * mismatches
        };
        let mut best = (current, score(current));
        for r in 0..k {
            if r == current {
                continue;
            }
            let s = score(r) - ROTATION_SWITCH_PENALTY;
            if s > best.1 {
                best = (r, s);
            }
        }
        best
    };
    // as `extend`: the rotations matched up to the best score, and the last
    // chunk of the run in this direction
    let extend = |indices: &mut dyn Iterator<Item = usize>,
                  matched: &mut Vec<usize>,
                  rotations: &[u8],
                  packed_rotations: &[u128]|
     -> Option<usize> {
        matched.clear();
        let (mut rotation, mut score, mut best, mut best_len) = (0, 0i64, 0i64, 0usize);
        let mut last = None;
        for j in indices {
            let (r, chunk_score) = best_rotation(j, rotations, packed_rotations, rotation);
            rotation = r;
            matched.push(r);
            score += chunk_score;
            if score > best {
                best = score;
                best_len = matched.len();
                last = Some(j);
            }
            if best - score > x_drop {
                break;
            }
        }
        matched.truncate(best_len);
        last
    };

    let mut previous_end = 0;
    let mut i = 0;
    while i + 1 < n {
        let (a, b) = (chunk(i), chunk(i + 1));
        if a[0] != b[0] || a != b {
            i += 1;
            continue;
        }
        for r in 0..k {
            let rot = &mut rotations[r * k..(r + 1) * k];
            rot[..k - r].copy_from_slice(&a[r..]);
            rot[k - r..].copy_from_slice(&a[..r]);
            if packed {
                packed_rotations[r] = pack(rot);
            }
        }
        let end = extend(&mut (i + 2..n), &mut forward, &rotations, &packed_rotations);
        let start = extend(
            &mut (previous_end..i).rev(),
            &mut backward,
            &rotations,
            &packed_rotations,
        );
        let start = start.unwrap_or(i);
        let end = end.map_or(i + 2, |e| e + 1);

        // as `consensus`: per position of the seed, the commonest base over
        // the run's chunks (undoing each chunk's rotation), unless it is no
        // commoner than the seed's base
        let rotation_of = |c: usize| -> usize {
            if c < i {
                backward[i - 1 - c]
            } else if c < i + 2 {
                0
            } else {
                forward[c - i - 2]
            }
        };
        let mut counts = [0u32; 128];
        for (p, u) in unit.iter_mut().enumerate() {
            for c in start..end {
                let r = rotation_of(c);
                counts[chunk(c)[(p + k - r) % k] as usize] += 1;
            }
            let seed_base = a[p];
            let seed_count = counts[seed_base as usize];
            // the commonest base, lowest byte on ties (rust-bio's HashMap
            // order made ties above the seed's count random)
            let (base, count) = counts
                .iter()
                .enumerate()
                .fold((0usize, 0u32), |best, (b, &n)| if n > best.1 { (b, n) } else { best });
            *u = if count > seed_count { base as u8 } else { seed_base };
            for c in start..end {
                let r = rotation_of(c);
                counts[chunk(c)[(p + k - r) % k] as usize] = 0;
            }
        }
        emit(start, end, &unit);
        previous_end = end;
        i = end;
    }
}

/// The most abundant candidate telomeric repeat unit across a set of named
/// sequences, as `tidk explore` would report it, along with the proportion
/// of its copies on the minor strand. Names must be unique, as runs are
/// merged by name.
pub(crate) fn top_unit(
    sequences: &[(String, &[u8])],
    minimum: usize,
    maximum: usize,
    threshold: usize,
    error_tolerant: bool,
) -> Result<Option<(String, Option<f64>)>> {
    Ok(top_units(sequences, minimum, maximum, threshold, error_tolerant)?
        .into_iter()
        .next())
}

/// As [`top_unit()`], every candidate unit, most abundant first.
pub(crate) fn top_units(
    sequences: &[(String, &[u8])],
    minimum: usize,
    maximum: usize,
    threshold: usize,
    error_tolerant: bool,
) -> Result<Vec<(String, Option<f64>)>> {
    let mut positions = RepeatPositions::new();
    for (name, sequence) in sequences {
        for length in minimum..=maximum {
            // find_runs, then filter_by_frequency
            positions.add(&mut filtered_runs(sequence, length, name, error_tolerant, threshold).1);
        }
    }
    let estimates = get_telomeric_repeat_estimates(&mut positions)?;
    Ok(estimates
        .into_iter()
        .map(|e| {
            let minor = e.minor_strand_proportion();
            (e.unit, minor)
        })
        .collect())
}

/// check if a sequence looks like it is not
/// a telomeric repeat
fn check_telomeric_repeat(sequence: &str) -> bool {
    let repeat_period = check_repeats(sequence);
    repeat_period < REPEAT_PERIOD_THRESHOLD
}

/// Which way round a run reads, relative to its canonical unit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum Orientation {
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
pub(crate) fn orientation(sequence: &str, unit: &str) -> Orientation {
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
        let record = crate::fasta::Record::with_attrs("id1", None, genome.as_bytes());
        split_seq_by_distance(&record, DIST_FROM_CHROM_END, genome.len())
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
        chunk_fasta(left, CHUNK_LENGTH, None, "".into())
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
        chunk_fasta(left, CHUNK_LENGTH, None, "".into())
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
        calculate_indexes(chunks, CHUNK_LENGTH, None, "test".into(), 0).unwrap()
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
        let chunks = chunk_fasta(seq.as_bytes(), CHUNK_LENGTH, None, "".into());
        let positions = calculate_indexes(chunks, CHUNK_LENGTH, None, "test".into(), 0).unwrap();
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
        assert_eq!(consensus(&chunks, &[0, 0, 2, 2], b"TTAGGG"), "TTAGGG");
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

    #[test]
    fn test_tolerant_extends_both_ways() {
        // errors throughout, so the first seed is not at the start of the run
        let mut seq = "TTAGGG".repeat(100).into_bytes();
        for i in (0..300).step_by(18) {
            seq[i + 1] = b'C';
        }
        let seq = String::from_utf8(seq).unwrap();
        let runs = tolerant(&seq, 6);
        assert_eq!(runs.len(), 1);
        assert_eq!((runs[0].start, runs[0].end), (0, 600));
    }

    #[test]
    fn test_tolerant_is_symmetric() {
        // the same telomere read either way round gives the same run length
        let mut seq = format!("{}{}", random_seq(3000, 12), "TTAGGG".repeat(300)).into_bytes();
        for i in (3000..4800).step_by(23) {
            seq[i] = b'A';
        }
        let forward = String::from_utf8(seq).unwrap();
        let reverse = crate::utils::reverse_complement(&forward);
        let length = |s: &str| {
            tolerant(s, 6)
                .iter()
                .map(|r| r.end - r.start)
                .max()
                .unwrap()
        };
        assert!(length(&forward).abs_diff(length(&reverse)) <= 12);
    }

    /// Sequences to compare the fast run scanners with the reference: random
    /// sequence with noisy tandem arrays of several units, mixed case, Ns.
    fn scanner_cases() -> Vec<Vec<u8>> {
        let mut state = 7u64;
        let mut next = move |n: u64| {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) % n
        };
        let units = [
            "TTAGG", "AACCCT", "TTTAGGG", "A", "AT", "AAT", "ACGTTGCAAGTC", "TTAGGTTAGGTTAGGT",
            "ACGTACGTAAGGTTCCAGT", "CCTAA",
        ];
        let mut cases = Vec::new();
        for c in 0..120 {
            let mut seq: Vec<u8> = Vec::new();
            for _ in 0..(1 + next(6)) {
                // random stretch
                for _ in 0..next(300) {
                    seq.push(b"ACGT"[next(4) as usize]);
                }
                // an array with substitutions, indels, case changes and Ns
                let unit = units[next(units.len() as u64) as usize].as_bytes();
                let copies = 2 + next(60);
                let rate = [0, 2, 5, 10, 20][next(5) as usize];
                let lower = next(4) == 0;
                for _ in 0..copies {
                    for &b in unit {
                        match next(100) {
                            r if r < rate / 3 => seq.push(b"ACGT"[next(4) as usize]),
                            r if r < 2 * rate / 3 => {}
                            r if r < rate => {
                                seq.push(b);
                                seq.push(b"ACGTN"[next(5) as usize]);
                            }
                            _ => seq.push(if lower { b.to_ascii_lowercase() } else { b }),
                        }
                    }
                    if next(10) == 0 {
                        // switch case within an array
                        seq.extend(unit.iter().map(|b| b.to_ascii_lowercase()));
                    }
                }
            }
            if c % 7 == 0 {
                seq.extend_from_slice(b"NNNNNNNNNN");
            }
            cases.push(seq);
        }
        cases
    }

    fn runs_key(runs: &[RepeatPosition]) -> Vec<(usize, usize, String)> {
        runs.iter()
            .map(|r| (r.start, r.end, r.sequence.clone()))
            .collect()
    }

    #[test]
    fn test_scan_runs_matches_reference() {
        for (c, seq) in scanner_cases().iter().enumerate() {
            for k in 1..=20 {
                for tolerant in [true, false] {
                    let fast = find_runs(seq, k, "x", tolerant);
                    let reference = find_runs_reference(seq, k, "x", tolerant);
                    assert_eq!(
                        runs_key(&fast),
                        runs_key(&reference),
                        "case {c}, k {k}, tolerant {tolerant}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_filtered_runs_match_reference() {
        for (c, seq) in scanner_cases().iter().enumerate() {
            for k in [2, 3, 5, 6, 7, 12, 17] {
                for frequency in [0, 1, 3, 10] {
                    // tolerant: tolerant_runs, then the filter
                    let reference = tolerant_runs_reference(seq, k, "x");
                    let any = !reference.is_empty();
                    let filtered = RepeatPositions(reference).filter_by_frequency(frequency);
                    let (fast_any, fast) = filtered_runs(seq, k, "x", true, frequency);
                    assert_eq!(fast_any, any, "case {c}, k {k}");
                    assert_eq!(runs_key(&fast), runs_key(&filtered.0), "case {c}, k {k}");
                    // exact: chunk_fasta and calculate_indexes
                    let indexes = chunk_fasta(seq, k, None, "x".into());
                    let reference = calculate_indexes(indexes, k, None, "x".into(), frequency);
                    let (fast_any, fast) = filtered_runs(seq, k, "x", false, frequency);
                    assert_eq!(fast_any, reference.is_some(), "case {c}, k {k}");
                    assert_eq!(
                        runs_key(&fast),
                        runs_key(&reference.map(|r| r.0).unwrap_or_default()),
                        "case {c}, k {k}"
                    );
                }
            }
        }
    }

    #[test]
    fn test_smallest_period_matches_check_repeats() {
        for s in ["A", "AA", "AB", "ABA", "ABAB", "AACCT", "AAAAAAAA", "ATATATAT", "AATAATAAT", "ABCABCAB"] {
            assert_eq!(smallest_period(s.as_bytes()), check_repeats(s), "{s}");
        }
    }
}
