//! A tandem repeat finder for periods from 1 bp to a few kb, for
//! characterising chromosome ends (`tidk caps`).
//!
//! Copies of a tandem array share k-mers at a lag of one period. For each
//! position, the nearest previous occurrence of its k-mer gives a candidate
//! lag; lags that recur are candidate periods. For each candidate period,
//! positions whose k-mer recurs one period back (within a tolerance for
//! indels) are joined into arrays. Arrays found at different periods are
//! then chosen greedily by support, preferring shorter periods on ties.
//! Exact k-mer matches tolerate divergent copies: with k = 10, copies 85%
//! identical still share about 20% of k-mers.

use std::cell::RefCell;
use std::collections::{HashMap, HashSet};
use std::hash::BuildHasherDefault;

/// A set of k-mer codes, with a fast hash (nothing depends on its order).
pub(crate) type CodeSet = HashSet<u32, BuildHasherDefault<crate::anchor::KmerHasher>>;

/// k-mer length for finding lags.
const K: usize = 10;
/// k-mer length for comparing arrays (families, telomeric units).
pub(crate) const FAMILY_K: usize = 8;
/// A candidate period needs this many positions recurring at that lag.
const MIN_LAG_SUPPORT: u32 = 12;
/// At most this many candidate periods per sequence.
const MAX_CANDIDATES: usize = 40;
/// Arrays shorter than this are not reported.
const MIN_ARRAY_LENGTH: usize = 40;
/// Arrays need at least this many copies.
const MIN_COPIES: f64 = 1.9;
/// Above this period, arrays need [`MIN_COPIES_LONG`] copies: two copies of a
/// dispersed repeat a few kb apart otherwise pass as a long-period array.
pub(crate) const LONG_PERIOD: usize = 2000;
const MIN_COPIES_LONG: f64 = 3.0;
/// Above [`LONG_PERIOD`], matched positions in one array may be at most this
/// far apart (below it, up to one period).
const MAX_GAP_LONG: usize = 2000;
/// Above [`DENSITY_PERIOD`], matched positions must make up at least
/// [`MIN_DENSITY`] of the array after its first copy: otherwise chance k-mer
/// matches up to a period apart chain into an "array" across a whole window
/// of repeat-rich sequence. Copies 85% identical share about 20% of their
/// k-mers. In the end and interior windows of 40 assemblies, 944 of 966
/// arrays of 500-2,000 bp periods that were no more self-identical at their
/// period than at other lags were below this density, and 2 of 308 that were.
/// Below 500 bp the gaps allowed are short, chance arrays are few, and real
/// short arrays can be this sparse.
const DENSITY_PERIOD: usize = 500;
const MIN_DENSITY: f64 = 0.05;

/// The copies an array of period `d` needs.
fn min_copies(d: usize) -> f64 {
    if d > LONG_PERIOD {
        MIN_COPIES_LONG
    } else {
        MIN_COPIES
    }
}
/// The longest period an array can have in a sequence of `window` bp when
/// lags up to `max_period` are searched: the array must fit its copies in the
/// window, so above [`LONG_PERIOD`] this is a third of the window (6,666 bp
/// in 20 kb), whatever `max_period` is.
pub(crate) fn longest_period(max_period: usize, window: usize) -> usize {
    let long = (window as f64 / MIN_COPIES_LONG) as usize;
    let by_window = if long > LONG_PERIOD {
        long
    } else {
        ((window as f64 / MIN_COPIES) as usize).min(LONG_PERIOD)
    };
    max_period.min(by_window)
}
/// Previous occurrences followed per position when testing a lag
/// (low-complexity sequence can have many within one period).
const MAX_CHAIN: usize = 16;
/// Up to this period, the consensus is the commonest exact word.
const WORD_CONSENSUS_MAX_PERIOD: usize = 50;
/// An array whose period is a multiple of an overlapping array's is dropped
/// if the shorter period has at least this fraction of its support.
const MULTIPLE_SUPPORT: f64 = 0.7;
/// An array overlapping a better one by more than this fraction is dropped.
const MAX_OVERLAP: f64 = 0.5;

/// One tandem array, in sequence coordinates (0-based, half open).
#[derive(Debug, Clone)]
pub(crate) struct Array {
    pub(crate) start: usize,
    pub(crate) end: usize,
    pub(crate) period: usize,
    /// Positions recurring at one period back: the evidence for the array.
    pub(crate) support: usize,
    /// Majority consensus of the copies, ungapped, as read in the sequence.
    pub(crate) consensus: Vec<u8>,
    /// [`FAMILY_K`]-mers of the consensus in tandem, as read. The consensus
    /// averages out divergence between copies, so arrays of one family share
    /// most of these wherever they are; but indels blur an ungapped
    /// consensus of long units, hence also `repeated`.
    pub(crate) kmers: CodeSet,
    /// [`FAMILY_K`]-mers seen at least twice in the array, as read.
    pub(crate) repeated: CodeSet,
}

impl Array {
    pub(crate) fn len(&self) -> usize {
        self.end - self.start
    }
}

fn base_code(b: u8) -> Option<u32> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// 2-bit codes of the k-mer starting at each position (None if it has an N).
pub(crate) fn kmer_codes(seq: &[u8], k: usize) -> Vec<Option<u32>> {
    let mask = if k == 16 {
        u32::MAX
    } else {
        (1u32 << (2 * k)) - 1
    };
    let mut out = vec![None; seq.len().saturating_sub(k - 1)];
    let mut code = 0u32;
    let mut valid = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        match base_code(b) {
            Some(c) => {
                code = ((code << 2) | c) & mask;
                valid += 1;
            }
            None => valid = 0,
        }
        if i + 1 >= k && valid >= k {
            out[i + 1 - k] = Some(code);
        }
    }
    out
}

/// Reverse complement of a 2-bit k-mer code.
pub(crate) fn revcomp_code(mut code: u32, k: usize) -> u32 {
    let mut rc = 0u32;
    for _ in 0..k {
        rc = (rc << 2) | (3 - (code & 3));
        code >>= 2;
    }
    rc
}

/// The lesser of a k-mer code and its reverse complement.
pub(crate) fn canonical_code(code: u32, k: usize) -> u32 {
    code.min(revcomp_code(code, k))
}

/// Indel tolerance around a lag.
fn tolerance(period: usize) -> usize {
    period / 50
}

thread_local! {
    /// For each [`K`]-mer code, its last position in the sequence being
    /// scanned (`u32::MAX`: none yet); left empty after each scan.
    static LAST: RefCell<Vec<u32>> = const { RefCell::new(Vec::new()) };
}

/// Tandem arrays in a sequence, non-overlapping, sorted by start.
pub(crate) fn find_arrays(seq: &[u8], max_period: usize) -> Vec<Array> {
    let codes = kmer_codes(seq, K);
    if codes.is_empty() {
        return Vec::new();
    }

    // previous occurrence of each position's k-mer, and a histogram of the
    // lags to the nearest one
    let mut prev = vec![usize::MAX; codes.len()];
    // lags are shorter than the sequence, so the tables stop there
    let limit = max_period.min(codes.len());
    let mut hist = vec![0u32; limit + 1];
    let mut record = |i: usize, j: usize| {
        prev[i] = j;
        let d = i - j;
        if d <= max_period {
            hist[d] += 1;
        }
    };
    if codes.len() < u32::MAX as usize {
        LAST.with(|last| {
            let mut last = last.borrow_mut();
            if last.is_empty() {
                last.resize(1 << (2 * K), u32::MAX);
            }
            for (i, c) in codes.iter().enumerate() {
                if let Some(c) = c {
                    let slot = &mut last[*c as usize];
                    if *slot != u32::MAX {
                        record(i, *slot as usize);
                    }
                    *slot = i as u32;
                }
            }
            for c in codes.iter().flatten() {
                last[*c as usize] = u32::MAX;
            }
        });
    } else {
        let mut last: HashMap<u32, usize> = HashMap::new();
        for (i, c) in codes.iter().enumerate() {
            if let Some(c) = c {
                if let Some(j) = last.insert(*c, i) {
                    record(i, j);
                }
            }
        }
    }

    // candidate periods: best supported lags, with their indel tolerance
    let mut cumulative = vec![0u32; limit + 1];
    for d in 1..=limit {
        cumulative[d] = cumulative[d - 1] + hist[d];
    }
    // (no lag beyond `limit` has any support)
    let support = |d: usize| -> u32 {
        let t = tolerance(d);
        let (lo, hi) = (
            d.saturating_sub(t).max(1),
            (d + t).min(max_period).min(limit),
        );
        cumulative[hi] - cumulative[lo - 1]
    };
    let mut lags: Vec<(usize, u32)> = (1..=limit)
        .filter(|&d| hist[d] > 0)
        .map(|d| (d, support(d)))
        .collect();
    lags.sort_unstable_by_key(|&(d, s)| (std::cmp::Reverse(s), d));
    let mut candidates: Vec<usize> = Vec::new();
    for (d, s) in lags {
        if candidates.len() >= MAX_CANDIDATES || s < MIN_LAG_SUPPORT {
            break;
        }
        if candidates
            .iter()
            .all(|&c| c.abs_diff(d) > tolerance(c.max(d)))
        {
            candidates.push(d);
        }
    }

    let arrays: Vec<Array> = candidates
        .iter()
        .flat_map(|&d| arrays_at_period(seq, &prev, d))
        .collect();

    // a multiple of a shorter period over the same copies is the same array
    let overlap = |a: &Array, b: &Array| a.end.min(b.end).saturating_sub(a.start.max(b.start));
    let is_multiple = |long: &Array, short: &Array| {
        let k = (long.period as f64 / short.period as f64).round() as usize;
        k >= 2
            && long.period.abs_diff(k * short.period) <= tolerance(long.period).max(1)
            && overlap(long, short) as f64 >= MAX_OVERLAP * long.len().min(short.len()) as f64
            && short.support as f64 >= MULTIPLE_SUPPORT * long.support as f64
    };
    let keep: Vec<bool> = arrays
        .iter()
        .map(|a| {
            !arrays
                .iter()
                .any(|b| b.period < a.period && is_multiple(a, b))
        })
        .collect();
    let mut arrays: Vec<Array> = arrays
        .into_iter()
        .zip(keep)
        .filter_map(|(a, k)| k.then_some(a))
        .collect();

    // greedy choice: best supported first, shorter periods on ties
    arrays.sort_unstable_by_key(|a| (std::cmp::Reverse(a.support), a.period, a.start));
    let mut chosen: Vec<Array> = Vec::new();
    for a in arrays {
        let overlap: usize = chosen
            .iter()
            .map(|c| c.end.min(a.end).saturating_sub(c.start.max(a.start)))
            .sum();
        if (overlap as f64) <= MAX_OVERLAP * a.len() as f64 {
            chosen.push(a);
        }
    }
    chosen.sort_unstable_by_key(|a| a.start);
    chosen
}

/// Arrays with period `d`: runs of positions whose k-mer recurs `d` (± the
/// indel tolerance) back, joined across gaps of up to max(d, 2K). `prev`
/// links each position to the previous occurrence of its k-mer.
fn arrays_at_period(seq: &[u8], prev: &[usize], d: usize) -> Vec<Array> {
    let t = tolerance(d);
    let (lo, hi) = (d.saturating_sub(t).max(1), d + t);
    let matched = |i: usize| -> bool {
        let mut j = prev[i];
        for _ in 0..MAX_CHAIN {
            if j == usize::MAX || i - j > hi {
                return false;
            }
            if i - j >= lo {
                return true;
            }
            j = prev[j];
        }
        false
    };
    let max_gap = if d > LONG_PERIOD {
        MAX_GAP_LONG
    } else {
        d.max(2 * K)
    };
    let mut out = Vec::new();
    let mut run: Option<(usize, usize, usize)> = None; // first, last, count
    for i in lo..prev.len() {
        if !matched(i) {
            continue;
        }
        run = match run {
            Some((a, b, n)) if i - b <= max_gap => Some((a, i, n + 1)),
            other => {
                if let Some(r) = other {
                    out.extend(make_array(seq, d, r));
                }
                Some((i, i, 1))
            }
        };
    }
    if let Some(r) = run {
        out.extend(make_array(seq, d, r));
    }
    out
}

/// An array from a run of matched positions (first, last, count), if long
/// enough. The first copy has nothing one period back, so the array starts a
/// period before the first matched position.
fn make_array(seq: &[u8], d: usize, (first, last, n): (usize, usize, usize)) -> Option<Array> {
    let start = first.saturating_sub(d);
    let end = (last + K).min(seq.len());
    let len = end - start;
    if len < MIN_ARRAY_LENGTH || (len as f64) < min_copies(d) * d as f64 {
        return None;
    }
    // the copies after the first are the matched part
    if d > DENSITY_PERIOD && (n as f64) < MIN_DENSITY * (len - d) as f64 {
        return None;
    }
    let region = &seq[start..end];
    let consensus = consensus(region, d);
    let kmers = unit_kmers(&String::from_utf8_lossy(&consensus));
    Some(Array {
        start,
        end,
        period: d,
        support: n,
        consensus,
        kmers,
        repeated: repeated_kmers(region),
    })
}

/// The consensus unit of period `d`. For short periods, the commonest exact
/// d-mer, which indels don't disturb (e.g. plant telomeres mixing TTTAGGG
/// and TTTTAGGG); otherwise the column-wise majority of the copies, ungapped.
fn consensus(region: &[u8], d: usize) -> Vec<u8> {
    // uppercase words fit 2 bits a base; lowercase ones are counted apart,
    // so are left to the slices below
    if d <= WORD_CONSENSUS_MAX_PERIOD && !region.iter().any(|b| b"acgt".contains(b)) {
        if let Some(w) = commonest_word(region, d) {
            return w;
        }
    } else if d <= WORD_CONSENSUS_MAX_PERIOD {
        let mut counts: HashMap<&[u8], usize> = HashMap::new();
        for w in region.windows(d) {
            if w.iter().all(|&b| base_code(b).is_some()) {
                *counts.entry(w).or_default() += 1;
            }
        }
        if let Some((w, _)) = counts
            .into_iter()
            .max_by(|a, b| a.1.cmp(&b.1).then(b.0.cmp(a.0)))
        {
            return w.to_vec();
        }
    }
    (0..d)
        .map(|j| {
            let mut counts = [0usize; 4];
            for b in region.iter().skip(j).step_by(d) {
                if let Some(c) = base_code(*b) {
                    counts[c as usize] += 1;
                }
            }
            let best = (0..4).max_by_key(|&c| counts[c]).unwrap();
            if counts[best] == 0 {
                b'N'
            } else {
                b"ACGT"[best]
            }
        })
        .collect()
}

/// The commonest d-mer of A, C, G and T in a region (d <= 64), the
/// lexicographically first on ties, as `consensus` counts words.
fn commonest_word(region: &[u8], d: usize) -> Option<Vec<u8>> {
    let mask = if d >= 64 {
        u128::MAX
    } else {
        (1u128 << (2 * d)) - 1
    };
    let mut words: Vec<u128> = Vec::with_capacity(region.len());
    let (mut code, mut valid) = (0u128, 0usize);
    for &b in region {
        match base_code(b) {
            Some(c) => {
                code = ((code << 2) | c as u128) & mask;
                valid += 1;
            }
            None => valid = 0,
        }
        if valid >= d {
            words.push(code);
        }
    }
    words.sort_unstable();
    let (mut best, mut best_count) = (None, 0usize);
    for run in words.chunk_by(|a, b| a == b) {
        if run.len() > best_count {
            best = Some(run[0]);
            best_count = run.len();
        }
    }
    best.map(|w| {
        (0..d)
            .map(|j| b"ACGT"[((w >> (2 * (d - 1 - j))) & 3) as usize])
            .collect()
    })
}

/// [`FAMILY_K`]-mers occurring at least twice in a region.
fn repeated_kmers(region: &[u8]) -> CodeSet {
    const _: () = assert!(2 * FAMILY_K <= 32);
    let mut seen = vec![0u64; (1usize << (2 * FAMILY_K)).div_ceil(64)];
    let mut repeated = CodeSet::default();
    for c in kmer_codes(region, FAMILY_K).into_iter().flatten() {
        let (word, bit) = ((c >> 6) as usize, c & 63);
        if seen[word] >> bit & 1 == 1 {
            repeated.insert(c);
        } else {
            seen[word] |= 1 << bit;
        }
    }
    repeated
}

/// [`FAMILY_K`]-mers of a unit repeated in tandem, e.g. a telomeric unit.
pub(crate) fn unit_kmers(unit: &str) -> CodeSet {
    let copies = FAMILY_K.div_ceil(unit.len().max(1)) + 2;
    kmer_codes(unit.repeat(copies).as_bytes(), FAMILY_K)
        .into_iter()
        .flatten()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_seq(len: usize, mut state: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(state >> 62) as usize]
            })
            .collect()
    }

    /// `copies` copies of `unit`, each base substituted with probability `div`.
    fn diverged(unit: &[u8], copies: usize, div: f64, mut state: u64) -> Vec<u8> {
        let mut out = Vec::new();
        for _ in 0..copies {
            for &b in unit {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let r = (state >> 11) as f64 / (1u64 << 53) as f64;
                out.push(if r < div {
                    b"ACGT"[((state >> 3) % 4) as usize]
                } else {
                    b
                });
            }
        }
        out
    }

    #[test]
    fn test_kmer_codes_and_revcomp() {
        let codes = kmer_codes(b"ACGTNACGTA", 4);
        assert_eq!(codes[0], Some(0b00011011));
        assert_eq!(codes[1], None);
        assert_eq!(codes[5], Some(0b00011011));
        // ACGT is its own reverse complement
        assert_eq!(revcomp_code(0b00011011, 4), 0b00011011);
        // AAAA -> TTTT
        assert_eq!(revcomp_code(0, 4), 0b11111111);
    }

    #[test]
    fn test_finds_satellite_in_random_sequence() {
        let unit = random_seq(176, 1);
        let mut seq = random_seq(3000, 2);
        seq.extend(diverged(&unit, 20, 0.08, 3));
        seq.extend(random_seq(3000, 4));
        let arrays = find_arrays(&seq, 2000);
        let a = arrays.iter().max_by_key(|a| a.len()).unwrap();
        assert!(a.period.abs_diff(176) <= 3, "period {}", a.period);
        assert!(a.start.abs_diff(3000) < 200 && a.end.abs_diff(3000 + 176 * 20) < 200);
    }

    #[test]
    fn test_long_unit_needs_three_copies() {
        // a 2,150 bp unit in 8 diverged copies is an array at its period
        let unit = random_seq(2150, 11);
        let mut seq = random_seq(1000, 12);
        seq.extend(diverged(&unit, 8, 0.05, 13));
        seq.extend(random_seq(1000, 14));
        let arrays = find_arrays(&seq, 10_000);
        let a = arrays.iter().max_by_key(|a| a.len()).unwrap();
        assert!(a.period.abs_diff(2150) <= 43, "period {}", a.period);
        // two copies of a 3 kb dispersed repeat 6 kb apart are not
        let element = random_seq(3000, 15);
        let mut seq = random_seq(2000, 16);
        seq.extend(diverged(&element, 1, 0.02, 17));
        seq.extend(random_seq(3000, 18));
        seq.extend(diverged(&element, 1, 0.02, 19));
        seq.extend(random_seq(2000, 20));
        let arrays = find_arrays(&seq, 10_000);
        assert!(
            arrays.iter().all(|a| a.period < 5000),
            "{:?}",
            arrays.iter().map(|a| a.period).collect::<Vec<_>>()
        );
    }

    /// Random sequence with `at` of its bases A or T.
    fn at_rich(len: usize, at: f64, mut state: u64) -> Vec<u8> {
        (0..len)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let r = (state >> 11) as f64 / (1u64 << 53) as f64;
                if r < at / 2.0 {
                    b'A'
                } else if r < at {
                    b'T'
                } else if r < at + (1.0 - at) / 2.0 {
                    b'C'
                } else {
                    b'G'
                }
            })
            .collect()
    }

    #[test]
    fn test_no_chance_arrays_of_1_to_2_kb_periods() {
        // AT-rich sequence has chance 10-mer matches at any lag; up to a
        // period apart they must not chain into arrays of 0.5-2 kb periods,
        // whichever lags are searched
        for seed in 40..46 {
            let seq = at_rich(20_000, 0.8, seed);
            for max_period in [2000, 6666, 10_000] {
                let arrays = find_arrays(&seq, max_period);
                assert!(
                    arrays.iter().all(|a| a.period <= DENSITY_PERIOD),
                    "seed {seed}, max period {max_period}: {:?}",
                    arrays
                        .iter()
                        .map(|a| (a.period, a.len()))
                        .collect::<Vec<_>>()
                );
            }
        }
    }

    #[test]
    fn test_sparse_array_needs_density() {
        // a 1,200 bp unit in 6 copies 85% identical is an array (about 20%
        // of its 10-mers recur a period back)
        let unit = random_seq(1200, 50);
        let mut seq = at_rich(4000, 0.8, 51);
        seq.extend(diverged(&unit, 6, 0.15, 52));
        seq.extend(at_rich(4000, 0.8, 53));
        for max_period in [2000, 6666, 10_000] {
            let arrays = find_arrays(&seq, max_period);
            let long: Vec<_> = arrays
                .iter()
                .filter(|a| a.period > DENSITY_PERIOD)
                .collect();
            assert_eq!(long.len(), 1, "max period {max_period}");
            assert!(long[0].period.abs_diff(1200) <= 24);
            assert!(long[0].start.abs_diff(4000) < 300 && long[0].end.abs_diff(11_200) < 300);
        }
    }

    #[test]
    fn test_longest_period_is_a_third_of_the_window() {
        assert_eq!(longest_period(10_000, 20_000), 6666);
        assert_eq!(longest_period(10_000, 30_000), 10_000);
        assert_eq!(longest_period(10_000, 100_000), 10_000);
        assert_eq!(longest_period(2000, 20_000), 2000);
        // short windows: two copies of a unit up to 2 kb
        assert_eq!(longest_period(10_000, 5000), 2000);
        assert_eq!(longest_period(10_000, 1900), 1000);
        // a 20 kb window of one unit: found up to the limit, not above it
        for (period, found) in [(6000, true), (6600, true), (7000, false), (9000, false)] {
            let seq = diverged(&random_seq(period, 31), 20_000 / period + 1, 0.02, 32);
            let arrays = find_arrays(&seq[..20_000], 10_000);
            let long = arrays
                .iter()
                .any(|a| a.period.abs_diff(period) <= period / 50);
            assert_eq!(long, found, "period {period}");
            assert!(arrays
                .iter()
                .all(|a| a.period <= longest_period(10_000, 20_000)));
        }
    }

    #[test]
    fn test_no_long_arrays_in_at_rich_random_sequence() {
        // 70% AT, as in fly genomes: chance k-mer matches must not chain into
        // a long-period array
        let mut state = 21u64;
        let seq: Vec<u8> = (0..20_000)
            .map(|_| {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let r = (state >> 11) as f64 / (1u64 << 53) as f64;
                if r < 0.35 {
                    b'A'
                } else if r < 0.7 {
                    b'T'
                } else if r < 0.85 {
                    b'C'
                } else {
                    b'G'
                }
            })
            .collect();
        let arrays = find_arrays(&seq, 10_000);
        assert!(
            arrays.iter().all(|a| a.period <= LONG_PERIOD),
            "{:?}",
            arrays.iter().map(|a| a.period).collect::<Vec<_>>()
        );
    }

    #[test]
    fn test_finds_telomere_and_microsatellite() {
        let mut seq = b"CCCTAA".repeat(300);
        seq.extend(random_seq(2000, 5));
        seq.extend(b"AT".repeat(50));
        seq.extend(random_seq(2000, 6));
        let arrays = find_arrays(&seq, 2000);
        assert_eq!(arrays[0].period, 6);
        assert!(arrays[0].len() >= 1790);
        assert!(arrays.iter().any(|a| a.period == 2 && a.len() >= 90));
    }

    #[test]
    fn test_no_arrays_in_random_sequence() {
        let arrays = find_arrays(&random_seq(20_000, 7), 2000);
        let bp: usize = arrays.iter().map(|a| a.len()).sum();
        assert!(bp < 500, "{bp} bp of arrays in random sequence");
    }

    #[test]
    fn test_long_period() {
        let unit = random_seq(1200, 8);
        let mut seq = random_seq(1000, 9);
        seq.extend(diverged(&unit, 5, 0.05, 10));
        let arrays = find_arrays(&seq, 2000);
        assert!(arrays
            .iter()
            .any(|a| a.period.abs_diff(1200) <= 24 && a.len() > 4000));
    }

    #[test]
    fn test_family_kmers_shared_between_copies() {
        let unit = random_seq(150, 11);
        let a = find_arrays(&diverged(&unit, 10, 0.05, 12), 2000);
        let b = find_arrays(&diverged(&unit, 10, 0.05, 13), 2000);
        let (ka, kb) = (&a[0].repeated, &b[0].repeated);
        let shared = ka.intersection(kb).count() as f64;
        assert!(shared / (ka.len() + kb.len()) as f64 * 2.0 > 0.5);
    }
}
