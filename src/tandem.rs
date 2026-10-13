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

use std::collections::{HashMap, HashSet};

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
    pub(crate) kmers: HashSet<u32>,
    /// [`FAMILY_K`]-mers seen at least twice in the array, as read.
    pub(crate) repeated: HashSet<u32>,
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

/// Tandem arrays in a sequence, non-overlapping, sorted by start.
pub(crate) fn find_arrays(seq: &[u8], max_period: usize) -> Vec<Array> {
    let codes = kmer_codes(seq, K);
    if codes.is_empty() {
        return Vec::new();
    }

    // previous occurrence of each position's k-mer, and a histogram of the
    // lags to the nearest one
    let mut last: HashMap<u32, usize> = HashMap::new();
    let mut prev = vec![usize::MAX; codes.len()];
    let mut hist = vec![0u32; max_period + 1];
    for (i, c) in codes.iter().enumerate() {
        if let Some(c) = c {
            if let Some(j) = last.insert(*c, i) {
                prev[i] = j;
                let d = i - j;
                if d <= max_period {
                    hist[d] += 1;
                }
            }
        }
    }

    // candidate periods: best supported lags, with their indel tolerance
    let support = |d: usize| -> u32 {
        let t = tolerance(d);
        hist[d.saturating_sub(t).max(1)..=(d + t).min(max_period)]
            .iter()
            .sum()
    };
    let mut lags: Vec<usize> = (1..=max_period).filter(|&d| hist[d] > 0).collect();
    lags.sort_unstable_by_key(|&d| (std::cmp::Reverse(support(d)), d));
    let mut candidates: Vec<usize> = Vec::new();
    for d in lags {
        if candidates.len() >= MAX_CANDIDATES || support(d) < MIN_LAG_SUPPORT {
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
    let max_gap = d.max(2 * K);
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
    if len < MIN_ARRAY_LENGTH || (len as f64) < MIN_COPIES * d as f64 {
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
    if d <= WORD_CONSENSUS_MAX_PERIOD {
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

/// [`FAMILY_K`]-mers occurring at least twice in a region.
fn repeated_kmers(region: &[u8]) -> HashSet<u32> {
    let mut seen = HashSet::new();
    let mut repeated = HashSet::new();
    for c in kmer_codes(region, FAMILY_K).into_iter().flatten() {
        if !seen.insert(c) {
            repeated.insert(c);
        }
    }
    repeated
}

/// [`FAMILY_K`]-mers of a unit repeated in tandem, e.g. a telomeric unit.
pub(crate) fn unit_kmers(unit: &str) -> HashSet<u32> {
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
