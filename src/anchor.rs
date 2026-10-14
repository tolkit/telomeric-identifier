//! Anchoring long reads to chromosome ends, for `tidk caps --reads`.
//!
//! The terminal `bait` bp of every end, read from the tip inwards, are
//! indexed by their k-mers, keeping only k-mers that occur once across all
//! baits. Each read is scanned at every [`STRIDE`]th position; hits to one
//! end and strand that agree on a diagonal (read position minus bait
//! position, within [`DIAGONAL_TOLERANCE`]) anchor the read to that end. The
//! read position of the tip is extrapolated from the anchor nearest the tip.
//! Every anchored read is returned with its span along the bait; a read
//! covering the tip (starting no more than [`TIP`] bp inside it) overhangs
//! it by the bases beyond: the chromosome continues past the assembly end if
//! many reads do. No alignment is done, so the tip is only
//! placed to within the read's indels between it and the nearest anchor.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// k-mer length for anchoring.
const K: usize = 25;
/// Read positions looked up.
const STRIDE: usize = 8;
/// Anchors on one diagonal agree to within this many bp.
const DIAGONAL_TOLERANCE: i64 = 200;
/// A read needs this many anchors on its diagonal...
const MIN_ANCHORS: usize = 10;
/// ...spanning at least this many bp of the bait.
const MIN_SPAN: usize = 2000;
/// A read starting at most this far inside the tip still covers it.
pub(crate) const TIP: i64 = 200;

/// A fast hasher for k-mer codes (multiplicative, as in FxHash).
#[derive(Default)]
pub(crate) struct KmerHasher(u64);

impl Hasher for KmerHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
        }
    }
    fn write_u64(&mut self, n: u64) {
        self.0 = (self.0.rotate_left(5) ^ n).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}

type KmerMap<V> = HashMap<u64, V, BuildHasherDefault<KmerHasher>>;

/// Where a k-mer is in the baits.
#[derive(Clone, Copy)]
struct Site {
    end: u32,
    pos: u32,
    /// The k-mer reads forward (not reverse complemented) in the bait.
    forward: bool,
}

/// The k-mer index of all end baits.
pub(crate) struct Baits {
    index: KmerMap<Option<Site>>,
    /// Per end, the distance from the tip to its first unique k-mer.
    pub(crate) first_anchor: Vec<Option<usize>>,
}

/// A read anchored to an end.
pub(crate) struct Hit {
    pub(crate) end: usize,
    /// Bases of the read beyond the tip (negative: the read starts inside).
    pub(crate) overhang: i64,
    /// Bait position of the read's inner end (from the tip).
    pub(crate) span_end: i64,
    /// The overhang, read in the bait's orientation (outside towards the tip),
    /// if at least `keep` bp.
    pub(crate) sequence: Option<Vec<u8>>,
}

fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Canonical k-mers of a sequence: (position, canonical code, reads forward).
fn kmers(seq: &[u8], stride: usize) -> impl Iterator<Item = (usize, u64, bool)> + '_ {
    let mask = (1u64 << (2 * K)) - 1;
    let shift = 2 * (K - 1);
    let (mut fwd, mut rev, mut valid) = (0u64, 0u64, 0usize);
    seq.iter().enumerate().filter_map(move |(i, &b)| {
        match code(b) {
            Some(c) => {
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << shift);
                valid += 1;
            }
            None => valid = 0,
        }
        if valid < K {
            return None;
        }
        let start = i + 1 - K;
        if !start.is_multiple_of(stride) {
            None
        } else if fwd <= rev {
            Some((start, fwd, true))
        } else {
            Some((start, rev, false))
        }
    })
}

pub(crate) fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

impl Baits {
    /// Index baits, each read from the tip inwards.
    pub(crate) fn new(baits: &[&[u8]]) -> Self {
        let mut index: KmerMap<Option<Site>> = KmerMap::default();
        for (end, bait) in baits.iter().enumerate() {
            for (pos, c, forward) in kmers(bait, 1) {
                index
                    .entry(c)
                    .and_modify(|s| *s = None)
                    .or_insert(Some(Site {
                        end: end as u32,
                        pos: pos as u32,
                        forward,
                    }));
            }
        }
        index.retain(|_, s| s.is_some());
        let mut first_anchor = vec![None; baits.len()];
        for s in index.values().flatten() {
            let f = &mut first_anchor[s.end as usize];
            if f.is_none_or(|p| (s.pos as usize) < p) {
                *f = Some(s.pos as usize);
            }
        }
        Baits {
            index,
            first_anchor,
        }
    }

    /// Anchor a read to an end. Overhangs of at least `keep` bp are returned
    /// with their sequence.
    pub(crate) fn anchor(&self, read: &[u8], keep: usize) -> Option<Hit> {
        // hits per (end, same strand): (diagonal, bait position)
        let mut hits: HashMap<(u32, bool), Vec<(i64, u32)>> = HashMap::new();
        let len = read.len() as i64;
        for (i, c, forward) in kmers(read, STRIDE) {
            if let Some(Some(s)) = self.index.get(&c) {
                let same = forward == s.forward;
                // read coordinate along the bait's orientation
                let r = if same {
                    i as i64
                } else {
                    len - i as i64 - K as i64
                };
                hits.entry((s.end, same))
                    .or_default()
                    .push((r - s.pos as i64, s.pos));
            }
        }
        let ((end, same), mut h) = hits.into_iter().max_by_key(|(k, v)| (v.len(), *k))?;
        if h.len() < MIN_ANCHORS {
            return None;
        }
        h.sort_unstable();
        let median = h[h.len() / 2].0;
        h.retain(|&(d, _)| (d - median).abs() <= DIAGONAL_TOLERANCE);
        let (lo, hi) = h
            .iter()
            .fold((u32::MAX, 0), |(lo, hi), &(_, p)| (lo.min(p), hi.max(p)));
        if h.len() < MIN_ANCHORS || ((hi - lo) as usize) < MIN_SPAN {
            return None;
        }
        // the diagonal at the anchor nearest the tip
        let overhang = h.iter().min_by_key(|&&(_, p)| p).unwrap().0;
        let sequence = (overhang >= keep as i64).then(|| {
            let n = overhang as usize;
            if same {
                read[..n].to_vec()
            } else {
                revcomp(&read[read.len() - n..]).to_ascii_uppercase()
            }
        });
        Some(Hit {
            end: end as usize,
            overhang,
            span_end: len - overhang,
            sequence: sequence.map(|s| s.to_ascii_uppercase()),
        })
    }
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

    #[test]
    fn test_overhanging_read_both_strands() {
        let bait0 = random_seq(30_000, 1);
        let bait1 = random_seq(30_000, 2);
        let baits = Baits::new(&[&bait0, &bait1]);
        // 3 kb beyond the tip of end 1, then 8 kb of the bait
        let beyond = random_seq(3000, 3);
        let read = [beyond.clone(), bait1[..8000].to_vec()].concat();
        let hit = baits.anchor(&read, 1000).unwrap();
        assert_eq!((hit.end, hit.overhang), (1, 3000));
        assert_eq!(hit.sequence.unwrap(), beyond);
        // the same read reverse complemented
        let hit = baits.anchor(&revcomp(&read), 1000).unwrap();
        assert_eq!((hit.end, hit.overhang), (1, 3000));
        assert_eq!(hit.sequence.unwrap(), beyond);
    }

    #[test]
    fn test_read_inside_bait_is_not_a_tip_read() {
        let bait = random_seq(30_000, 4);
        let baits = Baits::new(&[&bait]);
        let hit = baits.anchor(&bait[5000..15_000], 1000).unwrap();
        assert_eq!((hit.overhang, hit.span_end), (-5000, 15_000));
        assert!(hit.overhang < -TIP);
        // starting just inside the tip still covers it
        let hit = baits.anchor(&bait[100..10_000], 1000).unwrap();
        assert_eq!(hit.overhang, -100);
        assert!(hit.sequence.is_none());
    }

    #[test]
    fn test_kmers_shared_between_ends_are_not_anchors() {
        let shared = random_seq(5000, 5);
        let a = [shared.clone(), random_seq(25_000, 6)].concat();
        let b = [shared, random_seq(25_000, 7)].concat();
        let baits = Baits::new(&[&a, &b]);
        assert!((4975..=5000).contains(&baits.first_anchor[0].unwrap()));
        // a read beyond the tip of end 0 is still placed, from anchors past the shared part
        let read = [random_seq(2000, 8), a[..10_000].to_vec()].concat();
        let hit = baits.anchor(&read, 1000).unwrap();
        assert_eq!((hit.end, hit.overhang), (0, 2000));
    }
}
