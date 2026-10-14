//! `tidk caps`: what is at the chromosome ends, with or without a canonical
//! telomere.
//!
//! For each sequence, a window at both ends (read from the tip inwards, so
//! all ends are comparable) and a few interior windows as a baseline:
//!
//! 1. the telomere call at each end, as `tidk ends` makes it, if a unit is
//!    given or can be discovered;
//! 2. tandem arrays in every window ([`crate::tandem`]), grouped into
//!    families across the whole assembly by shared k-mers (rotation and
//!    strand free), and the family covering most of each window;
//! 3. end enrichment: the family that is top at most ends, how often it is
//!    top in interior windows, and whether its arrays read the same way
//!    round from the tip at every end.
//!
//! An assembly's ends are classed `telomeric` or `terminal_satellite` when
//! one family is top at >= half the ends, covers >= 20% of those windows
//! (median) and is >= 2x as frequent at ends as in the interior; otherwise
//! `none`.
//!
//! Ends are flagged when an assembly gap lies near the tip, and (with
//! `--reads`) when read depth just inside the tip is low. Two more
//! end-structure measures: whether the two ends of each sequence share their
//! top family (against random pairing), and segments copied from the tip
//! between ends (non-tandem sequence collinear from both tips, as left by
//! break-induced replication or ALT-like end copying).

use crate::anchor::{self, Baits};
use crate::ends::{self, EndsConfig, Side, Status};
use crate::tandem::{self, Array, FAMILY_K};
use crate::{explore, open_fasta_reader, utils};
use anyhow::{bail, Context, Result};
use clap::crate_version;
use rayon::prelude::*;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::{create_dir_all, File};
use std::hash::BuildHasherDefault;
use std::io::{LineWriter, Write};
use std::path::PathBuf;

/// Arrays are in the same family if their k-mer sets have at least this
/// Jaccard similarity.
const FAMILY_JACCARD: f64 = 0.3;
/// A family is telomeric if it holds at least this fraction of the unit's k-mers.
const TELOMERIC_SHARE: f64 = 0.5;
/// An array is telomeric (family 0) if at least this fraction of its k-mers
/// are the unit's, so variant-rich telomeres stay with the telomere.
const TELOMERIC_COVER: f64 = 0.3;
/// End-enrichment thresholds for a cap (see the module docs).
const CAP_END_FRACTION: f64 = 0.5;
const CAP_MEDIAN_COVER: f64 = 0.2;
const CAP_ENRICHMENT: f64 = 2.0;
/// Shared non-tandem sequence: k-mer length, bp shared for a partner, and
/// partners for a window to count as sharing. As for caps, at least half the
/// ends must share, at least [`CAP_ENRICHMENT`]x the interior fraction.
const SHARED_K: usize = 16;
const MIN_SHARED: usize = 2000;
const SHARED_MIN_PARTNERS: usize = 2;
const SHARED_END_FRACTION: f64 = 0.5;
/// Tip copies: shared k-mers collinear from both tips (same orientation read
/// from the tip), in non-tandem coordinates (tandem arrays, which differ in
/// length between copies, removed): starting within TIP_COPY_START bp of
/// both, gaps up to TIP_COPY_GAP bp; a partner needs TIP_COPY_MIN bp of k-mer
/// cover at >= TIP_COPY_DENSITY of the chain's extent (~ 97% identity).
const TIP_COPY_START: usize = 5000;
const TIP_COPY_GAP: usize = 5000;
const TIP_COPY_MIN: usize = 5000;
const TIP_COPY_DENSITY: f64 = 0.5;
/// Telomere calls: minimum length and offset from the tip, as in `tidk ends`.
const TELOMERE_MIN_LENGTH: usize = 200;
const TELOMERE_MAX_OFFSET: usize = 1000;
/// Reads overhanging the tip by at least this many bp continue the chromosome.
const MIN_OVERHANG: usize = 1000;
/// An end needs this many tip reads for a call.
const MIN_TIP_READS: usize = 3;
/// Read depth is measured over [DEPTH_AT, DEPTH_AT + DEPTH_SPAN) bp inside
/// the tip; an end below LOW_DEPTH x the median over ends is flagged.
const DEPTH_AT: i64 = 2000;
const DEPTH_SPAN: i64 = 1000;
const LOW_DEPTH: f64 = 0.25;
/// A run of at least this many Ns is an assembly gap; an end window with one
/// within GAP_NEAR_TIP bp of the tip is flagged (its tip is a separate contig).
const MIN_GAP_RUN: usize = 10;
const GAP_NEAR_TIP: usize = 5000;
/// Overhangs analysed for tandem repeats per end.
const OVERHANGS_PER_END: usize = 50;
/// An overhang's top family must cover this fraction of it.
const OVERHANG_MIN_COVER: f64 = 0.2;
// unit discovery at the ends, as in `tidk ends`
const DISCOVERY_MINIMUM: usize = 5;
const DISCOVERY_MAXIMUM: usize = 12;
const DISCOVERY_THRESHOLD: usize = 20;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Kind {
    Start,
    End,
    Interior,
}

impl Kind {
    fn name(&self) -> &'static str {
        match self {
            Kind::Start => "start",
            Kind::End => "end",
            Kind::Interior => "interior",
        }
    }

    fn is_end(&self) -> bool {
        *self != Kind::Interior
    }
}

/// One window. End windows are read from the tip inwards.
struct Window {
    id: String,
    sequence_length: usize,
    kind: Kind,
    /// Position of the window in the sequence.
    position: usize,
    seq: Vec<u8>,
    telomere: Option<ends::EndCall>,
}

/// Results for one window.
struct WindowResult {
    arrays: Vec<Array>,
    families: Vec<usize>,
    top: Option<TopFamily>,
    tandem_cover: f64,
    /// Non-tandem sequence shared with windows of the same kind elsewhere.
    shared: Shared,
}

/// Sequence a window shares with windows of the same kind (ends with ends,
/// interior with interior) on other sequences, tandem arrays masked.
#[derive(Default, Clone)]
struct Shared {
    /// Windows sharing at least [`MIN_SHARED`] bp.
    partners: usize,
    /// Shared bp with those partners, in the same and opposite orientation.
    same: usize,
    opposite: usize,
    /// End windows: partners sharing a segment copied from the tip, and the
    /// longest such segment (bp from the tip).
    tip_partners: usize,
    tip_extent: usize,
}

struct TopFamily {
    family: usize,
    cover: f64,
    period: usize,
    /// The family's arrays read the same way round as its representative.
    forward: bool,
}

struct Family {
    period: usize,
    consensus: String,
    /// Canonical k-mers of the representative's consensus.
    kmers: HashSet<u32>,
    /// Canonical k-mers repeated in the representative.
    repeated: HashSet<u32>,
    /// k-mers of the representative array as read, for orientation.
    forward: HashSet<u32>,
    telomeric: bool,
}

/// The entry point for `tidk caps`.
pub fn caps(matches: &clap::ArgMatches) -> Result<()> {
    let input_fasta = matches
        .get_one::<PathBuf>("fasta")
        .expect("errored by clap");
    let string = matches.get_one::<String>("string");
    let window = *matches
        .get_one::<usize>("window")
        .expect("defaulted by clap");
    let n_interior = *matches
        .get_one::<usize>("interior")
        .expect("defaulted by clap");
    let max_period = *matches
        .get_one::<usize>("max-period")
        .expect("defaulted by clap");
    let min_sequence_length = *matches
        .get_one::<usize>("min-sequence-length")
        .expect("defaulted by clap");
    let output = matches
        .get_one::<PathBuf>("output")
        .expect("errored by clap");
    let outdir = matches
        .get_one::<PathBuf>("dir")
        .expect("defaulted by clap");
    let pattern = |name: &str| -> Result<Option<Regex>> {
        matches
            .get_one::<String>(name)
            .map(|p| Regex::new(p).with_context(|| format!("invalid --{name} regular expression")))
            .transpose()
    };
    let include = pattern("include")?;
    let exclude = pattern("exclude")?;

    let reads: Vec<PathBuf> = matches
        .get_many::<PathBuf>("reads")
        .map(|r| r.cloned().collect())
        .unwrap_or_default();
    let bait = *matches.get_one::<usize>("bait").expect("defaulted by clap");
    let threads = *matches
        .get_one::<usize>("threads")
        .expect("defaulted by clap");
    if threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .context("could not set up the thread pool")?;
    }

    // windows; telomeres are called later, once the unit is known. With
    // reads, also the terminal `bait` bp of each end, from the tip inwards,
    // in the order of the end windows.
    let mut windows = Vec::new();
    let mut baits: Vec<Vec<u8>> = Vec::new();
    for result in open_fasta_reader(input_fasta)?.records() {
        let record = result?;
        let id = record.id().to_owned();
        let len = record.seq().len();
        if len < min_sequence_length
            || len < 2 * window
            || include.as_ref().is_some_and(|r| !r.is_match(&id))
            || exclude.as_ref().is_some_and(|r| r.is_match(&id))
        {
            continue;
        }
        let seq = record.seq().to_ascii_uppercase();
        let mut push = |kind, position, s: Vec<u8>| {
            windows.push(Window {
                id: id.clone(),
                sequence_length: len,
                kind,
                position,
                seq: s,
                telomere: None,
            })
        };
        push(Kind::Start, 0, seq[..window].to_vec());
        push(Kind::End, len - window, seq[len - window..].to_vec());
        if !reads.is_empty() {
            let b = bait.min(len / 2);
            baits.push(seq[..b].to_vec());
            baits.push(anchor::revcomp(&seq[len - b..]));
        }
        if len >= (n_interior + 2) * window {
            for i in 1..=n_interior {
                let mid = len * i / (n_interior + 1);
                push(
                    Kind::Interior,
                    mid - window / 2,
                    seq[mid - window / 2..mid + window / 2].to_vec(),
                );
            }
        }
    }
    let n_ends = windows.iter().filter(|w| w.kind.is_end()).count();
    if n_ends == 0 {
        bail!(
            "No sequences of at least {}bp in the input (after --include/--exclude).",
            min_sequence_length.max(2 * window)
        );
    }

    // telomere calls on the ends as they are in the sequence
    let (unit, unit_source) = match string {
        Some(s) => (Some(ends::parse_unit(s)?), "given"),
        None => (discover_unit(&windows)?, "discovered"),
    };
    if let Some(u) = &unit {
        eprintln!("[+]\tTelomeric repeat unit ({unit_source}): {u}");
        let config = EndsConfig {
            unit: u.clone(),
            min_length: TELOMERE_MIN_LENGTH,
            max_offset: TELOMERE_MAX_OFFSET,
            error_tolerant: true,
        };
        windows.par_iter_mut().for_each(|w| {
            w.telomere = match w.kind {
                Kind::Start => Some(ends::call_end(&w.seq, Side::Start, w.position, &config)),
                Kind::End => Some(ends::call_end(&w.seq, Side::End, w.position, &config)),
                Kind::Interior => None,
            };
        });
    } else {
        eprintln!("[+]\tNo telomeric repeat unit given or found at the ends; no telomere calls");
    }
    // from here on, end windows read from the tip inwards
    windows.par_iter_mut().for_each(|w| {
        if w.kind == Kind::End {
            w.seq = revcomp(&w.seq);
        }
    });
    eprintln!(
        "[+]\tFinding tandem repeats in {} end and {} interior windows of {window}bp",
        n_ends,
        windows.len() - n_ends
    );

    let (results, mut builder, summary) = analyse(&windows, unit.as_deref(), max_period);

    let gap_ends = windows
        .iter()
        .filter(|w| w.kind.is_end() && gaps(&w.seq).1.is_some_and(|g| g < GAP_NEAR_TIP))
        .count();
    if gap_ends > 0 {
        eprintln!("[+]\t{gap_ends} ends have an assembly gap within {GAP_NEAR_TIP}bp of the tip");
    }

    create_dir_all(outdir)?;
    let prefix = format!("{}/{}", outdir.display(), output.display());
    let reads_summary = if reads.is_empty() {
        serde_json::Value::Null
    } else {
        let ends: Vec<&Window> = windows.iter().filter(|w| w.kind.is_end()).collect();
        let per_end = end_reads(&reads, &baits, &ends, &mut builder, max_period, &prefix)?;
        let s = summarise_reads(&per_end, &builder.families);
        eprintln!(
            "[+]	Reads: {} ends open, {} mixed, {} closed, {} without reads. Commonest overhang family: {} ({} bp period), top at {} open ends",
            s["open"], s["mixed"], s["closed"], s["no_reads"],
            s["overhang_family_consensus_short"], s["overhang_family_period"], s["overhang_family_open_ends"]
        );
        s
    };
    let families = &builder.families;
    write_windows(&windows, &results, families, &format!("{prefix}.caps.tsv"))?;
    write_families(
        &windows,
        &results,
        families,
        &format!("{prefix}.caps.families.tsv"),
    )?;
    let json = serde_json::json!({
        "tidk_version": crate_version!(),
        "window": window,
        "interior_windows_per_sequence": n_interior,
        "max_period": max_period,
        "unit": unit,
        "unit_source": if unit.is_some() { unit_source } else { "none" },
        "sequences": windows.iter().filter(|w| w.kind == Kind::Start).count(),
        "telomeres": telomere_counts(&windows),
        "ends_gap_near_tip": gap_ends,
        "summary": summary,
        "reads": reads_summary,
    });
    let mut f = File::create(format!("{prefix}.caps.json"))?;
    serde_json::to_writer_pretty(&mut f, &json)?;
    writeln!(f)?;
    eprintln!("[+]\tWritten {prefix}.caps.tsv, .caps.families.tsv and .caps.json");
    eprintln!(
        "[+]\tClass: {}. Top family at ends: {} ({} bp period), top at {}/{} ends, {:.0}% of interior windows",
        summary["class"],
        summary["top_family_consensus_short"],
        summary["top_family_period"],
        summary["top_family_ends"],
        n_ends,
        summary["top_family_interior_fraction"].as_f64().unwrap_or(0.0) * 100.0
    );
    Ok(())
}

/// What the reads say about one end.
struct EndReads {
    id: String,
    kind: Kind,
    tip_reads: usize,
    /// Reads placed over the tip that leave the bait before it (not tip reads).
    diverged: usize,
    overhang_reads: usize,
    median_overhang: Option<i64>,
    state: &'static str,
    first_anchor: Option<usize>,
    /// Reads spanning [DEPTH_AT, DEPTH_AT + DEPTH_SPAN) inside the tip, and
    /// whether that is low for the assembly.
    depth: usize,
    low_coverage: bool,
    /// Overhangs analysed, their median tandem cover, and the family top in
    /// most of them (with how many).
    analysed: usize,
    tandem_cover: Option<f64>,
    family: Option<(usize, usize)>,
}

/// Anchor reads to the ends, call each end, and find the tandem repeat
/// families of the overhangs (added to `builder`). Kept overhangs are
/// written to `{prefix}.caps.overhang.fa`, and per end results to
/// `{prefix}.caps.reads.tsv`.
fn end_reads(
    reads: &[PathBuf],
    baits: &[Vec<u8>],
    ends: &[&Window],
    builder: &mut FamilyBuilder,
    max_period: usize,
    prefix: &str,
) -> Result<Vec<EndReads>> {
    let index = Baits::new(&baits.iter().map(|b| b.as_slice()).collect::<Vec<_>>());
    eprintln!("[+]	Anchoring reads to {} ends", baits.len());
    let mut overhangs: Vec<Vec<i64>> = vec![Vec::new(); baits.len()];
    let mut depth = vec![0usize; baits.len()];
    let mut diverged = vec![0usize; baits.len()];
    let mut kept: Vec<Vec<(String, Vec<u8>)>> = vec![Vec::new(); baits.len()];
    let mut n_reads = 0usize;
    for path in reads {
        let hits: Vec<(String, anchor::Hit)> = open_fasta_reader(path)?
            .records()
            .par_bridge()
            .filter_map(|r| {
                let r = r.ok()?;
                index
                    .anchor(r.seq(), MIN_OVERHANG)
                    .map(|h| (r.id().to_owned(), h))
            })
            .collect();
        for (id, h) in hits {
            if -h.overhang <= DEPTH_AT && h.span_end >= DEPTH_AT + DEPTH_SPAN {
                depth[h.end] += 1;
            }
            if h.overhang < -anchor::TIP {
                continue;
            }
            if h.diverged {
                diverged[h.end] += 1;
                continue;
            }
            n_reads += 1;
            overhangs[h.end].push(h.overhang);
            if let Some(s) = h.sequence {
                if kept[h.end].len() < OVERHANGS_PER_END {
                    kept[h.end].push((id, s));
                }
            }
        }
    }
    eprintln!(
        "[+]	{n_reads} reads cover a tip; {} diverge from the assembly before it",
        diverged.iter().sum::<usize>()
    );

    // tandem repeats in the kept overhangs, in the same families as the windows
    let flat: Vec<(usize, &[u8])> = kept
        .iter()
        .enumerate()
        .flat_map(|(e, k)| k.iter().map(move |(_, s)| (e, s.as_slice())))
        .collect();
    let arrays: Vec<Vec<Array>> = flat
        .par_iter()
        .map(|(_, s)| tandem::find_arrays(s, max_period))
        .collect();
    let seqs: Vec<&[u8]> = flat.iter().map(|(_, s)| *s).collect();
    let assignment = builder.assign_all(&arrays, &seqs);
    let mut tops: Vec<Vec<(Option<usize>, f64)>> = vec![Vec::new(); baits.len()];
    for ((arrays, fams), (e, s)) in arrays.into_iter().zip(assignment).zip(&flat) {
        let r = window_result(arrays, fams, s.len(), &builder.families);
        let top = r
            .top
            .filter(|t| t.cover >= OVERHANG_MIN_COVER)
            .map(|t| t.family);
        tops[*e].push((top, r.tandem_cover));
    }

    let mut fa = LineWriter::new(File::create(format!("{prefix}.caps.overhang.fa"))?);
    for (e, k) in kept.iter().enumerate() {
        for (id, s) in k {
            writeln!(
                fa,
                ">{id} {}:{} overhang={}",
                ends[e].id,
                ends[e].kind.name(),
                s.len()
            )?;
            fa.write_all(s)?;
            writeln!(fa)?;
        }
    }

    let median_depth = median(depth.iter().map(|&d| d as f64).collect());
    let per_end: Vec<EndReads> = (0..baits.len())
        .map(|e| {
            let ov = &mut overhangs[e];
            ov.sort_unstable();
            let n_over = ov.iter().filter(|&&o| o >= MIN_OVERHANG as i64).count();
            let state = if ov.len() < MIN_TIP_READS {
                "no_reads"
            } else if 2 * n_over >= ov.len() {
                "open"
            } else if 10 * n_over <= ov.len() {
                "closed"
            } else {
                "mixed"
            };
            let mut counts: HashMap<usize, usize> = HashMap::new();
            for f in tops[e].iter().filter_map(|t| t.0) {
                *counts.entry(f).or_default() += 1;
            }
            EndReads {
                id: ends[e].id.clone(),
                kind: ends[e].kind,
                tip_reads: ov.len(),
                diverged: diverged[e],
                overhang_reads: n_over,
                median_overhang: (!ov.is_empty()).then(|| ov[ov.len() / 2]),
                state,
                first_anchor: index.first_anchor[e],
                depth: depth[e],
                low_coverage: (depth[e] as f64) < LOW_DEPTH * median_depth,
                analysed: tops[e].len(),
                tandem_cover: (!tops[e].is_empty())
                    .then(|| median(tops[e].iter().map(|t| t.1).collect())),
                family: counts
                    .into_iter()
                    .max_by_key(|&(f, n)| (n, std::cmp::Reverse(f))),
            }
        })
        .collect();

    let mut f = LineWriter::new(File::create(format!("{prefix}.caps.reads.tsv"))?);
    writeln!(
        f,
        "id\tside\ttip_reads\toverhang_reads\tmedian_overhang\tstate\tfirst_anchor\toverhangs_analysed\toverhang_tandem_cover\toverhang_family\toverhang_family_reads\toverhang_family_period\toverhang_family_telomeric\toverhang_family_consensus\tdepth_{DEPTH_AT}\tlow_coverage\tdiverged_reads"
    )?;
    let na = |x: Option<String>| x.unwrap_or_else(|| "NA".into());
    for r in &per_end {
        let fam = match r.family {
            Some((fam, n)) => {
                let x = &builder.families[fam];
                format!(
                    "{fam}\t{n}\t{}\t{}\t{}",
                    x.period,
                    x.telomeric,
                    x.consensus.chars().take(200).collect::<String>()
                )
            }
            None => "NA\t0\tNA\tNA\tNA".into(),
        };
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{fam}\t{}\t{}\t{}",
            r.id,
            r.kind.name(),
            r.tip_reads,
            r.overhang_reads,
            na(r.median_overhang.map(|m| m.to_string())),
            r.state,
            na(r.first_anchor.map(|m| m.to_string())),
            r.analysed,
            na(r.tandem_cover.map(|c| format!("{c:.3}"))),
            r.depth,
            r.low_coverage,
            r.diverged,
        )?;
    }
    eprintln!("[+]	Written {prefix}.caps.reads.tsv and .caps.overhang.fa");
    Ok(per_end)
}

/// End states, and the overhang family top at most open ends: the repeat
/// beyond the assembly ends (the telomere in assemblies that trimmed it).
fn summarise_reads(per_end: &[EndReads], families: &[Family]) -> serde_json::Value {
    let count = |s: &str| per_end.iter().filter(|r| r.state == s).count();
    let open: Vec<&EndReads> = per_end.iter().filter(|r| r.state == "open").collect();
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for r in &open {
        if let Some((f, _)) = r.family {
            *counts.entry(f).or_default() += 1;
        }
    }
    let top = counts
        .into_iter()
        .max_by_key(|&(f, n)| (n, std::cmp::Reverse(f)));
    let covers: Vec<f64> = open.iter().filter_map(|r| r.tandem_cover).collect();
    let consensus = top
        .map(|(f, _)| families[f].consensus.clone())
        .unwrap_or_default();
    serde_json::json!({
        "ends": per_end.len(),
        "open": open.len(),
        "mixed": count("mixed"),
        "closed": count("closed"),
        "no_reads": count("no_reads"),
        "median_depth_near_tip": median(per_end.iter().map(|r| r.depth as f64).collect()),
        "low_coverage_ends": per_end.iter().filter(|r| r.low_coverage).count(),
        "tip_reads": per_end.iter().map(|r| r.tip_reads).sum::<usize>(),
        "diverged_reads": per_end.iter().map(|r| r.diverged).sum::<usize>(),
        "open_median_overhang": median(open.iter().filter_map(|r| r.median_overhang.map(|m| m as f64)).collect()),
        "open_overhang_tandem_cover_median": if covers.is_empty() { serde_json::Value::Null } else { median(covers).into() },
        "overhang_family": top.map(|(f, _)| f),
        "overhang_family_open_ends": top.map(|(_, n)| n).unwrap_or(0),
        "overhang_family_families_at_open_ends": open.iter().filter_map(|r| r.family.map(|f| f.0)).collect::<HashSet<_>>().len(),
        "overhang_family_period": top.map(|(f, _)| families[f].period).unwrap_or(0),
        "overhang_family_telomeric": top.map(|(f, _)| families[f].telomeric),
        "overhang_family_consensus": consensus,
        "overhang_family_consensus_short": consensus.chars().take(40).collect::<String>(),
    })
}

fn results_arrays(results: &[WindowResult]) -> Vec<&[Array]> {
    results.iter().map(|r| r.arrays.as_slice()).collect()
}

/// A k-mer of mostly one base, or of fewer than three bases.
fn low_complexity(mut code: u32, k: usize) -> bool {
    let mut counts = [0usize; 4];
    for _ in 0..k {
        counts[(code & 3) as usize] += 1;
        code >>= 2;
    }
    counts.iter().filter(|&&c| c > 0).count() < 3 || counts.iter().any(|&c| c * 4 > k * 3)
}

/// Non-tandem sequence shared between windows (e.g. the HeT-A/TART
/// retrotransposon arrays at *Drosophila* ends). Tandem arrays are masked;
/// for each window, the bp covered by [`SHARED_K`]-mers it shares with each
/// window of the same kind on another sequence. Partners share at least
/// [`MIN_SHARED`] bp. Orientation counts whether shared k-mers read the same
/// way round in both windows; end windows read from the tip, so arrays
/// oriented consistently towards the tip read the same way.
fn shared_sequence(windows: &[Window], arrays: &[&[Array]]) -> Vec<Shared> {
    // per window: (position, canonical k-mer, reads forward), and each
    // position's rank among unmasked bases (non-tandem coordinates)
    let (kmers, ranks): (Vec<Vec<WindowKmer>>, Vec<Vec<u32>>) = windows
        .par_iter()
        .zip(arrays)
        .map(|(w, arrays)| {
            let mut s = w.seq.clone();
            for a in arrays.iter() {
                s[a.start..a.end].fill(b'N');
            }
            let mut rank = Vec::with_capacity(s.len());
            let mut r = 0u32;
            for &b in &s {
                rank.push(r);
                r += u32::from(b != b'N');
            }
            let ks = tandem::kmer_codes(&s, SHARED_K)
                .into_iter()
                .enumerate()
                .filter_map(|(i, c)| {
                    let c = c?;
                    let rc = tandem::revcomp_code(c, SHARED_K);
                    (!low_complexity(c, SHARED_K)).then_some((i, c.min(rc), c <= rc))
                })
                .collect();
            (ks, rank)
        })
        .unzip();
    // k-mer -> (window, reads forward, first position), once per window
    type Index = HashMap<u32, Vec<(u32, bool, u32)>, BuildHasherDefault<anchor::KmerHasher>>;
    let mut index: Index = HashMap::default();
    for (w, ks) in kmers.iter().enumerate() {
        let mut seen = HashSet::new();
        for &(pos, c, f) in ks {
            if seen.insert(c) {
                index.entry(c).or_default().push((w as u32, f, pos as u32));
            }
        }
    }
    (0..windows.len())
        .into_par_iter()
        .map(|a| {
            let wa = &windows[a];
            // per partner: covered bp, end of the last k-mer counted, same, opposite
            let mut cover: HashMap<u32, [usize; 4]> = HashMap::new();
            // end windows: same-orientation hits per partner, for tip copies,
            // as (rank in a, rank in b, position in a, position in b)
            let mut hits: HashMap<u32, Vec<TipHit>> = HashMap::new();
            for &(pos, c, f) in &kmers[a] {
                for &(b, fb, pb) in index.get(&c).into_iter().flatten() {
                    let wb = &windows[b as usize];
                    if b as usize == a || wb.kind.is_end() != wa.kind.is_end() || wb.id == wa.id {
                        continue;
                    }
                    let e = cover.entry(b).or_default();
                    let add = (pos + SHARED_K).saturating_sub(pos.max(e[1]));
                    e[0] += add;
                    e[1] = pos + SHARED_K;
                    e[if f == fb { 2 } else { 3 }] += add;
                    if wa.kind.is_end() && f == fb {
                        let pb = pb as usize;
                        hits.entry(b).or_default().push(TipHit {
                            ra: ranks[a][pos] as usize,
                            rb: ranks[b as usize][pb] as usize,
                            pa: pos,
                            pb,
                        });
                    }
                }
            }
            let mut shared =
                cover
                    .values()
                    .filter(|e| e[0] >= MIN_SHARED)
                    .fold(Shared::default(), |s, e| Shared {
                        partners: s.partners + 1,
                        same: s.same + e[2],
                        opposite: s.opposite + e[3],
                        ..s
                    });
            for h in hits.values() {
                if let Some(extent) = tip_copy(h) {
                    shared.tip_partners += 1;
                    shared.tip_extent = shared.tip_extent.max(extent);
                }
            }
            shared
        })
        .collect()
}

/// A window k-mer: position, canonical code, reads forward.
type WindowKmer = (usize, u32, bool);

/// A k-mer shared by two end windows in the same orientation: its rank among
/// unmasked (non-tandem) bases and its position, in each window.
struct TipHit {
    ra: usize,
    rb: usize,
    pa: usize,
    pb: usize,
}

/// The extent (bp from the tip) of a segment copied from the tip between two
/// end windows, from their shared k-mers in order of window a. Chains run in
/// non-tandem coordinates, so tandem arrays of different length in the two
/// copies cost nothing: increasing in both, starting within
/// [`TIP_COPY_START`] non-tandem bp of both tips, with gaps up to
/// [`TIP_COPY_GAP`]. A chain counts if its k-mer cover is >= [`TIP_COPY_MIN`]
/// and >= [`TIP_COPY_DENSITY`] of its non-tandem extent.
fn tip_copy(hits: &[TipHit]) -> Option<usize> {
    const MAX_CHAINS: usize = 8;
    struct Chain {
        ra: usize,
        rb: usize,
        cover: usize,
        extent: usize,
    }
    let mut chains: Vec<Chain> = Vec::new();
    for h in hits {
        let mut extended = false;
        for c in chains.iter_mut() {
            if h.ra > c.ra
                && h.rb > c.rb
                && h.ra - c.ra <= TIP_COPY_GAP
                && h.rb - c.rb <= TIP_COPY_GAP
            {
                c.cover += (h.ra - c.ra).min(SHARED_K);
                (c.ra, c.rb) = (h.ra, h.rb);
                c.extent = h.pa.min(h.pb) + SHARED_K;
                extended = true;
            }
        }
        if !extended
            && h.ra <= TIP_COPY_START
            && h.rb <= TIP_COPY_START
            && chains.len() < MAX_CHAINS
        {
            chains.push(Chain {
                ra: h.ra,
                rb: h.rb,
                cover: SHARED_K,
                extent: h.pa.min(h.pb) + SHARED_K,
            });
        }
    }
    chains
        .into_iter()
        .filter(|c| {
            c.cover >= TIP_COPY_MIN
                && c.cover as f64 >= TIP_COPY_DENSITY * (c.ra.min(c.rb) + SHARED_K) as f64
        })
        .map(|c| c.extent)
        .max()
}

/// Tandem arrays, families and the end-enrichment summary for a set of
/// windows (end windows read from the tip inwards).
fn analyse(
    windows: &[Window],
    unit: Option<&str>,
    max_period: usize,
) -> (Vec<WindowResult>, FamilyBuilder, serde_json::Value) {
    let arrays: Vec<Vec<Array>> = windows
        .par_iter()
        .map(|w| tandem::find_arrays(&w.seq, max_period))
        .collect();
    let mut builder = FamilyBuilder::new(unit);
    let seqs: Vec<&[u8]> = windows.iter().map(|w| w.seq.as_slice()).collect();
    let assignment = builder.assign_all(&arrays, &seqs);
    let results: Vec<WindowResult> = arrays
        .into_iter()
        .zip(assignment)
        .zip(windows)
        .map(|((arrays, fams), w)| window_result(arrays, fams, w.seq.len(), &builder.families))
        .collect();
    let shared = shared_sequence(windows, &results_arrays(&results));
    let mut results = results;
    for (r, sh) in results.iter_mut().zip(shared) {
        r.shared = sh;
    }
    let summary = summarise(windows, &results, &builder.families);
    (results, builder, summary)
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// The most abundant candidate telomeric repeat unit at the ends, if any.
fn discover_unit(windows: &[Window]) -> Result<Option<String>> {
    let ends: Vec<(String, &[u8])> = windows
        .iter()
        .filter(|w| w.kind.is_end())
        .map(|w| (format!("{}:{}", w.id, w.kind.name()), w.seq.as_slice()))
        .collect();
    Ok(explore::top_unit(
        &ends,
        DISCOVERY_MINIMUM,
        DISCOVERY_MAXIMUM,
        DISCOVERY_THRESHOLD,
        true,
    )?
    .map(|(u, _)| u))
}

/// Families of tandem arrays, built up array by array.
struct FamilyBuilder {
    families: Vec<Family>,
    /// Canonical consensus k-mers to families.
    index: HashMap<u32, Vec<usize>>,
    /// Canonical repeated k-mers to families.
    index_repeated: HashMap<u32, Vec<usize>>,
    unit_kmers: HashSet<u32>,
}

fn canonical(set: &HashSet<u32>) -> HashSet<u32> {
    set.iter()
        .map(|&c| tandem::canonical_code(c, FAMILY_K))
        .collect()
}

/// Jaccard similarity of a k-mer set to each family sharing any k-mer.
fn similar(
    set: &HashSet<u32>,
    index: &HashMap<u32, Vec<usize>>,
    size: impl Fn(usize) -> usize,
) -> HashMap<usize, f64> {
    let mut shared: HashMap<usize, usize> = HashMap::new();
    for c in set {
        for &f in index.get(c).into_iter().flatten() {
            *shared.entry(f).or_default() += 1;
        }
    }
    shared
        .into_iter()
        .map(|(f, s)| (f, s as f64 / (set.len() + size(f) - s) as f64))
        .collect()
}

impl FamilyBuilder {
    /// With a unit, family 0 is the telomere.
    fn new(unit: Option<&str>) -> Self {
        let mut b = FamilyBuilder {
            families: Vec::new(),
            index: HashMap::new(),
            index_repeated: HashMap::new(),
            unit_kmers: HashSet::new(),
        };
        if let Some(u) = unit {
            let forward = tandem::unit_kmers(u);
            let canon = canonical(&forward);
            for &c in &canon {
                b.index.entry(c).or_default().push(0);
                b.index_repeated.entry(c).or_default().push(0);
            }
            b.families.push(Family {
                period: u.len(),
                consensus: u.to_string(),
                kmers: canon.clone(),
                repeated: canon.clone(),
                forward,
                telomeric: true,
            });
            b.unit_kmers = canon;
        }
        b
    }

    /// Is an array mostly made of the unit's k-mers (variants and all)?
    fn is_telomeric(&self, a: &Array, seq: &[u8]) -> bool {
        if self.unit_kmers.is_empty() {
            return false;
        }
        let codes = tandem::kmer_codes(&seq[a.start..a.end], FAMILY_K);
        let hits = codes
            .iter()
            .flatten()
            .filter(|&&c| {
                self.unit_kmers
                    .contains(&tandem::canonical_code(c, FAMILY_K))
            })
            .count();
        !codes.is_empty() && hits as f64 >= TELOMERIC_COVER * codes.len() as f64
    }

    /// The family of an array from `seq`, founding a new one if none is similar.
    fn assign(&mut self, a: &Array, seq: &[u8]) -> usize {
        if self.is_telomeric(a, seq) {
            return 0;
        }
        let canon = canonical(&a.kmers);
        let repeated = canonical(&a.repeated);
        // the better of the two similarities, consensus and repeated k-mers
        let families = &self.families;
        let mut best = similar(&canon, &self.index, |f| families[f].kmers.len());
        for (f, j) in similar(&repeated, &self.index_repeated, |f| {
            families[f].repeated.len()
        }) {
            let e = best.entry(f).or_default();
            *e = e.max(j);
        }
        if let Some((f, _)) = best
            .into_iter()
            .filter(|&(_, j)| j >= FAMILY_JACCARD)
            .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)))
        {
            return f;
        }
        let f = self.families.len();
        for &c in &canon {
            self.index.entry(c).or_default().push(f);
        }
        for &c in &repeated {
            self.index_repeated.entry(c).or_default().push(f);
        }
        let telomeric = !self.unit_kmers.is_empty()
            && canon.intersection(&self.unit_kmers).count() as f64
                >= TELOMERIC_SHARE * self.unit_kmers.len() as f64;
        self.families.push(Family {
            period: a.period,
            consensus: utils::lex_min(&String::from_utf8_lossy(&a.consensus)),
            kmers: canon,
            repeated,
            forward: a.kmers.clone(),
            telomeric,
        });
        f
    }

    /// Assign the arrays of a set of sequences, largest arrays first, so each
    /// family is represented by its largest array. Returns, per sequence, the
    /// family of each array.
    fn assign_all(&mut self, arrays: &[Vec<Array>], seqs: &[&[u8]]) -> Vec<Vec<usize>> {
        let mut order: Vec<(usize, usize)> = arrays
            .iter()
            .enumerate()
            .flat_map(|(w, a)| (0..a.len()).map(move |i| (w, i)))
            .collect();
        order.sort_unstable_by_key(|&(w, i)| (std::cmp::Reverse(arrays[w][i].len()), w, i));
        let mut assignment: Vec<Vec<usize>> = arrays.iter().map(|a| vec![0; a.len()]).collect();
        for (w, i) in order {
            assignment[w][i] = self.assign(&arrays[w][i], seqs[w]);
        }
        assignment
    }
}

/// Does an array read the same way round as its family's representative?
fn reads_forward(a: &Array, family: &Family) -> bool {
    let fwd = a
        .kmers
        .iter()
        .filter(|c| family.forward.contains(c))
        .count();
    let rev = a
        .kmers
        .iter()
        .filter(|&&c| family.forward.contains(&tandem::revcomp_code(c, FAMILY_K)))
        .count();
    fwd >= rev
}

fn window_result(
    arrays: Vec<Array>,
    fams: Vec<usize>,
    window_len: usize,
    families: &[Family],
) -> WindowResult {
    let mut bp: HashMap<usize, usize> = HashMap::new();
    for (a, &f) in arrays.iter().zip(&fams) {
        *bp.entry(f).or_default() += a.len();
    }
    let top = bp
        .iter()
        .max_by_key(|(f, b)| (**b, std::cmp::Reverse(**f)))
        .map(|(&f, &b)| {
            // the family's longest array in the window gives period and orientation
            let longest = arrays
                .iter()
                .zip(&fams)
                .filter(|(_, &g)| g == f)
                .map(|(a, _)| a)
                .max_by_key(|a| a.len())
                .unwrap();
            TopFamily {
                family: f,
                cover: (b as f64 / window_len as f64).min(1.0),
                period: longest.period,
                forward: reads_forward(longest, &families[f]),
            }
        });
    let tandem_cover =
        (arrays.iter().map(|a| a.len()).sum::<usize>() as f64 / window_len as f64).min(1.0);
    WindowResult {
        arrays,
        families: fams,
        top,
        tandem_cover,
        shared: Shared::default(),
    }
}

/// Assembly gaps (runs of >= [`MIN_GAP_RUN`] Ns) in a window: total gap bp
/// and the start of the first gap (from the tip, for end windows).
fn gaps(seq: &[u8]) -> (usize, Option<usize>) {
    let (mut total, mut first, mut run) = (0, None, 0);
    for (i, &b) in seq.iter().chain(std::iter::once(&b'A')).enumerate() {
        if b == b'N' {
            run += 1;
        } else {
            if run >= MIN_GAP_RUN {
                total += run;
                first.get_or_insert(i - run);
            }
            run = 0;
        }
    }
    (total, first)
}

fn median(mut v: Vec<f64>) -> f64 {
    if v.is_empty() {
        return 0.0;
    }
    v.sort_unstable_by(|a, b| a.total_cmp(b));
    v[v.len() / 2]
}

/// The end-enrichment summary and class.
fn summarise(
    windows: &[Window],
    results: &[WindowResult],
    families: &[Family],
) -> serde_json::Value {
    let ends: Vec<usize> = (0..windows.len())
        .filter(|&i| windows[i].kind.is_end())
        .collect();
    let interior: Vec<usize> = (0..windows.len())
        .filter(|&i| !windows[i].kind.is_end())
        .collect();
    let top_of = |i: usize| results[i].top.as_ref().map(|t| t.family);

    let mut counts: HashMap<usize, usize> = HashMap::new();
    for &i in &ends {
        if let Some(f) = top_of(i) {
            *counts.entry(f).or_default() += 1;
        }
    }
    let top = counts
        .iter()
        .max_by_key(|(f, n)| (**n, std::cmp::Reverse(**f)))
        .map(|(&f, &n)| (f, n));
    let frac = |n: usize, d: usize| if d == 0 { 0.0 } else { n as f64 / d as f64 };
    // shared non-tandem sequence
    let sharing = |ix: &[usize]| {
        frac(
            ix.iter()
                .filter(|&&i| results[i].shared.partners >= SHARED_MIN_PARTNERS)
                .count(),
            ix.len(),
        )
    };
    let (sh_end, sh_int) = (sharing(&ends), sharing(&interior));
    let shared = sh_end >= SHARED_END_FRACTION && sh_end >= CAP_ENRICHMENT * sh_int;
    let (same, opposite) = ends.iter().fold((0, 0), |(s, o), &i| {
        (s + results[i].shared.same, o + results[i].shared.opposite)
    });
    // segments copied from the tip between ends
    let tip_ends: Vec<usize> = ends
        .iter()
        .copied()
        .filter(|&i| results[i].shared.tip_partners > 0)
        .collect();
    // do the two ends of a sequence share their top family? Expected from
    // family frequencies among ends with a top family (random pairing)
    let mut by_seq: HashMap<&str, [Option<usize>; 2]> = HashMap::new();
    for &i in &ends {
        let k = usize::from(windows[i].kind == Kind::End);
        by_seq.entry(windows[i].id.as_str()).or_default()[k] = top_of(i);
    }
    let both: Vec<(usize, usize)> = by_seq
        .values()
        .filter_map(|p| Some((p[0]?, p[1]?)))
        .collect();
    let mut fam_counts: HashMap<usize, usize> = HashMap::new();
    for &(a, b) in &both {
        *fam_counts.entry(a).or_default() += 1;
        *fam_counts.entry(b).or_default() += 1;
    }
    let n_top = 2 * both.len();
    let p_same = if n_top < 2 {
        0.0
    } else {
        fam_counts
            .values()
            .map(|&c| (c * (c - 1)) as f64)
            .sum::<f64>()
            / (n_top * (n_top - 1)) as f64
    };
    let shared_json = serde_json::json!({
        "shared_ends": shared,
        "ends_shared_fraction": sh_end,
        "interior_shared_fraction": sh_int,
        "ends_shared_same_orientation": if same + opposite == 0 { serde_json::Value::Null } else { frac(same, same + opposite).into() },
        "tip_copied_ends": tip_ends.len(),
        "tip_copied_pairs": tip_ends.iter().map(|&i| results[i].shared.tip_partners).sum::<usize>() / 2,
        "tip_copied_max_extent": tip_ends.iter().map(|&i| results[i].shared.tip_extent).max().unwrap_or(0),
        "sequences_both_ends_with_top_family": both.len(),
        "sequences_same_top_family_both_ends": both.iter().filter(|(a, b)| a == b).count(),
        "sequences_same_top_family_expected": p_same * both.len() as f64,
    });
    let ends_cover = median(ends.iter().map(|&i| results[i].tandem_cover).collect());
    let interior_cover = median(interior.iter().map(|&i| results[i].tandem_cover).collect());

    let Some((f, n)) = top else {
        return merge(
            serde_json::json!({
                "class": if shared { "shared_complex" } else { "none" },
                "top_family": null,
                "top_family_consensus_short": "",
                "top_family_period": 0,
                "top_family_ends": 0,
                "top_family_interior_fraction": 0.0,
                "ends_tandem_cover_median": ends_cover,
                "interior_tandem_cover_median": interior_cover,
            }),
            shared_json,
        );
    };
    let at: Vec<usize> = ends
        .iter()
        .copied()
        .filter(|&i| top_of(i) == Some(f))
        .collect();
    let f_end = frac(n, ends.len());
    let f_int = frac(
        interior.iter().filter(|&&i| top_of(i) == Some(f)).count(),
        interior.len(),
    );
    let cover = median(
        at.iter()
            .map(|&i| results[i].top.as_ref().unwrap().cover)
            .collect(),
    );
    let forward = at
        .iter()
        .filter(|&&i| results[i].top.as_ref().unwrap().forward)
        .count();
    let cap =
        f_end >= CAP_END_FRACTION && cover >= CAP_MEDIAN_COVER && f_end >= CAP_ENRICHMENT * f_int;
    let class = match (cap, families[f].telomeric, shared) {
        (true, true, _) => "telomeric",
        (true, false, _) => "terminal_satellite",
        (false, _, true) => "shared_complex",
        (false, _, false) => "none",
    };
    let consensus = &families[f].consensus;
    merge(
        serde_json::json!({
            "class": class,
            "top_family": f,
            "top_family_period": families[f].period,
            "top_family_consensus": consensus,
            "top_family_consensus_short": consensus.chars().take(40).collect::<String>(),
            "top_family_telomeric": families[f].telomeric,
            "top_family_ends": n,
            "top_family_end_fraction": f_end,
            "top_family_interior_fraction": f_int,
            "top_family_median_cover": cover,
            "top_family_orientation_consistency": frac(forward.max(n - forward), n),
            "ends": ends.len(),
            "interior_windows": interior.len(),
            "ends_tandem_cover_median": ends_cover,
            "interior_tandem_cover_median": interior_cover,
        }),
        shared_json,
    )
}

/// The fields of two JSON objects in one.
fn merge(mut a: serde_json::Value, b: serde_json::Value) -> serde_json::Value {
    if let (Some(a), serde_json::Value::Object(b)) = (a.as_object_mut(), b) {
        a.extend(b);
    }
    a
}

fn telomere_counts(windows: &[Window]) -> serde_json::Value {
    let mut counts: HashMap<&str, usize> = HashMap::new();
    for w in windows {
        if let Some(t) = &w.telomere {
            *counts.entry(t.status.name()).or_default() += 1;
        }
    }
    serde_json::json!(counts)
}

fn write_windows(
    windows: &[Window],
    results: &[WindowResult],
    families: &[Family],
    path: &str,
) -> Result<()> {
    let mut f = LineWriter::new(File::create(path)?);
    writeln!(
        f,
        "id\tsequence_length\tkind\twindow_start\ttelomere\ttelomere_length\ttandem_cover\tarrays\tshared_partners\ttop_family\ttop_family_period\ttop_family_cover\ttop_family_telomeric\ttop_family_orientation\ttop_family_consensus\tgap_bp\tfirst_gap\ttip_copy_partners\ttip_copy_extent"
    )?;
    for (w, r) in windows.iter().zip(results) {
        let (status, tlen) = match &w.telomere {
            Some(t) => (
                t.status.name().to_string(),
                t.telomere
                    .as_ref()
                    .filter(|_| t.status != Status::Absent)
                    .map(|x| (x.end - x.start).to_string())
                    .unwrap_or("NA".into()),
            ),
            None => ("NA".into(), "NA".into()),
        };
        let top = match &r.top {
            Some(t) => {
                let fam = &families[t.family];
                format!(
                    "{}\t{}\t{:.3}\t{}\t{}\t{}",
                    t.family,
                    t.period,
                    t.cover,
                    fam.telomeric,
                    if !w.kind.is_end() {
                        "NA"
                    } else if t.forward {
                        "+"
                    } else {
                        "-"
                    },
                    fam.consensus.chars().take(200).collect::<String>()
                )
            }
            None => "NA\tNA\t0\tNA\tNA\tNA".into(),
        };
        let (gap_bp, first_gap) = gaps(&w.seq);
        let first_gap = first_gap.map_or("NA".into(), |g| g.to_string());
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{status}\t{tlen}\t{:.3}\t{}\t{}\t{top}\t{gap_bp}\t{first_gap}\t{}\t{}",
            w.id,
            w.sequence_length,
            w.kind.name(),
            w.position,
            r.tandem_cover,
            r.arrays.len(),
            r.shared.partners,
            r.shared.tip_partners,
            r.shared.tip_extent,
        )?;
    }
    Ok(())
}

/// Families top at, or present in, at least one end window, by ends where top.
fn write_families(
    windows: &[Window],
    results: &[WindowResult],
    families: &[Family],
    path: &str,
) -> Result<()> {
    // per family: ends top, interior top, ends present, interior present, ends top forward
    let mut stats = vec![[0usize; 5]; families.len()];
    for (w, r) in windows.iter().zip(results) {
        let e = w.kind.is_end();
        if let Some(t) = &r.top {
            stats[t.family][if e { 0 } else { 1 }] += 1;
            if e && t.forward {
                stats[t.family][4] += 1;
            }
        }
        for &fam in r.families.iter().collect::<HashSet<_>>() {
            stats[fam][if e { 2 } else { 3 }] += 1;
        }
    }
    let mut order: Vec<usize> = (0..families.len())
        .filter(|&f| stats[f][0] + stats[f][2] > 0)
        .collect();
    order.sort_unstable_by_key(|&f| {
        (
            std::cmp::Reverse(stats[f][0]),
            std::cmp::Reverse(stats[f][2]),
            f,
        )
    });
    let n_ends = windows.iter().filter(|w| w.kind.is_end()).count();
    let n_int = windows.len() - n_ends;
    let mut out = LineWriter::new(File::create(path)?);
    writeln!(
        out,
        "family\tperiod\ttelomeric\tends_top\tinterior_top\tends_present\tinterior_present\tends\tinterior_windows\tends_top_same_orientation\tconsensus"
    )?;
    for f in order {
        let s = stats[f];
        writeln!(
            out,
            "{f}\t{}\t{}\t{}\t{}\t{}\t{}\t{n_ends}\t{n_int}\t{}\t{}",
            families[f].period,
            families[f].telomeric,
            s[0],
            s[1],
            s[2],
            s[3],
            s[4].max(s[0] - s[4]),
            families[f].consensus
        )?;
    }
    Ok(())
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

    /// Windows for `n` chromosomes, the ends from `cap(chromosome)` followed
    /// by random sequence, read from the tip inwards; one random interior window each.
    fn windows(n: usize, cap: impl Fn(usize) -> Vec<u8>) -> Vec<Window> {
        windows_ends(n, &cap, &cap)
    }

    /// As [`windows`], with different caps at the start and end of each chromosome.
    fn windows_ends(
        n: usize,
        start: impl Fn(usize) -> Vec<u8>,
        end: impl Fn(usize) -> Vec<u8>,
    ) -> Vec<Window> {
        let mut out = Vec::new();
        for c in 0..n {
            for (k, kind) in [Kind::Start, Kind::End, Kind::Interior]
                .into_iter()
                .enumerate()
            {
                let mut seq = match kind {
                    Kind::Start => start(c),
                    Kind::End => end(c),
                    Kind::Interior => Vec::new(),
                };
                let fill = 20_000 - seq.len();
                seq.extend(random_seq(fill, (c * 3 + k) as u64 + 100));
                out.push(Window {
                    id: format!("chr{c}"),
                    sequence_length: 1_000_000,
                    kind,
                    position: 0,
                    seq,
                    telomere: None,
                });
            }
        }
        out
    }

    #[test]
    fn test_telomeric_ends() {
        let w = windows(6, |_| b"CCCTAA".repeat(1000));
        let (_, _, s) = analyse(&w, Some("AACCCT"), 2000);
        assert_eq!(s["class"], "telomeric");
        assert_eq!(s["top_family_ends"], 12);
        assert_eq!(s["top_family_orientation_consistency"], 1.0);
    }

    #[test]
    fn test_terminal_satellite() {
        // the same 180 bp satellite at every end, slightly different copies
        let unit = random_seq(180, 1);
        let w = windows(6, |c| {
            let mut s = Vec::new();
            for i in 0..30 {
                let mut copy = unit.clone();
                copy[(c * 7 + i * 13) % 180] = b'A';
                s.extend(copy);
            }
            s
        });
        let (_, _, s) = analyse(&w, Some("AACCCT"), 2000);
        assert_eq!(s["class"], "terminal_satellite");
        assert_eq!(s["top_family_period"], 180);
        assert_eq!(s["top_family_ends"], 12);
    }

    #[test]
    fn test_variant_rich_telomere_is_one_family() {
        // plant-like telomeres mixing CCCTAAA and CCCCTAAA copies
        let w = windows(6, |c| {
            (0..800)
                .flat_map(|i| {
                    if (i * 7 + c) % 3 == 0 {
                        b"CCCCTAAA".to_vec()
                    } else {
                        b"CCCTAAA".to_vec()
                    }
                })
                .collect()
        });
        let (_, b, s) = analyse(&w, Some("AAACCCT"), 2000);
        assert_eq!(s["class"], "telomeric");
        assert_eq!(s["top_family"], 0);
        assert_eq!(s["top_family_ends"], 12);
        assert_eq!(b.families[0].consensus, "AAACCCT");
    }

    #[test]
    fn test_shared_non_tandem_sequence() {
        // a 6 kb non-tandem element, in the same orientation, at every end
        let element = random_seq(6000, 3);
        let w = windows(6, |_| element.clone());
        let (r, _, s) = analyse(&w, None, 2000);
        assert_eq!(s["class"], "shared_complex");
        assert_eq!(s["ends_shared_fraction"], 1.0);
        assert_eq!(s["interior_shared_fraction"], 0.0);
        assert!(s["ends_shared_same_orientation"].as_f64().unwrap() > 0.99);
        assert_eq!(r[0].shared.partners, 10);
    }

    /// `copies` copies of a random unit, each with one substitution.
    fn satellite(seed: u64, period: usize, copies: usize, salt: usize) -> Vec<u8> {
        let unit = random_seq(period, seed);
        (0..copies)
            .flat_map(|i| {
                let mut copy = unit.clone();
                copy[(salt * 7 + i * 13) % period] = b'A';
                copy
            })
            .collect()
    }

    #[test]
    fn test_end_pairing() {
        // one satellite at every start, another at every end (as spider
        // centromeric vs distal ends): the two ends never share a family
        let w = windows_ends(
            6,
            |c| satellite(1, 180, 30, c),
            |c| satellite(2, 230, 30, c),
        );
        let (_, _, s) = analyse(&w, None, 2000);
        assert_eq!(s["sequences_both_ends_with_top_family"], 6);
        assert_eq!(s["sequences_same_top_family_both_ends"], 0);
        let expected = s["sequences_same_top_family_expected"].as_f64().unwrap();
        assert!((expected - 6.0 * 60.0 / 132.0).abs() < 1e-9);
        // a telomere at every end: all pairs share it
        let w = windows(6, |_| b"CCCTAA".repeat(1000));
        let (_, _, s) = analyse(&w, Some("AACCCT"), 2000);
        assert_eq!(s["sequences_same_top_family_both_ends"], 6);
    }

    #[test]
    fn test_tip_copy() {
        // the start of chr0 copied onto the end of chr3: 12 kb from the tip,
        // ~0.5% substitutions, and an array of different length in each copy
        let segment = random_seq(12_000, 7);
        let copy = |len: usize, salt: usize| {
            let mut s = segment[..6000].to_vec();
            s.extend(satellite(8, 170, len, salt));
            s.extend(&segment[6000..]);
            for i in (salt * 37..s.len()).step_by(211) {
                s[i] = b"ACGT"[(i + 1) % 4];
            }
            s
        };
        let w = windows_ends(
            6,
            |c| {
                if c == 0 {
                    copy(10, 1)
                } else {
                    random_seq(3000, 60 + c as u64)
                }
            },
            |c| {
                if c == 3 {
                    copy(14, 2)
                } else {
                    random_seq(3000, 80 + c as u64)
                }
            },
        );
        let (r, _, s) = analyse(&w, None, 2000);
        assert_eq!(s["tip_copied_ends"], 2);
        assert_eq!(s["tip_copied_pairs"], 1);
        assert_eq!(
            (r[0].shared.tip_partners, r[10].shared.tip_partners),
            (1, 1)
        );
        let extent = s["tip_copied_max_extent"].as_u64().unwrap();
        assert!((11_000..=14_500).contains(&extent), "{extent}");
        // shared but not from the tip: the element 8 kb inside one end
        let element = random_seq(8000, 9);
        let w = windows_ends(
            6,
            |c| {
                if c == 0 {
                    element.clone()
                } else {
                    random_seq(3000, 60 + c as u64)
                }
            },
            |c| {
                if c == 3 {
                    [random_seq(8000, 10), element.clone()].concat()
                } else {
                    random_seq(3000, 80 + c as u64)
                }
            },
        );
        let (r, _, s) = analyse(&w, None, 2000);
        assert_eq!(r[0].shared.partners, 1);
        assert_eq!(s["tip_copied_ends"], 0);
    }

    #[test]
    fn test_gaps() {
        let mut seq = random_seq(10_000, 9);
        seq[3000..3005].fill(b'N'); // too short for a gap
        seq[4000..4100].fill(b'N');
        seq[9950..].fill(b'N');
        assert_eq!(gaps(&seq), (150, Some(4000)));
        assert_eq!(gaps(&random_seq(1000, 10)), (0, None));
    }

    #[test]
    fn test_no_cap() {
        // a different satellite at each end
        let w = windows(6, |c| random_seq(150, c as u64 + 50).repeat(20));
        let (_, _, s) = analyse(&w, None, 2000);
        assert_eq!(s["class"], "none");
    }

    #[test]
    fn test_satellite_also_interior_is_not_a_cap() {
        let unit = random_seq(200, 2);
        let mut w = windows(6, |_| unit.repeat(20));
        for x in w.iter_mut().filter(|x| x.kind == Kind::Interior) {
            x.seq[..4000].copy_from_slice(&unit.repeat(20));
        }
        let (_, _, s) = analyse(&w, None, 2000);
        assert_eq!(s["class"], "none");
        assert_eq!(s["top_family_interior_fraction"], 1.0);
    }
}
