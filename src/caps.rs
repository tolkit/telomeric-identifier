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

use crate::ends::{self, EndsConfig, Side, Status};
use crate::tandem::{self, Array, FAMILY_K};
use crate::{explore, open_fasta_reader, utils};
use anyhow::{bail, Context, Result};
use clap::crate_version;
use rayon::prelude::*;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::fs::{create_dir_all, File};
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
/// Telomere calls: minimum length and offset from the tip, as in `tidk ends`.
const TELOMERE_MIN_LENGTH: usize = 200;
const TELOMERE_MAX_OFFSET: usize = 1000;
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

    // windows; telomeres are called later, once the unit is known
    let mut windows = Vec::new();
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

    let (results, families, summary) = analyse(&windows, unit.as_deref(), max_period);

    create_dir_all(outdir)?;
    let prefix = format!("{}/{}", outdir.display(), output.display());
    write_windows(&windows, &results, &families, &format!("{prefix}.caps.tsv"))?;
    write_families(
        &windows,
        &results,
        &families,
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
        "summary": summary,
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

/// Tandem arrays, families and the end-enrichment summary for a set of
/// windows (end windows read from the tip inwards).
fn analyse(
    windows: &[Window],
    unit: Option<&str>,
    max_period: usize,
) -> (Vec<WindowResult>, Vec<Family>, serde_json::Value) {
    let arrays: Vec<Vec<Array>> = windows
        .par_iter()
        .map(|w| tandem::find_arrays(&w.seq, max_period))
        .collect();
    let unit_kmers: HashSet<u32> = unit
        .map(|u| {
            tandem::unit_kmers(u)
                .into_iter()
                .map(|c| tandem::canonical_code(c, FAMILY_K))
                .collect()
        })
        .unwrap_or_default();
    // arrays mostly made of the unit's k-mers, variants and all
    let telomeric: Vec<Vec<bool>> = arrays
        .iter()
        .zip(windows)
        .map(|(arrays, w)| {
            arrays
                .iter()
                .map(|a| {
                    let codes = tandem::kmer_codes(&w.seq[a.start..a.end], FAMILY_K);
                    let hits = codes
                        .iter()
                        .flatten()
                        .filter(|&&c| unit_kmers.contains(&tandem::canonical_code(c, FAMILY_K)))
                        .count();
                    !codes.is_empty() && hits as f64 >= TELOMERIC_COVER * codes.len() as f64
                })
                .collect()
        })
        .collect();
    let (families, assignment) = assign_families(&arrays, &telomeric, unit);
    let results: Vec<WindowResult> = arrays
        .into_iter()
        .zip(assignment)
        .zip(windows)
        .map(|((arrays, fams), w)| window_result(arrays, fams, w.seq.len(), &families))
        .collect();
    let summary = summarise(windows, &results, &families);
    (results, families, summary)
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

/// Group all arrays into families, largest arrays first, so each family is
/// represented by its largest array. Returns the families and, per window,
/// the family of each array.
fn assign_families(
    arrays: &[Vec<Array>],
    telomeric: &[Vec<bool>],
    unit: Option<&str>,
) -> (Vec<Family>, Vec<Vec<usize>>) {
    let mut order: Vec<(usize, usize)> = arrays
        .iter()
        .enumerate()
        .flat_map(|(w, a)| (0..a.len()).map(move |i| (w, i)))
        .collect();
    order.sort_unstable_by_key(|&(w, i)| (std::cmp::Reverse(arrays[w][i].len()), w, i));

    let canonical = |set: &HashSet<u32>| -> HashSet<u32> {
        set.iter()
            .map(|&c| tandem::canonical_code(c, FAMILY_K))
            .collect()
    };
    // Jaccard similarity to each family sharing any k-mer, via an index
    let similar = |set: &HashSet<u32>,
                   index: &HashMap<u32, Vec<usize>>,
                   size: &dyn Fn(usize) -> usize|
     -> HashMap<usize, f64> {
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
    };

    let mut families: Vec<Family> = Vec::new();
    let mut index: HashMap<u32, Vec<usize>> = HashMap::new();
    let mut index_repeated: HashMap<u32, Vec<usize>> = HashMap::new();
    let mut assignment: Vec<Vec<usize>> = arrays.iter().map(|a| vec![0; a.len()]).collect();
    // with a unit, family 0 is the telomere, holding every telomeric array
    let unit_kmers: HashSet<u32> = match unit {
        Some(u) => {
            let forward = tandem::unit_kmers(u);
            let canon = canonical(&forward);
            for &c in &canon {
                index.entry(c).or_default().push(0);
                index_repeated.entry(c).or_default().push(0);
            }
            families.push(Family {
                period: u.len(),
                consensus: u.to_string(),
                kmers: canon.clone(),
                repeated: canon.clone(),
                forward,
                telomeric: true,
            });
            canon
        }
        None => HashSet::new(),
    };
    for (w, i) in order {
        let a = &arrays[w][i];
        if telomeric[w][i] {
            assignment[w][i] = 0;
            continue;
        }
        let canon = canonical(&a.kmers);
        let repeated = canonical(&a.repeated);
        // the better of the two similarities, consensus and repeated k-mers
        let mut best: HashMap<usize, f64> = similar(&canon, &index, &|f| families[f].kmers.len());
        for (f, j) in similar(&repeated, &index_repeated, &|f| families[f].repeated.len()) {
            let e = best.entry(f).or_default();
            *e = e.max(j);
        }
        let best = best
            .into_iter()
            .filter(|&(_, j)| j >= FAMILY_JACCARD)
            .max_by(|a, b| a.1.total_cmp(&b.1).then(b.0.cmp(&a.0)));
        let f = match best {
            Some((f, _)) => f,
            None => {
                let f = families.len();
                for &c in &canon {
                    index.entry(c).or_default().push(f);
                }
                for &c in &repeated {
                    index_repeated.entry(c).or_default().push(f);
                }
                let telomeric = !unit_kmers.is_empty()
                    && canon.intersection(&unit_kmers).count() as f64
                        >= TELOMERIC_SHARE * unit_kmers.len() as f64;
                families.push(Family {
                    period: a.period,
                    consensus: utils::lex_min(&String::from_utf8_lossy(&a.consensus)),
                    kmers: canon,
                    repeated,
                    forward: a.kmers.clone(),
                    telomeric,
                });
                f
            }
        };
        assignment[w][i] = f;
    }
    (families, assignment)
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
    }
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
    let ends_cover = median(ends.iter().map(|&i| results[i].tandem_cover).collect());
    let interior_cover = median(interior.iter().map(|&i| results[i].tandem_cover).collect());

    let Some((f, n)) = top else {
        return serde_json::json!({
            "class": "none",
            "top_family": null,
            "top_family_consensus_short": "",
            "top_family_period": 0,
            "top_family_ends": 0,
            "top_family_interior_fraction": 0.0,
            "ends_tandem_cover_median": ends_cover,
            "interior_tandem_cover_median": interior_cover,
        });
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
    let class = match (cap, families[f].telomeric) {
        (true, true) => "telomeric",
        (true, false) => "terminal_satellite",
        (false, _) => "none",
    };
    let consensus = &families[f].consensus;
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
    })
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
        "id\tsequence_length\tkind\twindow_start\ttelomere\ttelomere_length\ttandem_cover\tarrays\ttop_family\ttop_family_period\ttop_family_cover\ttop_family_telomeric\ttop_family_orientation\ttop_family_consensus"
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
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{status}\t{tlen}\t{:.3}\t{}\t{top}",
            w.id,
            w.sequence_length,
            w.kind.name(),
            w.position,
            r.tandem_cover,
            r.arrays.len()
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
        let mut out = Vec::new();
        for c in 0..n {
            for (k, kind) in [Kind::Start, Kind::End, Kind::Interior]
                .into_iter()
                .enumerate()
            {
                let mut seq = if kind.is_end() { cap(c) } else { Vec::new() };
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
        let (_, f, s) = analyse(&w, Some("AAACCCT"), 2000);
        assert_eq!(s["class"], "telomeric");
        assert_eq!(s["top_family"], 0);
        assert_eq!(s["top_family_ends"], 12);
        assert_eq!(f[0].consensus, "AAACCCT");
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
