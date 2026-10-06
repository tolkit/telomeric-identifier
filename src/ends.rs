use crate::explore::{self, Orientation};
use crate::{open_fasta_reader, utils};
use anyhow::{bail, Result};
use clap::crate_version;
use rayon::prelude::*;
use std::fs::{create_dir_all, File};
use std::io::{LineWriter, Write};
use std::path::PathBuf;

/// Runs of the unit closer together than this (bp) are merged into one telomere.
const MERGE_GAP: usize = 100;
// when discovering the unit, explore the end windows over these kmer
// lengths, counting runs of more than this many copies.
const DISCOVERY_MINIMUM: usize = 5;
const DISCOVERY_MAXIMUM: usize = 12;
const DISCOVERY_THRESHOLD: usize = 20;
// as in `tidk explore`, warn if the discovered unit is this strand biased
const STRAND_BIAS_WARNING: f64 = 0.1;

/// Which end of a sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Start,
    End,
}

impl Side {
    fn name(&self) -> &'static str {
        match self {
            Side::Start => "start",
            Side::End => "end",
        }
    }

    /// The strand of a telomere that reads along the sequence at this end.
    /// The G-rich strand runs 5' to 3' towards the chromosome end, so it
    /// reads along the sequence at the end, and the C-rich strand at the start.
    fn expected_strand(&self) -> Strand {
        match self {
            Side::Start => Strand::CRich,
            Side::End => Strand::GRich,
        }
    }
}

/// Which strand of the telomere reads along the sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Strand {
    GRich,
    CRich,
    /// The unit has as many Cs as Gs, so there is no way to tell.
    Unknown,
}

impl Strand {
    fn name(&self) -> &'static str {
        match self {
            Strand::GRich => "G-rich",
            Strand::CRich => "C-rich",
            Strand::Unknown => "NA",
        }
    }
}

/// The call for one end of a sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Status {
    Present,
    Absent,
    /// A telomere is present, but the wrong strand reads along the sequence.
    /// Often a misjoin or an inverted end, so it doesn't count towards T2T.
    WrongStrand,
    /// A telomere on the expected strand, but further than `--max-offset`
    /// from the end, so there is sequence beyond it to check. Doesn't count
    /// towards T2T.
    NotTerminal,
}

impl Status {
    fn name(&self) -> &'static str {
        match self {
            Status::Present => "present",
            Status::Absent => "absent",
            Status::WrongStrand => "wrong_strand",
            Status::NotTerminal => "not_terminal",
        }
    }
}

/// A telomere found at one end of a sequence, in sequence coordinates.
#[derive(Debug, PartialEq, Eq)]
struct Telomere {
    start: usize,
    end: usize,
    strand: Strand,
    /// The telomere reaches the inner edge of the search window, so it may
    /// be longer than reported.
    window_limited: bool,
}

#[derive(Debug, PartialEq, Eq)]
struct EndCall {
    status: Status,
    telomere: Option<Telomere>,
}

/// Both end calls for a sequence.
struct SequenceEnds {
    id: String,
    length: usize,
    start: EndCall,
    end: EndCall,
}

impl SequenceEnds {
    fn count(&self, status: Status) -> usize {
        [&self.start, &self.end]
            .iter()
            .filter(|c| c.status == status)
            .count()
    }
}

/// The ends of one sequence, kept in memory instead of the whole sequence.
struct EndWindows {
    id: String,
    length: usize,
    start: Vec<u8>,
    end: Vec<u8>,
}

/// Settings for calling telomeres.
struct EndsConfig {
    /// The primitive canonical repeat unit.
    unit: String,
    min_length: usize,
    max_offset: usize,
    error_tolerant: bool,
}

/// The entry point for `tidk ends`.
pub fn ends(matches: &clap::ArgMatches) -> Result<()> {
    let input_fasta = matches
        .get_one::<PathBuf>("fasta")
        .expect("errored by clap");
    let string = matches.get_one::<String>("string");
    let window = *matches
        .get_one::<usize>("window")
        .expect("defaulted by clap");
    let min_length = *matches
        .get_one::<usize>("min-length")
        .expect("defaulted by clap");
    let max_offset = *matches
        .get_one::<usize>("max-offset")
        .expect("defaulted by clap");
    let min_sequence_length = *matches
        .get_one::<usize>("min-sequence-length")
        .expect("defaulted by clap");
    let error_tolerant = !matches.get_flag("exact");
    let output = matches
        .get_one::<PathBuf>("output")
        .expect("errored by clap");
    let outdir = matches
        .get_one::<PathBuf>("dir")
        .expect("defaulted by clap");

    if min_length > window {
        bail!("--min-length ({min_length}) can't be longer than --window ({window}).");
    }

    // only keep the ends of each sequence in memory
    let mut records = Vec::new();
    for result in open_fasta_reader(input_fasta)?.records() {
        let record = result?;
        let seq = record.seq();
        if seq.len() < min_sequence_length {
            continue;
        }
        let w = window.min(seq.len());
        records.push(EndWindows {
            id: record.id().to_owned(),
            length: seq.len(),
            start: seq[..w].to_ascii_uppercase(),
            end: seq[seq.len() - w..].to_ascii_uppercase(),
        });
    }
    if records.is_empty() {
        bail!("No sequences of at least {min_sequence_length}bp in the input.");
    }

    let unit = match string {
        Some(s) => {
            let s = s.to_ascii_uppercase();
            if s.is_empty() || !s.bytes().all(|b| b"ACGT".contains(&b)) {
                bail!("--string must only contain A, C, G and T, but was {s}.");
            }
            utils::primitive_telomere_unit(&s)
        }
        None => discover_unit(&records, error_tolerant)?,
    };
    eprintln!(
        "[+]\tCalling telomeres with repeat unit {unit} at both ends of {} sequences",
        records.len()
    );

    let config = EndsConfig {
        unit: unit.clone(),
        min_length,
        max_offset,
        error_tolerant,
    };
    let calls: Vec<SequenceEnds> = records
        .par_iter()
        .map(|r| SequenceEnds {
            id: r.id.clone(),
            length: r.length,
            start: call_end(&r.start, Side::Start, 0, &config),
            end: call_end(&r.end, Side::End, r.length - r.end.len(), &config),
        })
        .collect();

    create_dir_all(outdir)?;
    let prefix = format!("{}/{}", outdir.display(), output.display());
    write_tsv(&calls, &format!("{prefix}.ends.tsv"))?;
    write_bed(&calls, &format!("{prefix}.ends.bed"))?;
    write_summary(&calls, &config, window, &format!("{prefix}.ends.json"))?;
    eprintln!("[+]\tWritten {prefix}.ends.tsv, .bed and .json");

    print_summary(&calls);

    Ok(())
}

/// Find the most abundant candidate telomeric repeat unit in the end windows.
fn discover_unit(records: &[EndWindows], error_tolerant: bool) -> Result<String> {
    eprintln!(
        "[+]\tNo --string given, discovering the telomeric repeat unit from the sequence ends"
    );
    let windows: Vec<(String, &[u8])> = records
        .iter()
        .flat_map(|r| {
            [
                (format!("{}:start", r.id), r.start.as_slice()),
                (format!("{}:end", r.id), r.end.as_slice()),
            ]
        })
        .collect();
    let Some((unit, minor)) = explore::top_unit(
        &windows,
        DISCOVERY_MINIMUM,
        DISCOVERY_MAXIMUM,
        DISCOVERY_THRESHOLD,
        error_tolerant,
    )?
    else {
        bail!(
            "No candidate telomeric repeat found at the sequence ends. Supply one with --string."
        );
    };
    eprintln!("[+]\tMost abundant repeat unit at the sequence ends: {unit}. Use --string if this is not the telomeric repeat.");
    if let Some(p) = minor {
        if p < STRAND_BIAS_WARNING {
            eprintln!(
                "[!]\t{unit}: only {:.1}% of copies are on the minor strand. Telomeres should read C-rich at sequence starts and G-rich at ends, so this may not be the telomeric repeat.",
                p * 100.0
            );
        }
    }
    Ok(unit)
}

/// The strand the unit reads as, in its canonical orientation.
fn unit_strand(unit: &str) -> Strand {
    let c = unit.bytes().filter(|&b| b == b'C').count();
    let g = unit.bytes().filter(|&b| b == b'G').count();
    match c.cmp(&g) {
        std::cmp::Ordering::Greater => Strand::CRich,
        std::cmp::Ordering::Less => Strand::GRich,
        std::cmp::Ordering::Equal => Strand::Unknown,
    }
}

/// The strand a run reads as, given its orientation relative to the unit.
fn run_strand(orientation: Orientation, unit_strand: Strand) -> Strand {
    match (orientation, unit_strand) {
        (Orientation::Either, _) | (_, Strand::Unknown) => Strand::Unknown,
        (Orientation::Unit, s) => s,
        (Orientation::Revcomp, Strand::CRich) => Strand::GRich,
        (Orientation::Revcomp, Strand::GRich) => Strand::CRich,
    }
}

/// Runs of the unit in a window, merged where they are less than
/// [`MERGE_GAP`] apart, as (start, end, strand) sorted by start. A block's
/// strand is the one covering most of it.
fn telomere_blocks(window: &[u8], config: &EndsConfig) -> Vec<(usize, usize, Strand)> {
    let unit_strand = unit_strand(&config.unit);
    let mut runs: Vec<(usize, usize, Strand)> =
        explore::find_runs(window, config.unit.len(), "", config.error_tolerant)
            .into_iter()
            .filter(|r| utils::primitive_telomere_unit(&r.sequence) == config.unit)
            .map(|r| {
                let orientation = explore::orientation(&r.sequence, &config.unit);
                (r.start, r.end, run_strand(orientation, unit_strand))
            })
            .collect();
    runs.sort_unstable_by_key(|r| r.0);

    let strands = [Strand::GRich, Strand::CRich, Strand::Unknown];
    let index = |s: Strand| strands.iter().position(|&x| x == s).unwrap();
    let majority = |bp: [usize; 3]| strands[(0..3).max_by_key(|&i| bp[i]).unwrap()];

    let mut blocks = Vec::new();
    // (start, end, bp covered by each strand)
    let mut current: Option<(usize, usize, [usize; 3])> = None;
    for (start, end, strand) in runs {
        match current.as_mut() {
            Some((_, ce, bp)) if start <= *ce + MERGE_GAP => {
                *ce = (*ce).max(end);
                bp[index(strand)] += end - start;
            }
            _ => {
                if let Some((cs, ce, bp)) = current {
                    blocks.push((cs, ce, majority(bp)));
                }
                let mut bp = [0; 3];
                bp[index(strand)] = end - start;
                current = Some((start, end, bp));
            }
        }
    }
    if let Some((cs, ce, bp)) = current {
        blocks.push((cs, ce, majority(bp)));
    }
    blocks
}

/// Call the telomere at one end of a sequence, from the window at that end.
/// `window_start` is the position of the window in the sequence.
fn call_end(window: &[u8], side: Side, window_start: usize, config: &EndsConfig) -> EndCall {
    let absent = EndCall {
        status: Status::Absent,
        telomere: None,
    };
    // distance of a block from the end of the sequence it should be at
    let offset = |&(start, end, _): &(usize, usize, Strand)| match side {
        Side::Start => start,
        Side::End => window.len() - end,
    };
    // runs are found in chunks of the unit length from the start of the
    // window, so at the end side, trim the window so chunks finish exactly
    // at the sequence end
    let trim = match side {
        Side::Start => 0,
        Side::End => window.len() % config.unit.len(),
    };
    let on_expected_strand =
        |strand: Strand| strand == Strand::Unknown || strand == side.expected_strand();
    let blocks: Vec<_> = telomere_blocks(&window[trim..], config)
        .into_iter()
        .map(|(start, end, strand)| (start + trim, end + trim, strand))
        .filter(|b| b.1 - b.0 >= config.min_length)
        .collect();
    // the longest block close enough to the end, else the closest block on
    // the expected strand further in. A block further in on the wrong strand
    // is more likely interstitial than a telomere, so is ignored.
    let terminal = blocks
        .iter()
        .filter(|b| offset(b) <= config.max_offset)
        .max_by_key(|b| b.1 - b.0);
    let (status, &(start, end, strand)) = match terminal {
        Some(b) if on_expected_strand(b.2) => (Status::Present, b),
        Some(b) => (Status::WrongStrand, b),
        None => match blocks
            .iter()
            .filter(|b| on_expected_strand(b.2))
            .min_by_key(|b| offset(b))
        {
            Some(b) => (Status::NotTerminal, b),
            None => return absent,
        },
    };

    let window_limited = match side {
        Side::Start => end + MERGE_GAP >= window.len(),
        Side::End => start <= MERGE_GAP,
    };
    EndCall {
        status,
        telomere: Some(Telomere {
            start: window_start + start,
            end: window_start + end,
            strand,
            window_limited,
        }),
    }
}

fn write_tsv(calls: &[SequenceEnds], path: &str) -> Result<()> {
    let mut f = LineWriter::new(File::create(path)?);
    writeln!(
        f,
        "id\tsequence_length\tside\tstatus\ttelomere_start\ttelomere_end\ttelomere_length\toffset\tstrand\twindow_limited"
    )?;
    for c in calls {
        for (side, call) in [(Side::Start, &c.start), (Side::End, &c.end)] {
            let fields = match &call.telomere {
                Some(t) => {
                    let offset = match side {
                        Side::Start => t.start,
                        Side::End => c.length - t.end,
                    };
                    format!(
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        t.start,
                        t.end,
                        t.end - t.start,
                        offset,
                        t.strand.name(),
                        t.window_limited
                    )
                }
                None => "NA\tNA\tNA\tNA\tNA\tNA".to_string(),
            };
            writeln!(
                f,
                "{}\t{}\t{}\t{}\t{fields}",
                c.id,
                c.length,
                side.name(),
                call.status.name()
            )?;
        }
    }
    Ok(())
}

/// Telomeres found, as BED (0-based, half open), named by end and status.
fn write_bed(calls: &[SequenceEnds], path: &str) -> Result<()> {
    let mut f = LineWriter::new(File::create(path)?);
    for c in calls {
        for (side, call) in [(Side::Start, &c.start), (Side::End, &c.end)] {
            if let Some(t) = &call.telomere {
                writeln!(
                    f,
                    "{}\t{}\t{}\t{}_{}",
                    c.id,
                    t.start,
                    t.end,
                    side.name(),
                    call.status.name()
                )?;
            }
        }
    }
    Ok(())
}

fn write_summary(
    calls: &[SequenceEnds],
    config: &EndsConfig,
    window: usize,
    path: &str,
) -> Result<()> {
    let with = |n: usize| {
        calls
            .iter()
            .filter(|c| c.count(Status::Present) == n)
            .count()
    };
    let flagged = |status: Status| {
        flagged_ends(calls, status)
            .map(|(c, side, offset)| serde_json::json!({"id": c.id, "side": side.name(), "offset": offset}))
            .collect::<Vec<_>>()
    };
    let summary = serde_json::json!({
        "tidk_version": crate_version!(),
        "unit": config.unit,
        "window": window,
        "min_length": config.min_length,
        "max_offset": config.max_offset,
        "error_tolerant": config.error_tolerant,
        "sequences": calls.len(),
        "t2t": with(2),
        "one_end": with(1),
        "no_ends": with(0),
        "t2t_ids": calls
            .iter()
            .filter(|c| c.count(Status::Present) == 2)
            .map(|c| &c.id)
            .collect::<Vec<_>>(),
        "wrong_strand": flagged(Status::WrongStrand),
        "not_terminal": flagged(Status::NotTerminal),
    });
    let mut f = File::create(path)?;
    serde_json::to_writer_pretty(&mut f, &summary)?;
    writeln!(f)?;
    Ok(())
}

fn print_summary(calls: &[SequenceEnds]) {
    let with = |n: usize| {
        calls
            .iter()
            .filter(|c| c.count(Status::Present) == n)
            .count()
    };
    eprintln!(
        "[+]\tT2T: {}/{} sequences. One telomere: {}. No telomeres: {}.",
        with(2),
        calls.len(),
        with(1),
        with(0)
    );
    for (c, side, _) in flagged_ends(calls, Status::WrongStrand) {
        eprintln!(
            "[!]\t{} {}: telomere on the wrong strand, possibly a misjoin or an inverted end.",
            c.id,
            side.name()
        );
    }
    for (c, side, offset) in flagged_ends(calls, Status::NotTerminal) {
        eprintln!(
            "[!]\t{} {}: telomere {offset}bp from the end, check the sequence beyond it.",
            c.id,
            side.name()
        );
    }
}

/// Ends with a given status, with the distance of the telomere from the end.
fn flagged_ends(
    calls: &[SequenceEnds],
    status: Status,
) -> impl Iterator<Item = (&SequenceEnds, Side, usize)> {
    calls.iter().flat_map(move |c| {
        [(Side::Start, &c.start), (Side::End, &c.end)]
            .into_iter()
            .filter(move |(_, call)| call.status == status)
            .filter_map(move |(side, call)| {
                let t = call.telomere.as_ref()?;
                let offset = match side {
                    Side::Start => t.start,
                    Side::End => c.length - t.end,
                };
                Some((c, side, offset))
            })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config(unit: &str) -> EndsConfig {
        EndsConfig {
            unit: unit.into(),
            min_length: 500,
            max_offset: 1000,
            error_tolerant: true,
        }
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

    #[test]
    fn test_unit_strand() {
        assert_eq!(unit_strand("AACCCT"), Strand::CRich);
        assert_eq!(unit_strand("TTAGGG"), Strand::GRich);
        assert_eq!(unit_strand("AATT"), Strand::Unknown);
    }

    #[test]
    fn test_start_telomere() {
        let seq = format!("{}{}", "CCCTAA".repeat(200), random_seq(5000, 1));
        let call = call_end(seq.as_bytes(), Side::Start, 0, &config("AACCCT"));
        assert_eq!(call.status, Status::Present);
        let t = call.telomere.unwrap();
        assert_eq!(
            (t.start, t.strand, t.window_limited),
            (0, Strand::CRich, false)
        );
        assert!(t.end >= 1194);
    }

    #[test]
    fn test_end_telomere_coordinates() {
        let window = format!("{}{}", random_seq(5000, 2), "TTAGGG".repeat(200));
        let call = call_end(window.as_bytes(), Side::End, 100_000, &config("AACCCT"));
        assert_eq!(call.status, Status::Present);
        let t = call.telomere.unwrap();
        assert_eq!(t.end, 106_200);
        assert_eq!(t.strand, Strand::GRich);
    }

    #[test]
    fn test_wrong_strand() {
        // G-rich at the start of a sequence
        let seq = format!("{}{}", "TTAGGG".repeat(200), random_seq(5000, 3));
        let call = call_end(seq.as_bytes(), Side::Start, 0, &config("AACCCT"));
        assert_eq!(call.status, Status::WrongStrand);
    }

    #[test]
    fn test_not_terminal_when_too_far_from_end() {
        let seq = format!(
            "{}{}{}",
            random_seq(2000, 4),
            "CCCTAA".repeat(200),
            random_seq(2000, 8)
        );
        let call = call_end(seq.as_bytes(), Side::Start, 0, &config("AACCCT"));
        assert_eq!(call.status, Status::NotTerminal);
        // runs are found in chunks, so the start is only known to within a unit
        let start = call.telomere.unwrap().start;
        assert!((1994..=2000).contains(&start), "{start}");
    }

    #[test]
    fn test_inner_block_on_wrong_strand_is_absent() {
        let seq = format!(
            "{}{}{}",
            random_seq(2000, 9),
            "TTAGGG".repeat(200),
            random_seq(2000, 10)
        );
        let call = call_end(seq.as_bytes(), Side::Start, 0, &config("AACCCT"));
        assert_eq!(call.status, Status::Absent);
    }

    #[test]
    fn test_terminal_preferred_over_longer_inner_block() {
        let seq = format!(
            "{}{}{}",
            "CCCTAA".repeat(100),
            random_seq(2000, 11),
            "CCCTAA".repeat(500)
        );
        let call = call_end(seq.as_bytes(), Side::Start, 0, &config("AACCCT"));
        assert_eq!(call.status, Status::Present);
        assert_eq!(call.telomere.unwrap().start, 0);
    }

    #[test]
    fn test_absent_when_too_short() {
        let seq = format!("{}{}", "CCCTAA".repeat(50), random_seq(5000, 5));
        let call = call_end(seq.as_bytes(), Side::Start, 0, &config("AACCCT"));
        assert_eq!(call.status, Status::Absent);
    }

    #[test]
    fn test_other_repeats_ignored() {
        let seq = format!("{}{}", "CCTAA".repeat(300), random_seq(5000, 6));
        let call = call_end(seq.as_bytes(), Side::Start, 0, &config("AACCCT"));
        assert_eq!(call.status, Status::Absent);
    }

    #[test]
    fn test_blocks_merge_across_small_gaps() {
        let seq = format!(
            "{}{}{}",
            "CCCTAA".repeat(100),
            random_seq(50, 7),
            "CCCTAA".repeat(100)
        );
        let blocks = telomere_blocks(seq.as_bytes(), &config("AACCCT"));
        assert_eq!(blocks.len(), 1);
        assert_eq!(blocks[0].2, Strand::CRich);
    }

    #[test]
    fn test_window_limited() {
        let seq = "CCCTAA".repeat(1000);
        let call = call_end(seq.as_bytes(), Side::Start, 0, &config("AACCCT"));
        assert!(call.telomere.unwrap().window_limited);
    }
}
