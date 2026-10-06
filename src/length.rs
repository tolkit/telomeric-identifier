use crate::ends::{self, EndsConfig, Side, Status, Strand, Telomere, STRAND_BIAS_WARNING};
use crate::open_fasta_reader;
use anyhow::Result;
use clap::crate_version;
use rayon::prelude::*;
use std::fs::{create_dir_all, File};
use std::io::{LineWriter, Write};
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::channel;

/// A telomere at one end of a read.
struct ReadTelomere {
    /// Position of the read in the input, to keep output in input order.
    index: usize,
    id: String,
    read_length: usize,
    side: Side,
    status: Status,
    telomere: Telomere,
}

impl ReadTelomere {
    fn length(&self) -> usize {
        self.telomere.end - self.telomere.start
    }

    /// There is non-telomeric sequence inward of the telomere, so the read
    /// spans the whole telomere, assuming the read end is the chromosome end.
    /// Otherwise the length is a lower bound.
    fn anchored(&self) -> bool {
        !self.telomere.window_limited
    }
}

/// The entry point for `tidk length`.
pub fn length(matches: &clap::ArgMatches) -> Result<()> {
    let input_fasta = matches
        .get_one::<PathBuf>("fasta")
        .expect("errored by clap");
    let string = matches
        .get_one::<String>("string")
        .expect("errored by clap");
    let min_length = *matches
        .get_one::<usize>("min-length")
        .expect("defaulted by clap");
    let max_offset = *matches
        .get_one::<usize>("max-offset")
        .expect("defaulted by clap");
    let error_tolerant = !matches.get_flag("exact");
    let output = matches
        .get_one::<PathBuf>("output")
        .expect("errored by clap");
    let outdir = matches
        .get_one::<PathBuf>("dir")
        .expect("defaulted by clap");

    let config = EndsConfig {
        unit: ends::parse_unit(string)?,
        min_length,
        max_offset,
        error_tolerant,
    };
    eprintln!(
        "[+]\tMeasuring telomeres with repeat unit {} at the ends of reads",
        config.unit
    );

    let scanned = AtomicUsize::new(0);
    let (sender, receiver) = channel();
    open_fasta_reader(input_fasta)?
        .records()
        .enumerate()
        .par_bridge()
        .for_each_with(sender, |s, (index, record)| {
            let record = record.expect("[-]\tError during fasta record parsing.");
            scanned.fetch_add(1, Ordering::Relaxed);
            for t in read_telomeres(index, record.id(), record.seq(), &config) {
                s.send(t).expect("Did not send!");
            }
        });
    let mut telomeres: Vec<ReadTelomere> = receiver.into_iter().collect();
    telomeres.sort_by_key(|t| (t.index, t.side == Side::End));
    let scanned = scanned.into_inner();

    create_dir_all(outdir)?;
    let prefix = format!("{}/{}", outdir.display(), output.display());
    write_tsv(&telomeres, &format!("{prefix}.length.tsv"))?;
    let groups = strand_groups(&telomeres);
    write_summary(
        &telomeres,
        &groups,
        scanned,
        &config,
        &format!("{prefix}.length.json"),
    )?;
    eprintln!("[+]\tWritten {prefix}.length.tsv and .json");
    print_summary(&telomeres, &groups, scanned);

    Ok(())
}

/// Telomeres at either end of a read. A telomere further into the read than
/// `--max-offset` can't be at a chromosome end, so only telomeres at the
/// read ends are kept (`present` or `wrong_strand`).
///
/// The length is that of the uninterrupted tract at the read end. If the read
/// has more telomeric repeat on the same strand further in (an interrupted
/// telomere, or a read lying within the telomere), the tract is not anchored
/// and its length is a lower bound. A same strand telomere at the other end
/// of such a read is part of the same telomere, so is not reported as well.
fn read_telomeres(index: usize, id: &str, seq: &[u8], config: &EndsConfig) -> Vec<ReadTelomere> {
    let seq = seq.to_ascii_uppercase();
    let blocks = ends::telomere_blocks(&seq, config);
    if blocks.is_empty() {
        return vec![];
    }
    let mut telomeres: Vec<ReadTelomere> = [Side::Start, Side::End]
        .into_iter()
        .filter_map(|side| {
            let call = ends::classify_end(&blocks, seq.len(), side, 0, config);
            match call.status {
                Status::Present | Status::WrongStrand => Some(ReadTelomere {
                    index,
                    id: id.to_string(),
                    read_length: seq.len(),
                    side,
                    status: call.status,
                    telomere: call.telomere?,
                }),
                Status::NotTerminal | Status::Absent => None,
            }
        })
        .collect();

    // both ends on the same strand: keep the end where that strand is expected
    if telomeres.len() == 2 && telomeres[0].telomere.strand == telomeres[1].telomere.strand {
        telomeres.retain(|t| t.status == Status::Present);
    }
    for t in &mut telomeres {
        let interrupted = blocks.iter().any(|&(start, end, strand)| {
            strand == t.telomere.strand
                && end - start >= config.min_length
                && (start, end) != (t.telomere.start, t.telomere.end)
        });
        if interrupted {
            t.telomere.window_limited = true;
        }
    }
    telomeres
}

fn write_tsv(telomeres: &[ReadTelomere], path: &str) -> Result<()> {
    let mut f = LineWriter::new(File::create(path)?);
    writeln!(
        f,
        "read_id\tread_length\tside\tstatus\tstrand\ttelomere_start\ttelomere_end\ttelomere_length\tanchored"
    )?;
    for t in telomeres {
        writeln!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            t.id,
            t.read_length,
            t.side.name(),
            t.status.name(),
            t.telomere.strand.name(),
            t.telomere.start,
            t.telomere.end,
            t.length(),
            t.anchored()
        )?;
    }
    Ok(())
}

/// Summary of telomere lengths from reads on one strand.
struct StrandGroup {
    strand: Strand,
    reads: usize,
    /// Lengths from anchored reads, sorted.
    anchored_lengths: Vec<usize>,
}

impl StrandGroup {
    /// Nearest rank percentile of the anchored lengths.
    fn percentile(&self, p: f64) -> Option<usize> {
        let n = self.anchored_lengths.len();
        if n == 0 {
            return None;
        }
        let rank = ((p / 100.0) * n as f64).ceil().max(1.0) as usize;
        Some(self.anchored_lengths[rank.min(n) - 1])
    }

    fn mean(&self) -> Option<f64> {
        let n = self.anchored_lengths.len();
        (n > 0).then(|| self.anchored_lengths.iter().sum::<usize>() as f64 / n as f64)
    }
}

/// `present` telomeres grouped by strand. Strands are summarised separately
/// as one may be basecalled much worse than the other (e.g. the G-rich strand
/// in older ONT data).
fn strand_groups(telomeres: &[ReadTelomere]) -> Vec<StrandGroup> {
    [Strand::GRich, Strand::CRich, Strand::Unknown]
        .into_iter()
        .map(|strand| {
            let present: Vec<_> = telomeres
                .iter()
                .filter(|t| t.status == Status::Present && t.telomere.strand == strand)
                .collect();
            let mut anchored_lengths: Vec<_> = present
                .iter()
                .filter(|t| t.anchored())
                .map(|t| t.length())
                .collect();
            anchored_lengths.sort_unstable();
            StrandGroup {
                strand,
                reads: present.len(),
                anchored_lengths,
            }
        })
        .filter(|g| g.reads > 0)
        .collect()
}

fn write_summary(
    telomeres: &[ReadTelomere],
    groups: &[StrandGroup],
    scanned: usize,
    config: &EndsConfig,
    path: &str,
) -> Result<()> {
    let strands: serde_json::Map<String, serde_json::Value> = groups
        .iter()
        .map(|g| {
            (
                g.strand.name().to_string(),
                serde_json::json!({
                    "reads": g.reads,
                    "anchored_reads": g.anchored_lengths.len(),
                    "median": g.percentile(50.0),
                    "mean": g.mean(),
                    "p10": g.percentile(10.0),
                    "p90": g.percentile(90.0),
                    "max": g.anchored_lengths.last(),
                }),
            )
        })
        .collect();
    let summary = serde_json::json!({
        "tidk_version": crate_version!(),
        "unit": config.unit,
        "min_length": config.min_length,
        "max_offset": config.max_offset,
        "error_tolerant": config.error_tolerant,
        "reads_scanned": scanned,
        "wrong_strand": telomeres.iter().filter(|t| t.status == Status::WrongStrand).count(),
        "strands": strands,
    });
    let mut f = File::create(path)?;
    serde_json::to_writer_pretty(&mut f, &summary)?;
    writeln!(f)?;
    Ok(())
}

fn print_summary(telomeres: &[ReadTelomere], groups: &[StrandGroup], scanned: usize) {
    let fmt = |x: Option<usize>| x.map_or("NA".to_string(), |x| x.to_string());
    for g in groups {
        eprintln!(
            "[+]\t{} telomeres: {} reads ({} anchored), median length {}bp (10th-90th percentile {}-{}bp)",
            g.strand.name(),
            g.reads,
            g.anchored_lengths.len(),
            fmt(g.percentile(50.0)),
            fmt(g.percentile(10.0)),
            fmt(g.percentile(90.0))
        );
    }
    if groups.is_empty() {
        eprintln!("[-]\tNo telomeres found at read ends in {scanned} reads.");
    }
    let wrong = telomeres
        .iter()
        .filter(|t| t.status == Status::WrongStrand)
        .count();
    if wrong > 0 {
        eprintln!("[!]\t{wrong} telomeres at read ends on the wrong strand (chimeric reads, interstitial repeats or basecalling artefacts), not included in lengths.");
    }
    // a G-rich and a C-rich group, with one much smaller than the other
    let reads: Vec<_> = groups
        .iter()
        .filter(|g| g.strand != Strand::Unknown)
        .map(|g| g.reads)
        .collect();
    let total: usize = reads.iter().sum();
    let minor = if reads.len() == 2 {
        reads[0].min(reads[1])
    } else {
        0
    };
    if total > 0 && (minor as f64 / total as f64) < STRAND_BIAS_WARNING {
        eprintln!("[!]\tOnly {:.1}% of telomeric reads are on the minor strand. One strand may be basecalled badly; compare the strands before trusting the lengths.", minor as f64 / total as f64 * 100.0);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn config() -> EndsConfig {
        EndsConfig {
            unit: "AACCCT".into(),
            min_length: 200,
            max_offset: 200,
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
    fn test_read_from_chromosome_start() {
        let read = format!("{}{}", "CCCTAA".repeat(500), random_seq(5000, 1));
        let t = read_telomeres(0, "r", read.as_bytes(), &config());
        assert_eq!(t.len(), 1);
        assert_eq!((t[0].side, t[0].status), (Side::Start, Status::Present));
        assert_eq!(t[0].telomere.strand, Strand::CRich);
        assert!(t[0].anchored());
        assert!((2994..=3000).contains(&t[0].length()));
    }

    #[test]
    fn test_read_from_chromosome_end() {
        let read = format!("{}{}", random_seq(5000, 2), "TTAGGG".repeat(500));
        let t = read_telomeres(0, "r", read.as_bytes(), &config());
        assert_eq!(t.len(), 1);
        assert_eq!((t[0].side, t[0].status), (Side::End, Status::Present));
        assert_eq!(t[0].telomere.strand, Strand::GRich);
    }

    #[test]
    fn test_unanchored_read() {
        let read = "TTAGGG".repeat(1000);
        let t = read_telomeres(0, "r", read.as_bytes(), &config());
        assert!(t.iter().all(|t| !t.anchored()));
    }

    #[test]
    fn test_wrong_strand_read() {
        let read = format!("{}{}", "TTAGGG".repeat(500), random_seq(5000, 3));
        let t = read_telomeres(0, "r", read.as_bytes(), &config());
        assert_eq!(t.len(), 1);
        assert_eq!(t[0].status, Status::WrongStrand);
    }

    #[test]
    fn test_interstitial_ignored() {
        let read = format!(
            "{}{}{}",
            random_seq(3000, 4),
            "CCCTAA".repeat(500),
            random_seq(3000, 5)
        );
        assert!(read_telomeres(0, "r", read.as_bytes(), &config()).is_empty());
    }

    #[test]
    fn test_percentile() {
        let g = StrandGroup {
            strand: Strand::GRich,
            reads: 5,
            anchored_lengths: vec![100, 200, 300, 400, 500],
        };
        assert_eq!(g.percentile(50.0), Some(300));
        assert_eq!(g.percentile(10.0), Some(100));
        assert_eq!(g.percentile(90.0), Some(500));
        assert_eq!(g.mean(), Some(300.0));
    }

    #[test]
    fn test_read_within_interrupted_telomere() {
        // C-rich at both ends, interrupted in the middle
        let read = format!(
            "{}{}{}",
            "CCCTAA".repeat(300),
            random_seq(3000, 6),
            "CCCTAA".repeat(300)
        );
        let t = read_telomeres(0, "r", read.as_bytes(), &config());
        assert_eq!(t.len(), 1);
        assert_eq!((t[0].side, t[0].status), (Side::Start, Status::Present));
        assert!(!t[0].anchored());
    }

    #[test]
    fn test_telomere_continuing_after_a_gap_is_not_anchored() {
        let read = format!(
            "{}{}{}{}",
            "CCCTAA".repeat(300),
            random_seq(2000, 7),
            "CCCTAA".repeat(300),
            random_seq(5000, 8)
        );
        let t = read_telomeres(0, "r", read.as_bytes(), &config());
        assert_eq!(t.len(), 1);
        assert!(!t[0].anchored());
    }

    #[test]
    fn test_entirely_telomeric_read_counted_once() {
        let read = "CCCTAA".repeat(1000);
        let t = read_telomeres(0, "r", read.as_bytes(), &config());
        assert_eq!(t.len(), 1);
        assert_eq!(t[0].status, Status::Present);
    }
}
