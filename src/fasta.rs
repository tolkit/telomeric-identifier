//! A FASTA reader with the behaviour of `bio::io::fasta::Reader` (rust-bio
//! 2.0.3), which tidk used before, but faster: the input is read (and
//! decompressed) on its own thread, and lines are not copied twice.
//!
//! What it keeps from rust-bio, so the output of every subcommand is the
//! same, errors included:
//!
//! - a record starts at a line beginning with `>`; anything else where a
//!   record should start is the error "Expected > at record start.";
//! - the header, without `>` and trailing whitespace, splits at the first
//!   whitespace into the id and the description;
//! - sequence lines are added with their trailing whitespace (Unicode
//!   `char::is_whitespace`) removed, so blank lines add nothing and spaces
//!   inside a line are kept;
//! - a record with an empty id, no description and no sequence (a bare `>`)
//!   ends the records, as does the first error;
//! - every line must be UTF-8, or the error is the one
//!   `BufRead::read_line` gives ("stream did not contain valid UTF-8"); an
//!   I/O error while reading a line takes precedence over it;
//! - gzipped input (by the `.gz` extension) is read with `MultiGzDecoder` in
//!   reads of 8 KiB, as a `BufReader` did, so a corrupt stream fails at the
//!   same point.

use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, Read};
use std::path::Path;
use std::sync::mpsc::{sync_channel, Receiver};
use std::thread;

/// Bytes sent from the reading thread at a time.
const CHUNK: usize = 1 << 20;
/// Size of each read from a gzip decoder: the capacity of the `BufReader`
/// rust-bio read through. A decoder that hits corrupt data drops what it
/// decoded in that read, so this sets how much is read before the error.
const GZ_READ: usize = 8 * 1024;
/// Chunks in flight between the reading thread and the parser.
const IN_FLIGHT: usize = 4;

/// One FASTA record.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: Vec<u8>,
}

impl Record {
    /// A record from its parts.
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: &[u8]) -> Self {
        Record {
            id: id.to_owned(),
            desc: desc.map(|d| d.to_owned()),
            seq: seq.to_vec(),
        }
    }

    /// The id: the header up to the first whitespace.
    pub fn id(&self) -> &str {
        &self.id
    }

    /// The rest of the header, if any.
    pub fn desc(&self) -> Option<&str> {
        self.desc.as_deref()
    }

    /// The sequence, as in the file (case and all).
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    fn is_empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_none() && self.seq.is_empty()
    }
}

/// The error `BufRead::read_line` gives for a line that is not UTF-8,
/// taken from `read_line` itself so it is the same value.
fn invalid_utf8() -> io::Error {
    let mut s = String::new();
    (&b"\xff\n"[..])
        .read_line(&mut s)
        .expect_err("0xff is not UTF-8")
}

/// Is `b` whitespace by `char::is_whitespace`, for an ASCII byte?
fn ascii_whitespace(b: u8) -> bool {
    matches!(b, b'\t' | b'\n' | 0x0b | 0x0c | b'\r' | b' ')
}

/// Read a file (gzipped if it ends in `.gz`) on a thread, sending its bytes
/// in chunks, then the error that stopped it, if any.
fn spawn_reader(path: &Path) -> io::Result<Receiver<io::Result<Vec<u8>>>> {
    let file = File::open(path)?;
    let gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
    let (tx, rx) = sync_channel(IN_FLIGHT);
    thread::Builder::new()
        .name("fasta-reader".into())
        .spawn(move || {
            let (mut source, step): (Box<dyn Read>, usize) = if gz {
                (Box::new(MultiGzDecoder::new(file)), GZ_READ)
            } else {
                (Box::new(file), CHUNK)
            };
            loop {
                let mut chunk = vec![0u8; CHUNK];
                let mut filled = 0;
                let mut stop = None;
                while filled + step <= CHUNK {
                    match source.read(&mut chunk[filled..filled + step]) {
                        Ok(0) => {
                            stop = Some(None);
                            break;
                        }
                        Ok(n) => filled += n,
                        Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                        Err(e) => {
                            stop = Some(Some(e));
                            break;
                        }
                    }
                }
                chunk.truncate(filled);
                if filled > 0 && tx.send(Ok(chunk)).is_err() {
                    return;
                }
                match stop {
                    // end of input: dropping the sender tells the parser
                    Some(None) => return,
                    Some(Some(e)) => {
                        let _ = tx.send(Err(e));
                        return;
                    }
                    None => {}
                }
            }
        })?;
    Ok(rx)
}

/// A FASTA reader.
pub struct Reader {
    rx: Receiver<io::Result<Vec<u8>>>,
    /// Unparsed input is `data[pos..]`.
    data: Vec<u8>,
    pos: usize,
    /// The reading thread has finished (end of input or an error).
    done: bool,
    /// The line after the last record, a header (rust-bio's `line`).
    line: Vec<u8>,
}

/// A line: (start, end) in `data`, end after the newline if there is one.
type Span = (usize, usize);

impl Reader {
    /// Open a FASTA file, gzipped if its name ends in `.gz`.
    pub fn from_path<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        Ok(Reader {
            rx: spawn_reader(path.as_ref())?,
            data: Vec::new(),
            pos: 0,
            done: false,
            line: Vec::new(),
        })
    }

    /// The next line including its newline, empty at the end of the input,
    /// as `BufRead::read_line` reads it but without UTF-8 checks.
    fn read_line(&mut self) -> io::Result<Span> {
        let mut searched = self.pos;
        loop {
            if let Some(i) = memchr::memchr(b'\n', &self.data[searched..]) {
                let span = (self.pos, searched + i + 1);
                self.pos = span.1;
                return Ok(span);
            }
            searched = self.data.len();
            if self.done {
                let span = (self.pos, self.data.len());
                self.pos = span.1;
                return Ok(span);
            }
            // more input: drop what is parsed, then append the next chunk
            if self.pos > 0 {
                self.data.drain(..self.pos);
                searched -= self.pos;
                self.pos = 0;
            }
            match self.rx.recv() {
                Ok(Ok(chunk)) => {
                    if self.data.is_empty() {
                        self.data = chunk;
                    } else {
                        self.data.extend_from_slice(&chunk);
                    }
                }
                Ok(Err(e)) => {
                    self.done = true;
                    return Err(e);
                }
                Err(_) => self.done = true,
            }
        }
    }

    /// Read the next record into `record`, left empty at the end of the input.
    fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.id.clear();
        record.desc = None;
        record.seq.clear();
        if self.line.is_empty() {
            let (a, b) = self.read_line()?;
            if a == b {
                return Ok(());
            }
            self.line.extend_from_slice(&self.data[a..b]);
        }
        let header = std::str::from_utf8(&self.line).map_err(|_| invalid_utf8())?;
        if !header.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
            ));
        }
        let mut fields = header[1..].trim_end().splitn(2, char::is_whitespace);
        record.id = fields.next().map(|s| s.to_owned()).unwrap();
        record.desc = fields.next().map(|s| s.to_owned());
        self.line.clear();
        loop {
            let (a, b) = self.read_line()?;
            let line = &self.data[a..b];
            if line.is_empty() {
                break;
            }
            if line[0] == b'>' {
                std::str::from_utf8(line).map_err(|_| invalid_utf8())?;
                self.line.extend_from_slice(line);
                break;
            }
            if line.is_ascii() {
                let end = line
                    .iter()
                    .rposition(|&b| !ascii_whitespace(b))
                    .map_or(0, |i| i + 1);
                record.seq.extend_from_slice(&line[..end]);
            } else {
                let s = std::str::from_utf8(line).map_err(|_| invalid_utf8())?;
                record.seq.extend_from_slice(s.trim_end().as_bytes());
            }
        }
        Ok(())
    }

    /// The records, ending at the first error (which is returned) or at a
    /// record with no id, description or sequence.
    pub fn records(self) -> Records {
        Records {
            reader: self,
            error_has_occured: false,
        }
    }
}

/// An iterator over the records of a FASTA file.
pub struct Records {
    reader: Reader,
    error_has_occured: bool,
}

impl Iterator for Records {
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        if self.error_has_occured {
            return None;
        }
        let mut record = Record::default();
        match self.reader.read(&mut record) {
            Ok(()) if record.is_empty() => None,
            Ok(()) => Some(Ok(record)),
            Err(e) => {
                self.error_has_occured = true;
                Some(Err(e))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// Records read by this reader and by rust-bio from the same bytes.
    fn both(name: &str, bytes: &[u8]) -> (Vec<String>, Vec<String>) {
        let dir = std::env::temp_dir().join(format!("tidk-fasta-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::File::create(&path)
            .unwrap()
            .write_all(bytes)
            .unwrap();
        let fmt = |r: io::Result<(String, Option<String>, Vec<u8>)>| match r {
            Ok((id, desc, seq)) => format!("{id:?} {desc:?} {:?}", String::from_utf8_lossy(&seq)),
            Err(e) => format!("error {e:?}"),
        };
        let ours = Reader::from_path(&path)
            .unwrap()
            .records()
            .map(|r| fmt(r.map(|r| (r.id, r.desc, r.seq))))
            .collect();
        let bio = crate::open_fasta_reader(&path)
            .unwrap()
            .records()
            .map(|r| {
                fmt(r.map(|r| {
                    (
                        r.id().to_owned(),
                        r.desc().map(|d| d.to_owned()),
                        r.seq().to_vec(),
                    )
                }))
            })
            .collect();
        std::fs::remove_file(&path).unwrap();
        (ours, bio)
    }

    #[test]
    fn test_same_as_rust_bio() {
        let cases: &[&[u8]] = &[
            b"",
            b"\n",
            b">a\nACGT\n",
            b">a desc here\nAC\nGT\n>b\n\nTT \t\n",
            b">a\nACGT",
            b">only",
            b">\nACGT\n>b\nAA\n",
            b">a\nAA\n>\n>b\nCC\n",
            b"> lead\nAA\n",
            b">a\tb c\n  AA  \n\x0bCC\x0c\r\n",
            b">a\r\nAC\r\n\r\nGT\r\n",
            b"ACGT\n>a\nAA\n",
            b">a\nA\xc2\xa0\nC\xc2\xa0C\n",
            b">a\nAC\xffGT\n>b\nAA\n",
            b">\xff\nAA\n",
            b">a\nAA\n>\xffb\nCC\n",
            b">a\n;c\nAA\n",
        ];
        for (i, c) in cases.iter().enumerate() {
            let (ours, bio) = both(&format!("case{i}.fa"), c);
            assert_eq!(ours, bio, "case {i}: {:?}", String::from_utf8_lossy(c));
        }
    }

    #[test]
    fn test_lines_across_chunks() {
        // a sequence line longer than a chunk, and lines straddling chunk ends
        let mut f = b">long\n".to_vec();
        f.extend(std::iter::repeat_n(b'A', CHUNK + 17));
        f.extend_from_slice(b"\n>short desc\n");
        for i in 0..(CHUNK / 7 + 3) {
            f.extend_from_slice(if i % 2 == 0 { b"ACGTAC\n" } else { b"TTAG  \n" });
        }
        let (ours, bio) = both("chunks.fa", &f);
        assert_eq!(ours, bio);
    }
}
