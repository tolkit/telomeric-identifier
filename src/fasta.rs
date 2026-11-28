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
//!   same point. bgzip (BGZF) files are inflated block by block in parallel
//!   instead, each block checked against its CRC and size; at anything else
//!   (another kind of member, corruption, a truncated file) reading goes
//!   back to `MultiGzDecoder` from the start, skipping what was already
//!   sent, so it fails as before.

use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::thread;

/// Bytes sent from the reading thread at a time.
const CHUNK: usize = 1 << 20;
/// Size of each read from a gzip decoder: the capacity of the `BufReader`
/// rust-bio read through. A decoder that hits corrupt data drops what it
/// decoded in that read, so this sets how much is read before the error.
const GZ_READ: usize = 8 * 1024;
/// Chunks in flight between the reading thread and the parser.
const IN_FLIGHT: usize = 4;
/// For bgzip files: compressed bytes inflated at a time, and chunks in
/// flight (about two batches, so the next is inflated while one is parsed).
const BGZF_BATCH: usize = 2 << 20;
const BGZF_IN_FLIGHT: usize = 16;
/// Most threads inflating bgzip blocks.
const BGZF_THREADS: usize = 8;

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

/// Is this the start of a BGZF block (a gzip member with only the BC extra
/// field)? Returns the block length.
fn bgzf_block_len(h: &[u8]) -> Option<usize> {
    (h.len() >= 18 && h[..4] == [0x1f, 0x8b, 8, 4] && h[10..16] == [6, 0, b'B', b'C', 2, 0])
        .then(|| u16::from_le_bytes([h[16], h[17]]) as usize + 1)
}

/// Inflate one BGZF block, if it is whole and its CRC and size match.
fn inflate_block(block: &[u8]) -> Option<Vec<u8>> {
    let n = block.len();
    if n < 26 {
        return None;
    }
    let cdata = &block[18..n - 8];
    let crc = u32::from_le_bytes(block[n - 8..n - 4].try_into().unwrap());
    let size = u32::from_le_bytes(block[n - 4..].try_into().unwrap()) as usize;
    if size > 1 << 16 {
        return None;
    }
    let mut out = Vec::with_capacity(size + 1);
    let mut d = flate2::Decompress::new(false);
    match d.decompress_vec(cdata, &mut out, flate2::FlushDecompress::Finish) {
        Ok(flate2::Status::StreamEnd)
            if d.total_in() as usize == cdata.len() && out.len() == size => {}
        _ => return None,
    }
    let mut c = flate2::Crc::new();
    c.update(&out);
    (c.sum() == crc).then_some(out)
}

/// What the reading side gives the parser: input, the error that stopped
/// it, or the end of the input.
enum Msg {
    Data(Vec<u8>),
    Fail(io::Error),
    End,
}

/// How reading a bgzip file in parallel ended.
enum Bgzf {
    /// All of it was sent.
    Done,
    /// Something else came up after this many bytes were sent.
    Fallback(usize),
    /// The parser has gone.
    Disconnected,
}

/// Send the inflated contents of a bgzip file in chunks, inflating blocks
/// in parallel, as long as every block is whole and checks out.
fn read_bgzf(file: &mut File, tx: &SyncSender<Msg>, threads: usize) -> Bgzf {
    let mut buf: Vec<u8> = Vec::new();
    let mut eof = false;
    let mut sent = 0;
    let mut chunk: Vec<u8> = Vec::with_capacity(CHUNK);
    let mut first = true;
    loop {
        while !eof && buf.len() < BGZF_BATCH {
            let start = buf.len();
            buf.resize(start + BGZF_BATCH, 0);
            match file.read(&mut buf[start..]) {
                Ok(0) => {
                    buf.truncate(start);
                    eof = true;
                }
                Ok(n) => buf.truncate(start + n),
                Err(e) if e.kind() == io::ErrorKind::Interrupted => buf.truncate(start),
                Err(_) => return Bgzf::Fallback(sent),
            }
        }
        // whole blocks in the buffer
        let mut blocks = Vec::new();
        let mut pos = 0;
        while pos < buf.len() {
            match bgzf_block_len(&buf[pos..]) {
                Some(len) if pos + len <= buf.len() => {
                    blocks.push(pos..pos + len);
                    pos += len;
                }
                // a block continuing past the buffer: read more
                Some(_) if !eof => break,
                None if !eof && buf.len() - pos < 18 => break,
                _ => {
                    // not a whole BGZF block: leave it to MultiGzDecoder
                    pos = usize::MAX;
                    break;
                }
            }
        }
        if first && blocks.is_empty() {
            // an empty file, or not bgzip after all
            return Bgzf::Fallback(0);
        }
        first = false;
        let stop = pos == usize::MAX;
        // inflate, a contiguous share of blocks per thread (in this thread
        // if one can't be started)
        let per = blocks.len().div_ceil(threads).max(1);
        let buf_ref = &buf;
        let work = move |group: &[std::ops::Range<usize>]| -> Vec<Option<Vec<u8>>> {
            group
                .iter()
                .map(|r| inflate_block(&buf_ref[r.clone()]))
                .collect()
        };
        let inflated: Vec<Option<Vec<u8>>> = thread::scope(|s| {
            let parts: Vec<_> = blocks
                .chunks(per)
                .map(|group| {
                    thread::Builder::new()
                        .name("fasta-inflate".into())
                        .spawn_scoped(s, move || work(group))
                        .map_err(|_| work(group))
                })
                .collect();
            parts
                .into_iter()
                .flat_map(|p| match p {
                    Ok(h) => h.join().expect("a bgzip inflate thread panicked"),
                    Err(done) => done,
                })
                .collect()
        });
        for out in inflated {
            let Some(out) = out else {
                if !chunk.is_empty() && tx.send(Msg::Data(std::mem::take(&mut chunk))).is_err() {
                    return Bgzf::Disconnected;
                }
                return Bgzf::Fallback(sent);
            };
            sent += out.len();
            chunk.extend_from_slice(&out);
            if chunk.len() >= CHUNK {
                let full = std::mem::replace(&mut chunk, Vec::with_capacity(CHUNK));
                if tx.send(Msg::Data(full)).is_err() {
                    return Bgzf::Disconnected;
                }
            }
        }
        if !chunk.is_empty() && tx.send(Msg::Data(std::mem::take(&mut chunk))).is_err() {
            return Bgzf::Disconnected;
        }
        if stop {
            return Bgzf::Fallback(sent);
        }
        buf.drain(..pos);
        if eof && buf.is_empty() {
            return Bgzf::Done;
        }
    }
}

/// Reads input sequentially, in chunks: a plain file, or gzip through
/// `MultiGzDecoder` in 8 KiB reads, leaving out the first `skip` bytes.
struct Producer {
    source: Box<dyn Read + Send>,
    step: usize,
    skip: usize,
    /// What comes after the chunk last given (the end, or an error).
    pending: Option<Msg>,
}

impl Producer {
    fn new(file: File, gz: bool, skip: usize) -> Self {
        let (source, step): (Box<dyn Read + Send>, usize) = if gz {
            (Box::new(MultiGzDecoder::new(file)), GZ_READ)
        } else {
            (Box::new(file), CHUNK)
        };
        Producer {
            source,
            step,
            skip,
            pending: None,
        }
    }

    /// The next chunk, then the end or the error that stopped the input.
    fn next(&mut self) -> Msg {
        loop {
            if let Some(m) = self.pending.take() {
                return m;
            }
            let mut chunk = vec![0u8; CHUNK];
            let mut filled = 0;
            let mut stop = None;
            while filled + self.step <= CHUNK {
                match self.source.read(&mut chunk[filled..filled + self.step]) {
                    Ok(0) => {
                        stop = Some(Msg::End);
                        break;
                    }
                    Ok(n) => filled += n,
                    Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                    Err(e) => {
                        stop = Some(Msg::Fail(e));
                        break;
                    }
                }
            }
            chunk.truncate(filled);
            if self.skip >= filled {
                self.skip -= filled;
                chunk.clear();
            } else if self.skip > 0 {
                chunk.drain(..self.skip);
                self.skip = 0;
            }
            if !chunk.is_empty() {
                self.pending = stop;
                return Msg::Data(chunk);
            }
            if let Some(m) = stop {
                return m;
            }
        }
    }
}

/// Where the parser's input comes from.
enum Source {
    /// A reading thread.
    Thread(Receiver<Msg>),
    /// This thread, if a reading thread couldn't be started.
    Inline(Producer),
}

impl Source {
    fn next(&mut self) -> Msg {
        match self {
            Source::Thread(rx) => rx.recv().unwrap_or_else(|_| {
                // the thread ended without saying so (it panicked)
                Msg::Fail(io::Error::other(
                    "the FASTA reading thread stopped before the end of the input",
                ))
            }),
            Source::Inline(p) => p.next(),
        }
    }
}

/// Open a file (gzipped if it ends in `.gz`) to be read on its own thread,
/// bgzip blocks inflated in parallel.
fn open_source(path: &Path) -> io::Result<Source> {
    let mut file = File::open(path)?;
    let gz = path.extension().and_then(|e| e.to_str()) == Some("gz");
    // a bgzip file? Only regular files are looked at (and rewound); any
    // read error here is met again when reading.
    let regular = file.metadata().is_ok_and(|m| m.is_file());
    let mut bgzf = false;
    if gz && regular {
        let mut header = [0u8; 18];
        bgzf = file.read_exact(&mut header).is_ok() && bgzf_block_len(&header).is_some();
        file.seek(SeekFrom::Start(0))?;
    }
    // inflate threads: no more than the thread pool has (`caps -t`)
    let threads = rayon::current_num_threads().clamp(1, BGZF_THREADS);
    let (tx, rx) = sync_channel(if bgzf { BGZF_IN_FLIGHT } else { IN_FLIGHT });
    let spawned = thread::Builder::new()
        .name("fasta-reader".into())
        .spawn(move || {
            // bytes already sent from the bgzip blocks, not to send again
            let mut skip = 0;
            if bgzf {
                match read_bgzf(&mut file, &tx, threads) {
                    Bgzf::Done => {
                        let _ = tx.send(Msg::End);
                        return;
                    }
                    Bgzf::Disconnected => return,
                    Bgzf::Fallback(sent) => {
                        if let Err(e) = file.seek(SeekFrom::Start(0)) {
                            let _ = tx.send(Msg::Fail(e));
                            return;
                        }
                        skip = sent;
                    }
                }
            }
            let mut producer = Producer::new(file, gz, skip);
            loop {
                let m = producer.next();
                let last = !matches!(m, Msg::Data(_));
                if tx.send(m).is_err() || last {
                    return;
                }
            }
        });
    match spawned {
        Ok(_) => Ok(Source::Thread(rx)),
        // no thread: read here instead, from the start
        Err(_) => Ok(Source::Inline(Producer::new(File::open(path)?, gz, 0))),
    }
}

/// A FASTA reader.
pub struct Reader {
    source: Source,
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
            source: open_source(path.as_ref())?,
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
            match self.source.next() {
                Msg::Data(chunk) => {
                    if self.data.is_empty() {
                        self.data = chunk;
                    } else {
                        self.data.extend_from_slice(&chunk);
                    }
                }
                Msg::Fail(e) => {
                    self.done = true;
                    return Err(e);
                }
                Msg::End => self.done = true,
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
            return Err(io::Error::other("Expected > at record start."));
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
