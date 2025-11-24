use anyhow::Result;
use std::fmt::{self, Display};
use tabled::{
    settings::{
        object::{Columns, Rows},
        Disable, Modify, Panel, Width,
    },
    Table, Tabled,
};

use crate::build::{get_database_path, TelomereRepeatRow};

/// A telomeric repeat sequence, or sequences.
#[derive(Debug, Clone)]
pub struct Seq(pub Vec<String>);

impl Seq {
    /// Create a new empty sequence.
    pub fn new() -> Self {
        Self(vec![])
    }
    /// Push a new sequence to the current list.
    pub fn push(&mut self, seq: String) {
        self.0.push(seq);
    }
    /// Get the sequence corresponding to an index.
    pub fn get(&self, index: usize) -> Option<&String> {
        self.0.get(index)
    }

    /// Get the inner list of sequences.
    pub fn get_inner(&self) -> &Vec<String> {
        &self.0
    }
}

impl Default for Seq {
    fn default() -> Self {
        Self::new()
    }
}

impl Display for Seq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let inner = self.0.join(", ");
        write!(f, "{inner}")
    }
}

/// All the relevant information about a
/// telomeric repeat sequence.
#[derive(Debug, Clone, Tabled)]
pub struct TelomereSeq {
    #[tabled(rename = "Clade")]
    /// The clade a telomeric repeat belongs to.
    pub clade: String,
    #[tabled(rename = "Telomeric repeat units")]
    /// The actual telomeric repeat sequence(s).
    pub seq: Seq,
    /// How many different telomeric repeats counted
    /// for a clade.
    pub length: usize,
}

impl TelomereSeq {
    /// Create a new telomeric sequence.
    pub fn new(clade: String, seq: Seq) -> Self {
        Self {
            clade,
            seq,
            length: 0,
        }
    }

    /// Push a new sequence to the list of sequences.
    /// But only if the sequence is not already present.
    pub fn push(&mut self, seq: String) {
        if !self.seq.0.contains(&seq) {
            self.seq.push(seq);
        }
    }

    /// Set the clade
    pub fn set_clade(&mut self, clade: String) {
        self.clade = clade;
    }

    /// Set the length from the number of sequences.
    pub fn set_length(&mut self) {
        self.length = self.seq.0.len();
    }
}

/// Read from a csv file containing all the clades
/// and only return a list of clades.
pub fn get_clades() -> Result<Vec<String>> {
    // open from disk
    let path = get_database_path()?;
    let mut rdr = csv::Reader::from_path(path)?;

    let mut out = vec![];

    for result in rdr.deserialize() {
        let record: TelomereRepeatRow = result?;
        // just the orders
        let order = record.order;
        out.push(order);
    }

    // remove duplicates and empty strings
    out.dedup();
    out.retain(|e| !e.is_empty());

    Ok(out)
}

/// A function to get a telomeric repeat sequence
/// given a clade name.
pub fn return_telomere_sequence(clade: String) -> Result<TelomereSeq> {
    let path = get_database_path()?;
    // read the csv file
    let mut rdr = csv::Reader::from_path(path)?;

    // iterate over records, if they match the clade
    // push all the sequences into a TelomereSeq object

    let mut telomere_seq = TelomereSeq::new(clade.clone(), Seq::new());

    for result in rdr.deserialize() {
        let record: TelomereRepeatRow = result?;
        if record.order == clade {
            telomere_seq.push(record.telomeric_repeat);
        }
    }

    // set the length
    telomere_seq.set_length();

    // return the telomeric sequence
    Ok(telomere_seq)
}
// automated input end

/// Pretty print a table containing all the information about
/// telomeric repeats that we currently have.
pub fn print_table() -> Result<()> {
    let mut clade_vec = Vec::new();

    let clades = get_clades()?;
    for clade in clades {
        clade_vec.push(return_telomere_sequence(clade)?);
    }

    eprintln!(
        "{}",
        Table::new(&clade_vec)
            .with(
                Modify::new(Rows::new(1..clade_vec.len() - 1)).with(Width::wrap(30).keep_words(true))
            )
            .with(Disable::column(Columns::new(2..3)))
            .with(Panel::footer(
                "This table is created from a curated database of repeats. This database can be found in its raw form here: https://github.com/tolkit/telomeric-identifier/tree/main/clades/curated.csv"
            )).with(Width::wrap(60).keep_words(true))
    );

    Ok(())
}
