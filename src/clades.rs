use std::{
    fmt::{self, Display},
    str::FromStr,
};

use anyhow::Context;
use strum::IntoEnumIterator;
use strum_macros::{Display, EnumIter, EnumString, EnumVariantNames};

use tabled::{
    object::{Columns, Rows},
    Disable, Modify, Panel, Table, Tabled, Width,
};

/// A telomeric repeat sequence, or sequences.
#[derive(Debug, Clone)]
pub struct Seq<'a> {
    pub inner: &'a [&'a str],
}

impl Seq<'_> {
    /// Get the sequence corresponding to an index.
    pub fn get(&self, index: usize) -> Option<&&str> {
        self.inner.get(index)
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }
}

impl Display for Seq<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let inner = self.inner.join(", ");
        write!(f, "{}", inner)
    }
}

/// All the relevant information about a
/// telomeric repeat sequence.
#[derive(Debug, Clone, Tabled)]
pub struct TelomereSeq<'a> {
    #[tabled(rename = "Clade")]
    /// The clade a telomeric repeat belongs to.
    pub clade: Clade,
    #[tabled(rename = "Telomeric repeat units")]
    /// The actual telomeric repeat sequence(s).
    pub seq: Seq<'a>,
}

// automated input start

#[derive(Display, EnumString, EnumVariantNames, EnumIter, Debug, Clone)]
/// All the clades for which we have data.
pub enum Clade {
    Accipitriformes,
    Actiniaria,
    Anura,
    Apiales,
    Aplousobranchia,
    Asterales,
    Buxales,
    Caprimulgiformes,
    Carangiformes,
    Carcharhiniformes,
    Cardiida,
    Carnivora,
    Caryophyllales,
    Cheilostomatida,
    Chiroptera,
    Chlamydomonadales,
    Coleoptera,
    Crassiclitellata,
    Cypriniformes,
    Eucoccidiorida,
    Fabales,
    Fagales,
    Forcipulatida,
    Hemiptera,
    Heteronemertea,
    Hirudinida,
    Hymenoptera,
    Hypnales,
    Labriformes,
    Lamiales,
    Lepidoptera,
    Malpighiales,
    Myrtales,
    Odonata,
    Orthoptera,
    Pectinida,
    Perciformes,
    Phlebobranchia,
    Phyllodocida,
    Plecoptera,
    Pleuronectiformes,
    Poales,
    Rodentia,
    Rosales,
    Salmoniformes,
    Sapindales,
    Solanales,
    Symphypleona,
    Syngnathiformes,
    Trichoptera,
    Trochida,
    Venerida,
}

impl From<Clade> for Seq<'_> {
    fn from(value: Clade) -> Self {
        use Clade::*;
        let seq: &[_] = match value {
            Accipitriformes | Actiniaria | Anura | Aplousobranchia | Caprimulgiformes
            | Syngnathiformes | Trochida | Venerida | Carangiformes | Carcharhiniformes
            | Cardiida | Carnivora | Labriformes | Pectinida | Perciformes | Phlebobranchia
            | Phyllodocida | Pleuronectiformes | Rodentia | Salmoniformes => &["AACCCT"],

            Apiales | Asterales | Buxales | Caryophyllales | Eucoccidiorida | Fabales | Fagales
            | Hypnales | Lamiales | Malpighiales | Myrtales | Poales | Rosales | Sapindales => {
                &["AAACCCT"]
            }

            Cheilostomatida => &["AAACCCC"],

            Chiroptera | Chlamydomonadales | Cypriniformes | Forcipulatida | Heteronemertea
            | Hirudinida => &["AACCCT"],

            Coleoptera => &["AACCT", "ACCTG", "AACAGACCCG", "AACCC"],

            Crassiclitellata => &["AAGGAC"],

            Hemiptera => &["AAACCACCCT", "AACCATCCCT"],

            Hymenoptera => &[
                "AAACCC",
                "AAAGAACCT",
                "AACCCAGACGC",
                "AACCCGAACCT",
                "AACCCTGACGC",
                "AAAATTGTCCGTCC",
                "AACCC",
                "AACCCCAACCT",
                "AAATGTGGAGG",
                "AACCCAGACCC",
                "ACCCAG",
                "AACCCAGACCT",
                "AACCCT",
                "ACGGCAGCG",
                "AACCT",
            ],

            Lepidoptera | Orthoptera | Plecoptera | Symphypleona | Trichoptera => &["AACCT"],

            Odonata => &["AACCC"],

            Solanales => &["AACCCTG"],
        };
        Seq { inner: seq }
    }
}

/// A function to get a telomeric repeat sequence
/// given a clade name.
impl From<Clade> for TelomereSeq<'_> {
    fn from(value: Clade) -> Self {
        Self {
            clade: value.clone(),
            seq: value.into(),
        }
    }
}

impl FromStr for TelomereSeq<'_> {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Clade::from_str(s)
            .map(Into::into)
            .with_context(|| format!("{} is not yet accounted for in this pipeline.", s))
    }
}

/// Pretty print a table containing all the information about
/// telomeric repeats that we currently have.
pub fn print_table() {
    let clade_vec: Vec<TelomereSeq> = Clade::iter().map(Into::into).collect();
    eprintln!(
        "{}",
        Table::new(&clade_vec)
            .with(
                Modify::new(Rows::new(1..clade_vec.len() - 1)).with(Width::wrap(30).keep_words())
            )
            .with(Disable::column(Columns::new(2..3)))
            .with(Panel::footer(
                "This table is created from a curated database of repeats. This database can be found in its raw form here: https://github.com/tolkit/telomeric-identifier/tree/main/clades/curated.csv"
            )).with(Width::wrap(60).keep_words())
    );
}
