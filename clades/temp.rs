use std::{
    boxed::Box,
    fmt::{self, Display},
};
use tabled::{Disable, Footer, MaxWidth, Modify, Rows, Table, Tabled};

#[derive(Debug, Clone)]
pub struct Seq<'a>(pub Box<&'a [&'a str]>);

impl Seq<'_> {
    pub fn get(&self, index: usize) -> Option<&&str> {
        self.0.get(index)
    }
}

impl Display for Seq<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let inner = self.0.join(", ");
        write!(f, "{}", inner)
    }
}
#[derive(Debug, Clone, Tabled)]
pub struct TelomereSeq<'a> {
    #[tabled(rename = "Clade")]
    pub clade: &'a str,
    #[tabled(rename = "Telomeric repeat units")]
    pub seq: Seq<'a>,
    pub length: usize,
}

// automated input start
// automated input end

pub fn print_table() {
    let mut clade_vec = Vec::new();

    for clade in CLADES {
        clade_vec.push(return_telomere_sequence(clade));
    }

    eprintln!(
        "{}",
        Table::new(&clade_vec)
            .with(
                Modify::new(Rows::new(1..clade_vec.len() - 1))
                    .with(MaxWidth::wrapping(30).keep_words()),
            )
            .with(Disable::Column(2..3))
            .with(Footer(
                "This table is modified from \"A telomeric repeat database\"\nhttps://github.com/tolkit/a-telomeric-repeat-database"
            ))
            .to_string()
    );
}
