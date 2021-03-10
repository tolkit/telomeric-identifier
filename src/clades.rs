pub mod clades {
    use std::boxed::Box;

    pub const CLADES: &[&str] = &[
        "vertebrate",
        "ascidian",
        "echinodermata",
        "mollusca",
        "coleoptera",
        "hymenoptera",
        "lepidoptera",
        "megaloptera",
        "trichoptera",
        "neuroptera",
        "blattodea",
        "orthoptera",
        "nematoda",
        "amoeba",
        "plants",
        "ciliates",
    ];

    // all because there might be multiple telomere sequences...

    #[derive(Debug, Clone)]
    pub struct TelomereSeq<'a> {
        pub clade: &'a str,
        pub seq: Box<&'a [&'a str]>,
        pub length: usize,
    }
    // see http://telomerase.asu.edu/sequences_telomere.html

    pub fn return_telomere_sequence(clade: &str) -> TelomereSeq {
        let result = match clade {
            "vertebrate" => TelomereSeq {
                clade: "Vertebrate",
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "ascidian" => TelomereSeq {
                clade: "Ascidian",
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "echinodermata" => TelomereSeq {
                clade: "Echinodermata",
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "mollusca" => TelomereSeq {
                clade: "Mollusca",
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "coleoptera" => TelomereSeq {
                clade: "Coleoptera",
                seq: Box::new(&["TCAGG", "TTAGG"]),
                length: 2,
            },
            "hymenoptera" => TelomereSeq {
                clade: "Hymenoptera",
                // general, Bombus, Ectemnius, Vespula, Vespa
                seq: Box::new(&[
                    "TTAGG",
                    "TTAGGTTGGGG",
                    "TTAGGTCTGGG",
                    "TTGCGTCTGGG",
                    "TTGCGTCAGGG",
                ]),
                length: 5,
            },
            "lepidoptera" => TelomereSeq {
                clade: "Lepidoptera",
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "megaloptera" => TelomereSeq {
                clade: "Megaloptera",
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "trichoptera" => TelomereSeq {
                clade: "Trichoptera",
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "neuroptera" => TelomereSeq {
                clade: "Neuroptera",
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "blattodea" => TelomereSeq {
                clade: "Blattodea",
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "orthoptera" => TelomereSeq {
                clade: "Orthoptera",
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "nematoda" => TelomereSeq {
                clade: "Nematoda",
                seq: Box::new(&["TTAGGC", "TTGCA"]),
                length: 2,
            },
            "amoeba" => TelomereSeq {
                clade: "Amoeba",
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "plants" => TelomereSeq {
                clade: "Plants",
                seq: Box::new(&["TTTAGGG", "TTAGGG", "TTTTAGGG"]),
                length: 3,
            },
            "ciliates" => TelomereSeq {
                clade: "Ciliates",
                seq: Box::new(&["TTGGGG", "TTTTGGGG"]),
                length: 2,
            },
            _ => panic!("{} is not yet accounted for in this pipeline.", clade),
        };
        result
    }

    // 'pretty' print a simplified table from the telomerase sequence database.

    pub fn print_table() {
        let clades = CLADES;

        println!("\nTelomeric repeats by clade\nSee http://telomerase.asu.edu/sequences_telomere.html for a full list\n");
        println!("{:^14} | {:^30}", "Clade", "Telomeric repeat");

        for clade in clades {
            println!("{:-<50}", "");
            let clade_info = return_telomere_sequence(clade);
            println!("{:^14} | {:^30?}", clade_info.clade, clade_info.seq);
        }
        println!("\n");
    }
}
