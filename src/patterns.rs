pub mod patterns {
    use std::boxed::Box;
    // all because there might be multiple telomere sequences...
    pub struct TelomereSeq<'a> {
        pub seq: Box<&'a [&'a str]>,
        pub length: usize,
    }
    // see http://telomerase.asu.edu/sequences_telomere.html
    pub fn return_telomere_sequence(clade: &str) -> TelomereSeq {
        let result = match clade {
            "vertebrate" => TelomereSeq {
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "ascidian" => TelomereSeq {
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "echinodermata" => TelomereSeq {
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "mollusca" => TelomereSeq {
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "coleoptera" => TelomereSeq {
                seq: Box::new(&["TCAGG", "TTAGG"]),
                length: 2,
            },
            "hymenoptera" => TelomereSeq {
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "lepidoptera" => TelomereSeq {
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "megaloptera" => TelomereSeq {
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "trichoptera" => TelomereSeq {
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "neuroptera" => TelomereSeq {
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "blattodea" => TelomereSeq {
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "orthoptera" => TelomereSeq {
                seq: Box::new(&["TTAGG"]),
                length: 1,
            },
            "nematoda" => TelomereSeq {
                seq: Box::new(&["TTAGGC", "TTGCA"]),
                length: 2,
            },
            "amoeba" => TelomereSeq {
                seq: Box::new(&["TTAGGG"]),
                length: 1,
            },
            "plants" => TelomereSeq {
                seq: Box::new(&["TTTAGGG", "TTAGGG", "TTTTAGGG"]),
                length: 3,
            },
            "ciliates" => TelomereSeq {
                seq: Box::new(&["TTGGGG", "TTTTGGGG"]),
                length: 2,
            },
            _ => panic!("{} is not yet accounted for in this pipeline.", clade),
        };
        result
    }

    // create all kmers of length k
    pub const BASES: &'static [char] = &['A', 'G', 'C', 'T'];
    // enumerate all kmers
    // this function will rapidly produce *millions* of values.
    // limit to 12?
    pub fn enumerate_mers(length: usize) -> Vec<String> {
        if length == 1 {
            BASES.iter().map(|i| i.to_string()).collect()
        } else if length > 12 {
            panic!(
                "Do not set length of kmer higher than 12 as this is  > 16 million comparisons per sequence."
            );
        } else {
            let strings = enumerate_mers(length - 1);
            BASES
                .iter()
                .flat_map(|c| strings.iter().map(move |rest| format!("{}{}", c, rest)))
                .collect()
        }
    }
}
