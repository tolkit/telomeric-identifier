pub mod utils {
    use bio::pattern_matching::shift_and;
    pub struct Motifs {
        pub indexes: Vec<usize>,
        pub length: usize,
    }

    pub fn find_motifs(motif: &str, string: &str) -> Motifs {
        let matcher = shift_and::ShiftAnd::new(motif.as_bytes());
        let matches = &matcher.find_all(string.as_bytes()).collect::<Vec<usize>>();

        Motifs {
            indexes: matches.to_owned(),
            length: matches.len(),
        }
    }
    pub fn longest_repeat(indexes: Motifs) -> Vec<usize> {
        let mut res = Vec::new();
        for index in 1..indexes.length {
            if index + 1 < indexes.length {
                res.push(indexes.indexes[index + 1] - indexes.indexes[index])
            }
        }
        res
    }

    pub fn reverse_complement(dna: &str) -> String {
        let dna_chars = dna.chars();
        let mut revcomp = Vec::new();

        for base in dna_chars {
            revcomp.push(switch_base(base))
        }
        revcomp.as_mut_slice().reverse();
        revcomp.into_iter().collect()
    }

    fn switch_base(c: char) -> char {
        match c {
            'A' => 'T',
            'C' => 'G',
            'T' => 'A',
            'G' => 'C',
            'N' => 'N',
            _ => 'N',
        }
    }
}
