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

    // &str manipulation / rotation
    // see https://github.com/sebinsua/cracking-the-coding-interview/blob/9aed47bddaa0adbc527b650b8f28d0c751ce82e8/arrays-and-strings/src/string_rotation.rs
    // when there is an error/snp in the telomeric sequence, it causes a shift in the
    // repeat that is returned in 'explore' subcommand. Perhaps this info can be leveraged? Help.
    pub fn string_rotation(string1: &str, string2: &str) -> bool {
        if string1.len() != string2.len() {
            return false;
        }

        let first_character = string2.chars().nth(0).unwrap_or('\0');
        let rotation_idx = match string1.find(|chr| chr == first_character) {
            Some(idx) => idx,
            None => return false,
        };

        let rotated_string = format!("{}{}", &string1[rotation_idx..], &string1[..rotation_idx]);

        if string2.contains(&rotated_string) {
            return true;
        }

        false
    }
}
