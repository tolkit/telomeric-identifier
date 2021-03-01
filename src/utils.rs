pub mod utils {
    use bio::pattern_matching::bom::BOM;
    // I'd like to use shift_and, but may have to wait until next public release
    // of rust-bio (or hard code it here...)
    // see https://github.com/rust-bio/rust-bio/blob/master/src/pattern_matching/shift_and.rs
    use bio::pattern_matching::kmp::KMP;

    #[derive(Debug)]
    pub struct Motifs {
        pub indexes: Vec<usize>,
        pub length: usize,
    }

    pub fn find_motifs(motif: &str, string: &str) -> Motifs {
        let motif_length = motif.len();
        let matches: Vec<usize>;

        if motif_length < 65 {
            let matcher = KMP::new(motif.as_bytes());
            matches = matcher.find_all(string.as_bytes()).collect::<Vec<usize>>();
        } else {
            let matcher = BOM::new(motif.as_bytes());
            matches = matcher.find_all(string.as_bytes()).collect::<Vec<usize>>();
        }

        Motifs {
            indexes: matches.to_owned(),
            length: matches.len(),
        }
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

    pub fn remove_overlapping_indexes(indexes: Motifs, pattern_length: usize) -> Vec<usize> {
        let mut indexes = indexes.indexes;
        let mut index = 0;
        let mut vec_len;

        loop {
            vec_len = indexes.len();

            if indexes.is_empty() || index == 0 || index == vec_len - 1 {
                break;
            }
            while indexes[index + 1] < indexes[index] + pattern_length {
                indexes.remove(index + 1);
                index += 1;
            }
        }
        indexes
    }

    // &str rotation
    // see https://github.com/rrbonham96/rust-ctci/blob/a2540532b098a06c29f2a5f06f54fc5717fd7669/src/arrays_and_strings/is_rotation.rs
    // when there is an error/snp in the telomeric sequence, it causes a shift in the
    // repeat that is returned in 'explore' subcommand. This info can be leveraged to count consecutive sequences...
    pub fn string_rotation(s1: &str, s2: &str) -> bool {
        if s1.len() == s2.len() {
            let mut s2 = s2.to_string();
            s2.push_str(&s2.clone());
            return s2.contains(&s1);
        }
        false
    }
}
