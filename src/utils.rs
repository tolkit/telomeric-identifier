pub mod utils {
    use bio::pattern_matching::bom::BOM;
    // I'd like to use shift_and, but may have to wait until next public release
    // of rust-bio (or hard code it here...)
    // see https://github.com/rust-bio/rust-bio/blob/master/src/pattern_matching/shift_and.rs
    use bio::pattern_matching::kmp::KMP;
    use lexical_sort::{natural_lexical_cmp, StringSort};
    use std::cmp::min;

    // this does the hard lifting in `tidk search` and `tidk find`
    // take input putative telomeric repeat (motif) and search against
    // a dna sequence. Optimised for motif length.

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

    // calculate the reverse complement of a telomeric repeat

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

    // not sure if this function is necessary, but it looks at the indexes of the motifs
    // in the genome and removes indexes which occur consecutively less than the pattern length apart

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

    // string rotation algorithms
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

    // Booth's algorithm for the lexicographically minimal
    // rotation of a string. Should give us a canonical rotation
    // given a rotated string.
    // written by https://github.com/zimpha/algorithmic-library/blob/61e897983314033615bcd278d22a754bfc3c3f22/rust/src/strings/mod.rs

    fn minimal_rotation<T: Ord>(s: &[T]) -> usize {
        let n = s.len();
        let mut i = 0;
        let mut j = 1;
        loop {
            let mut k = 0;
            let mut ci = &s[i % n];
            let mut cj = &s[j % n];
            while k < n {
                ci = &s[(i + k) % n];
                cj = &s[(j + k) % n];
                if ci != cj {
                    break;
                }
                k += 1
            }
            if k == n {
                return min(i, j);
            }
            if ci > cj {
                i += k + 1;
                i += (i == j) as usize;
            } else {
                j += k + 1;
                j += (i == j) as usize;
            }
        }
    }

    // a wrapper for `minimal_rotation` which gives us a string back
    // we have two sequences we *know* are string rotations of one another
    // or string rotations of the reverse complement
    // so find the lexographically minimal version of a string/its reverse complement.

    pub fn lms(telomeric_repeat1: &str, telomeric_repeat2: &str) -> String {
        // get index of where to rotate
        // for forward
        let index_f = minimal_rotation(telomeric_repeat1.as_bytes());
        let index_r = minimal_rotation(telomeric_repeat2.as_bytes());
        // for reverse
        let telomeric_repeat1_r = reverse_complement(telomeric_repeat1);
        let telomeric_repeat2_r = reverse_complement(telomeric_repeat2);
        // reverse indexes
        let index_fr = minimal_rotation(telomeric_repeat1_r.as_bytes());
        let index_rr = minimal_rotation(telomeric_repeat2_r.as_bytes());

        // put 0 -> index at end
        let end_f = &telomeric_repeat1[0..index_f];
        let end_r = &telomeric_repeat2[0..index_r];
        let end_fr = &telomeric_repeat1_r[0..index_fr];
        let end_rr = &telomeric_repeat2_r[0..index_rr];

        // put index -> end at start
        let start_f = &telomeric_repeat1[index_f..];
        let start_r = &telomeric_repeat2[index_r..];
        let start_fr = &telomeric_repeat1_r[index_fr..];
        let start_rr = &telomeric_repeat2_r[index_rr..];

        // give us the string
        let lms_f = format!("{}{}", start_f, end_f);
        let lms_r = format!("{}{}", start_r, end_r);
        let lms_fr = format!("{}{}", start_fr, end_fr);
        let lms_rr = format!("{}{}", start_rr, end_rr);

        // now we have four strings, and have to report one
        let mut strings = vec![&lms_f, &lms_r, &lms_fr, &lms_rr];
        strings.string_sort_unstable(natural_lexical_cmp);
        strings[0].to_string()
    }

    // given a string (of DNA)
    // return the lexicographical minimal representation
    // accounting for both forward & rev comp.
    // I guess it will kind of be similar to the above function
    // except we will expose this to the API

    pub fn lex_min(dna_string: &str) -> String {
        // revcomp
        let dna_string_r = reverse_complement(dna_string);
        // index for forward and reverse
        let index_f = minimal_rotation(dna_string.as_bytes());
        let index_r = minimal_rotation(dna_string_r.as_bytes());
        // create the substrings
        // starts
        let start_f = &dna_string[index_f..];
        let start_r = &dna_string[index_r..];
        // ends
        let end_f = &dna_string[0..index_f];
        let end_r = &dna_string[0..index_r];
        // string
        let lms_f = format!("{}{}", start_f, end_f);
        let lms_r = format!("{}{}", start_r, end_r);
        let mut strings = vec![&lms_f, &lms_r];
        strings.string_sort_unstable(natural_lexical_cmp);
        strings[0].to_string()
    }
}
