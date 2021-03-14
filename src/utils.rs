pub mod utils {
    use bio::pattern_matching::bom::BOM;
    // I'd like to use shift_and, but may have to wait until next public release
    // of rust-bio (or hard code it here...)
    // see https://github.com/rust-bio/rust-bio/blob/master/src/pattern_matching/shift_and.rs
    use bio::pattern_matching::kmp::KMP;

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

    // string rotation
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

    //https://stackoverflow.com/questions/50380352/how-can-i-group-consecutive-integers-in-a-vector-in-rust
    // group consecutive numbers in a slice into their own slice in a vec.
    fn consecutive_slices(data: &[usize]) -> Vec<&[usize]> {
        let mut slice_start = 0;
        let mut result = Vec::new();
        for i in 1..data.len() {
            if data[i - 1] + 1 != data[i] {
                result.push(&data[slice_start..i]);
                slice_start = i;
            }
        }
        if data.len() > 0 {
            result.push(&data[slice_start..]);
        }
        result
    }

    // format a putative telomeric repeat, which could be any one of
    // several string rotations. Here I try to standardise the notation
    // by reporting the longest consecutive string of G's as the last part
    // of the telomeric repeat. This may be totally the wrong thing to do...
    pub fn format_telomeric_repeat(telomeric_repeat: String) -> String {
        // borrow to get &str
        let tr = &telomeric_repeat;

        // match the indices of C's and G's
        // result is a tuple, so collect the first element (indexes) only
        let c_m: Vec<usize> = tr.match_indices("C").map(|x| x.0).collect();
        let g_m: Vec<usize> = tr.match_indices("G").map(|x| x.0).collect();

        // get a vec of slices of consecutive matches
        let c_slices = consecutive_slices(&c_m);
        let g_slices = consecutive_slices(&g_m);

        // get the lengths of each of the slices in the vec.
        let c_slices_lens: Vec<_> = c_slices.iter().map(|e| e.len()).collect();
        let g_slices_lens: Vec<_> = g_slices.iter().map(|e| e.len()).collect();

        // now get the maximum length from each vec, as this will
        // determine what to do next.
        let max_cs_lens = match c_slices_lens.iter().max() {
            Some(x) => x,
            None => &0usize,
        };
        let max_gs_lens = match g_slices_lens.iter().max() {
            Some(x) => x,
            None => &0usize,
        };

        // if there are no G's or the C's outnumber the G's,
        // put the C's to the start and take the reverse complement.
        // and of course, the alternate arm requires the opposite conditions.
        if max_gs_lens == &0usize || max_cs_lens > max_gs_lens {
            // means we want the revcomp of this string, so put C's at the start of the string.
            let mut indices_lens = Vec::new();
            for i in c_slices {
                indices_lens.push((i, i.len()));
            }
            // because I know there should be at least one element, unwrap should be fine..?
            let max = indices_lens.iter().max_by_key(|i| i.1).unwrap();
            // need everything up until the first C, then shove it to the end.
            let index = max.0.first().unwrap();
            // if the index is at the start, we don't need to shuffle the string
            // is zero right or &0usize?
            if index == &0usize {
                return reverse_complement(tr);
            } else {
                let start = &tr[*index..];
                let end = &tr[..*index];
                let res = format!("{}{}", start, end);
                return reverse_complement(&res);
            }
        } else {
            // means G's go at the end of the string
            // means we want the revcomp of this string, so put C's at the start of the string.
            let mut indices_lens = Vec::new();
            for i in g_slices {
                indices_lens.push((i, i.len()));
            }
            // because I know there should be at least one element, unwrap should be fine..?
            let max = indices_lens.iter().max_by_key(|i| i.1).unwrap();
            // need everything up until the first C, then shove it to the end.
            let index = max.0.last().unwrap();
            // if the index is at the start, we don't need to shuffle the string
            if index == &tr.len() {
                return tr.to_owned();
            } else {
                let start = &tr[*index + 1..];
                let end = &tr[..*index];
                let res = format!("{}{}", start, end);
                return res;
            }
        }
    }
}
