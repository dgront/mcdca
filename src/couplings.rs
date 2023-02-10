use std::fmt;

use bioshell_core::sequence::Sequence;
use crate::EvolvingSequence;
// use crate::evolving_sequence::EvolvingSequence;

#[derive(Clone, Debug)]
pub struct Couplings {
    pub seq_length: usize,
    pub k_aa_types: usize,
    pub data: Vec<Vec<f32>>,
}

impl Couplings {
    /// Creates an empty Coupling instance i.e. none of the positions are coupled
    pub fn new(seq_len: usize, n_aa_types: usize) -> Couplings {
        let size: usize = seq_len * n_aa_types;
        let m = vec![vec![0.0; size]; size];
        Couplings { seq_length: seq_len, k_aa_types: n_aa_types, data: m }
    }

    pub fn normalize(&mut self, cnt: f32) {
        self.data.iter_mut().for_each(|el| el.iter_mut().for_each(|iel| *iel /= cnt))
    }

    /// Find the minimum coupling value
    pub fn min(&self) -> f32 {
        let mut m = self.data[0][0];
        let nk = self.seq_length * self.k_aa_types;
        for i in 0..nk {
            // use this: .iter().min().unwrap();
            for j in 0..nk {
                m = m.min(self.data[i][j])
            }
        }
        return m;
    }

    /// Find the maximum coupling value
    pub fn max(&self) -> f32 {
        let mut m = self.data[0][0];
        let nk = self.seq_length * self.k_aa_types;
        for i in 0..nk {
            for j in 0..nk {
                m = m.max(self.data[i][j])
            }
        }
        return m;
    }

    /// Sets all coupling values to 0.0
    pub fn clear(&mut self) {
        let nk = self.seq_length *self.k_aa_types;
        for i in 0..nk {
            for j in 0..nk {
                self.data[i][j] = 0.0;
            }
        }
    }
}

impl fmt::Display for Couplings {
    /// Creates a `String` representation of a given `Couplings`
    /// # Examples
    ///
    /// Create diagonal `Couplings` and turn it into a string
    ///
    /// ```rust
    /// use std::fmt::Write;
    /// // create an empty Couplings object: three-letter sequence with 4-letter alphabet (e.g. RNA)
    /// let mut cplngs = Couplings::new(3,4);
    /// init_couplings_diagonally(&cplngs);
    /// println!("{}", cplngs);
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let mut total: f32 = 0.0;
        let power = self.max().log10().floor() as f32;
        let factor:f32 = (10.0_f32).powf(power) as f32;
        for (i, row) in self.data.iter().enumerate() {
            for (j, val) in row.iter().enumerate() {
                write!(f, "{:5.2} ", val/factor);
                total += val;
                if j % self.k_aa_types == (self.k_aa_types - 1) { write!(f, "  "); }
            }
            writeln!(f, "");
            if i % self.k_aa_types == (self.k_aa_types - 1) { writeln!(f, "#"); }
        }
        writeln!(f, "# scaled by: {factor}");
        writeln!(f, "# total: {total}");
        Ok(())
    }
}




/// Initializes coupling from a given set of sequences.
///
pub fn counts_from_msa(system: &EvolvingSequence, msa: &Vec<Sequence>) -> Couplings {

    let n = system.seq_len();
    let k = system.res_mapping.size();
    let mut cplngs: Couplings = Couplings::new(n, k);

    let mut tmp_i: Vec<usize> = vec![0; n];       // Temporary indexes aa->cplngs matrix
    let mut n_seq: f32 = 0.0;
    for sequence in msa {
        if sequence.len() != n {
            eprintln!("\nSequence of incorrect length! Is: {}, should be: {}. The sequence skipped::\n {}\n",
                      sequence.len(), n, sequence);
            continue;
        }
        // --- for every position in a sequence, find the location of that AA in the big cplngs matrix
        // --- This is done in linear time here and called below in a square loop
        for idx in 0..n {
            tmp_i[idx] = idx * k + system.res_mapping.encode_letter(&sequence.char(idx)) as usize;
        }
        // --- Count each pair coincidence
        for i_pos in 1..n {
            for j_pos in 0..i_pos {
                cplngs.data[tmp_i[i_pos]][tmp_i[j_pos]] += 0.5;
                cplngs.data[tmp_i[j_pos]][tmp_i[i_pos]] += 0.5;
            }
        }
        n_seq += 1.0;
    }
    cplngs.normalize(n_seq as f32);   // --- times 2.0 because sequences are loaded twice: to upper and lower triangle of the matrix

    return cplngs;
}
