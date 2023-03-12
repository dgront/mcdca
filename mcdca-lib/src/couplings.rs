use std::fmt;
use std::slice::Iter;

use bioshell_core::sequence::Sequence;
use crate::evolving_sequence::SequenceEntry;
use crate::EvolvingSequence;
use crate::pseudocounts::Pseudocounts;

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
                if j % self.k_aa_types == (self.k_aa_types - 1) { write!(f, "  ").ok(); }
            }
            writeln!(f, "");
            if i % self.k_aa_types == (self.k_aa_types - 1) { writeln!(f, "#").ok(); }
        }
        writeln!(f, "# scaled by: {factor}").ok();
        writeln!(f, "# total: {total}").ok();
        Ok(())
    }
}

/// Initializes coupling from a given set of sequences.
///
pub fn counts_from_weighted_sequences(system: &EvolvingSequence, seq_pool: Iter<'_, SequenceEntry>) -> Couplings {

    let mut weights: Vec<f64> = Vec::new();
    let mut msa: Vec<Sequence> = Vec::new();

    let seq_id: String = "".parse().unwrap();
    for sequence in seq_pool {
        weights.push(sequence.energy);
        msa.push(Sequence::new(&seq_id, &sequence.sequence));
    }

    return counts_from_weighted_msa(system, &msa, &weights);
}

pub fn counts_from_weighted_msa(system: &EvolvingSequence, msa: &Vec<Sequence>, weights: &Vec<f64>) -> Couplings {
    let n = system.seq_len();
    let k = system.res_mapping.size();
    let mut cplngs: Couplings = Couplings::new(n, k);

    let mut tmp_i: Vec<usize> = vec![0; n];       // Temporary indexes aa->cplngs matrix
    let mut i_seq: usize = 0;
    let total_weight: f64 = weights.iter().sum();
    for sequence in msa {
        if sequence.len() != n {
            eprintln!("\nSequence of incorrect length! Is: {}, should be: {}. The sequence skipped::\n {}\n",
                      sequence.len(), n, sequence);
            continue;
        }
        // --- for every position in a sequence, find the location of that AA in the big cplngs matrix
        // --- This is done in linear time here and called below in a square loop
        for idx in 0..n {
            tmp_i[idx] = idx * k + system.res_mapping.type_to_index(&sequence.seq()[idx]) as usize;
        }
        // --- Count each pair coincidence
        let w = 0.5 * weights[i_seq] / total_weight;
        for i_pos in 1..n {
            for j_pos in 0..i_pos {
                cplngs.data[tmp_i[i_pos]][tmp_i[j_pos]] += w as f32;
                cplngs.data[tmp_i[j_pos]][tmp_i[i_pos]] += w as f32;
            }
        }
        i_seq += 1;
    }
    cplngs.normalize(i_seq as f32);   // --- times 2.0 because sequences are loaded twice: to upper and lower triangle of the matrix

    return cplngs;
}

/// Initializes coupling from a given set of aligned sequences.
///
pub fn counts_from_msa(system: &EvolvingSequence, msa: &Vec<Sequence>) -> Couplings {

    let weights = vec![1.0; msa.len()];
    return counts_from_weighted_msa(system, &msa, &weights);
}


/// Update the `couplings` matrix to lower the distance between expected and target observations
///
/// ```math
/// j^{k,l}_{A,B} = j^{k,l}_{A,B} - \text{c} * \frac{p^{k,l}_{A,B} - \hat{p}^{k,l}_{A,B}}{p^{k,l}_{A,B}}
/// ```
pub fn update_couplings(target_counts: &Couplings, observed_counts: &Couplings, pseudocounts: &Pseudocounts,
                        newton_step: f64, simulated_couplings: &mut Couplings) -> f64 {

    let mut error: f64 = 0.0;
    let k_aa = target_counts.k_aa_types;
    let n_res = target_counts.seq_length;
    for i_pos in 0..n_res {
        for i_aa in 0..k_aa {
            let i_ind = i_pos * k_aa + i_aa;
            for j_pos in 0..n_res {
                if i_pos == j_pos { continue }
                for j_aa in 0..k_aa {
                    let mut pseudo: f64 = pseudocounts.pseudo_fraction(i_pos, i_aa, j_pos, j_aa);
                    let j_ind = j_pos * k_aa + j_aa;
                    let target_val = target_counts.data[i_ind][j_ind] as f64;
                    let observed_val = observed_counts.data[i_ind][j_ind] as f64;
                    let mut delta = target_val - observed_val;
                    if observed_val.abs() < 1e-10 && target_val.abs() < 1e-10 {
                        continue;
                    }
                    error += delta*delta;
                    let step = newton_step * delta / (2.0 * observed_val + pseudo);
                    if delta.abs() > 10000.0 {
                        println!("{:3} {:2} {:3} {:2}  should be: {:5.3}  observed: {:5.3}, J: {:5.3} -> {:5.3} by {}, delta: {}, err: {}",
                                 i_pos, i_aa, j_pos, j_aa,
                                 target_val, observed_val,
                                 simulated_couplings.data[i_ind][j_ind], simulated_couplings.data[i_ind][j_ind] as f64 - step, step, delta, delta * delta);
                    }
                    simulated_couplings.data[i_ind][j_ind] -= step as f32;
                }
            }
        }
    }
    return error;
}