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
    /// // create an empty Couplings object
    /// let mut cplngs = Couplings::new(10,2);
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
                cplngs.data[tmp_i[i_pos]][tmp_i[j_pos]] += 1.0;
                cplngs.data[tmp_i[j_pos]][tmp_i[i_pos]] += 1.0;
            }
        }
        n_seq += 1.0;
    }
    cplngs.normalize(n_seq * n as f32);

    return cplngs;
}

/*
pub struct EvolvingSequence2 {
    pub total_energy: f32,
    pub sequence: Vec<u8>,
    index_to_aa: Vec<char>,
    aa_to_index_map: HashMap<char, u8>,
    pub cplngs: Couplings,
    pub observed: Couplings,
}

impl EvolvingSequence2 {
    pub fn new(starting_sequence: &String, aa_order: &str) -> EvolvingSequence2 {

        let n = starting_sequence.len();
        let k = aa_order.len();
        let sequence: Vec<u8> = Vec::new();
        let index_to_aa = aa_order.chars().collect();
        let mut aa_to_index_map: HashMap<char, u8> = HashMap::new();
        for (i, aai) in aa_order.chars().enumerate() {
            aa_to_index_map.insert(aai, i as u8);
        }
        let mut cplngs = Couplings::new(n, k);
        let observed = Couplings::new(n, k);
        let mut out = EvolvingSequence2 { total_energy: 0.0, sequence, index_to_aa, aa_to_index_map, cplngs, observed };
        out.sequence = out.encode_sequence(starting_sequence);
        out.total_energy = out.energy();

        return out;
    }

    pub fn seq_len(&self) -> usize { self.cplngs.seq_length }

    pub fn aa_cnt(&self) -> usize { self.cplngs.k_aa_types }

    /// Convert a given amino acid character to its internal index
    pub fn aa_to_index(&self, aa:&char) -> u8 {
        let aa_id: &u8 =  match self.aa_to_index_map.get(&aa) {
            Some(i) => { i },
            None => {
                eprintln!("Unknown amino acid symbol {}, converted to gap", &aa);
                &self.aa_to_index_map[&'-']
            }
        };

        *aa_id
    }

    pub fn decode_other(&self, system: &Vec<u8>) -> String {
        let mut buffer: Vec<char> = Vec::new();
        buffer.reserve(system.len());
        for aai in system {
            buffer.push(self.index_to_aa[*aai as usize]);
        }
        buffer.iter().collect()
    }


    pub fn decode_sequence(&self) -> String {
        self.decode_other(&self.sequence)
    }

    pub fn encode_sequence(&self, seq: &String) -> Vec<u8> {
        let mut buffer: Vec<u8> = Vec::new();
        buffer.reserve(seq.len());
        for aai in seq.chars() {
            buffer.push(self.aa_to_index(&aai));
        }
        return buffer;
    }

    pub fn energy(&self) -> f32 {

        let mut en:f32 = 0.0;
        let mut pos_i :usize = 0;
        for aa_i in self.sequence.iter() {
            let mut pos_j :usize = 0;
            for aa_j in self.sequence.iter() {
                en += self.cplngs.data[pos_i + *aa_i as usize][pos_j + *aa_j as usize];
                pos_j += self.aa_cnt();
            }
            pos_i += self.aa_cnt();
        }
        return en/2.0;
    }

    pub fn energy_by_pos(&self, pos: usize) -> f32 {
        let mut en: f32 = 0.0;
        let pos_i: usize = pos * self.aa_cnt();
        let aa_i = self.sequence[pos];
        let mut pos_j: usize = 0;
        for aa_j in self.sequence.iter() {
            en += self.cplngs.data[pos_j + *aa_j as usize][pos_i + aa_i as usize];
            pos_j += self.aa_cnt();
        }

        return en;
    }

    pub fn delta_energy(&self, pos: usize, new_aa: u8) -> f32 {

        let mut en: f32 = 0.0;
        let pos_i: usize = pos * self.aa_cnt();
        let pos_i_old = pos_i + self.sequence[pos] as usize;
        let pos_i_new = pos_i + new_aa as usize;
        let mut pos_j: usize = 0;
        for aa_j in self.sequence.iter() {
            en += self.cplngs.data[pos_j + *aa_j as usize][pos_i_new];
            en -= self.cplngs.data[pos_j + *aa_j as usize][pos_i_old];
            pos_j += self.aa_cnt();
        }

        return en;
    }

    pub fn energy_row(&self, pos: usize, energy: &mut Vec<f32>)  {

        for aa_i in 0..self.cplngs.k_aa_types { energy[aa_i] = 0.0; }
        let pos_i: usize = pos * self.aa_cnt();
        let mut pos_j: usize = 0;
        for aa_j in self.sequence.iter() {
            let row: &Vec<f32> = &self.cplngs.data[pos_j + *aa_j as usize];
            for aa_i in 0..self.aa_cnt() {
                energy[aa_i] += row[pos_i + aa_i];
            }
            pos_j += self.aa_cnt();
        }
    }
}
*/
