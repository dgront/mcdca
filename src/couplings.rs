use std::fmt;
use rand::Rng;
use std::collections::HashMap;

use bioshell_core::Sequence;

pub struct Couplings {
    pub n: usize,
    pub k: usize,
    pub data: Vec<Vec<f32>>,
}

impl Couplings {
    /// Creates an empty Coupling instance i.e. none of the positions are coupled
    pub fn new(seq_len: usize, n_aa_types: usize) -> Couplings {
        let size: usize = seq_len * n_aa_types;
        let m = vec![vec![0.0; size]; size];
        Couplings { n: seq_len, k: n_aa_types, data: m }
    }

    pub fn normalize(&mut self, cnt: f32) {
        self.data.iter_mut().for_each(|el| el.iter_mut().for_each(|iel| *iel /= cnt as f32))
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
        for (i, row) in self.data.iter().enumerate() {
            for (j, val) in row.iter().enumerate() {
                write!(f, "{:.4} ", val);
                if j % self.k == (self.k - 1) { write!(f, "  "); }
            }
            writeln!(f, "");
            if i % self.k == (self.k - 1) { writeln!(f, "#"); }
        }
        Ok(())
    }
}

/// Initializes coupling diagonally.
///
pub fn init_couplings_diagonally(cplngs: &mut Couplings) {

    let n = cplngs.n;
    let k = cplngs.k;
    for i in 1..n {
        let ii = i * k;
        for j in 0..k {
            cplngs.data[ii + j][ii - k + j] = -1.0;
            cplngs.data[ii - k + j][ii + j] = -1.0;
        }
    }
    for j in 0..k {
        cplngs.data[j][k + j] = -1.0;
        cplngs.data[k + j][j] = -1.0;
    }
}


// -----------------------------------------------------
pub struct EvolvingSequence {
    pub(crate) total_energy: f32,
    pub(crate) sequence: Vec<u8>,
    index_to_aa: Vec<char>,
    aa_to_index: HashMap<char, u8>,
    cplngs: Couplings,
}

impl EvolvingSequence {
    pub fn new(starting_sequence: &String, aa_order: &str) -> EvolvingSequence {

        let n = starting_sequence.len();
        let k = aa_order.len();
        let cplngs = Couplings::new(n, k);
        let sequence: Vec<u8> = Vec::new();
        let index_to_aa = aa_order.chars().collect();
        let mut aa_to_index: HashMap<char, u8> = HashMap::new();
        for (i, aai) in aa_order.chars().enumerate() {
            aa_to_index.insert(aai, i as u8);
        }
        let mut cplngs = Couplings::new(n, k);
        init_couplings_diagonally(&mut cplngs);
        let mut out = EvolvingSequence { total_energy: 0.0, sequence, index_to_aa, aa_to_index, cplngs };
        out.sequence = out.encode_sequence(starting_sequence);
        out.total_energy = out.energy();

        return out;
    }

    pub fn seq_len(&self) -> usize { self.cplngs.n }

    pub fn aa_cnt(&self) -> usize { self.cplngs.k }

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
            buffer.push(self.aa_to_index[&aai]);
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
}

/// Initializes coupling from a given set of sequences.
///
pub fn init_couplings_by_msa(system: &mut EvolvingSequence, msa: &Vec<Sequence>) {

    for sequence in msa {
        let n = system.seq_len();
        let k = system.aa_cnt();

        for i_pos in 1..n {
            let i_aa = sequence.char(i_pos);
            let ii = i_pos * k + system.aa_to_index[&i_aa] as usize;
            for j_pos in 0..i_pos {
                let j_aa: u8 = system.aa_to_index[&sequence.char(j_pos)];
                system.cplngs.data[ii][j_pos * k + j_aa as usize] += -1.0;
            }
        }
    }
    system.cplngs.normalize(msa.len() as f32);
}
