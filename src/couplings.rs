use std::fmt;
use std::collections::HashMap;

use bioshell_core::sequence::Sequence;

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
        self.data.iter_mut().for_each(|el| el.iter_mut().for_each(|iel| *iel /= cnt))
    }

    /// Sets all coupling values to 0.0
    pub fn clear(&mut self) {
        let nk = self.n*self.k;
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
    pub total_energy: f32,
    pub sequence: Vec<u8>,
    index_to_aa: Vec<char>,
    aa_to_index_map: HashMap<char, u8>,
    pub(crate) cplngs: Couplings,
}

impl EvolvingSequence {
    pub fn new(starting_sequence: &String, aa_order: &str) -> EvolvingSequence {

        let n = starting_sequence.len();
        let k = aa_order.len();
        let cplngs = Couplings::new(n, k);
        let sequence: Vec<u8> = Vec::new();
        let index_to_aa = aa_order.chars().collect();
        let mut aa_to_index_map: HashMap<char, u8> = HashMap::new();
        for (i, aai) in aa_order.chars().enumerate() {
            aa_to_index_map.insert(aai, i as u8);
        }
        let mut cplngs = Couplings::new(n, k);
        init_couplings_diagonally(&mut cplngs);
        let mut out = EvolvingSequence { total_energy: 0.0, sequence, index_to_aa, aa_to_index_map, cplngs };
        out.sequence = out.encode_sequence(starting_sequence);
        out.total_energy = out.energy();

        return out;
    }

    pub fn seq_len(&self) -> usize { self.cplngs.n }

    pub fn aa_cnt(&self) -> usize { self.cplngs.k }

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
}

/// Initializes coupling from a given set of sequences.
///
pub fn counts_from_msa(system: &EvolvingSequence, msa: &Vec<Sequence>) -> Couplings {

    if !check_msa(msa) { std::process::exit(1);}

    let n = system.seq_len();
    let k = system.aa_cnt();
    let mut cplngs: Couplings = Couplings::new(n, k);

    let mut tmp_i: Vec<usize> = vec![0; n];       // Temporary indexes aa->cplngs matrix
    for sequence in msa {
        // --- for every position in a sequence, find the location of that AA in the big cplngs matrix
        // --- This is done in linear time here and called below in a square loop
        for idx in 0..n {
            tmp_i[idx] = idx * k + system.aa_to_index(&sequence.char(idx)) as usize;
        }
        // --- Count each pair coincidence
        for i_pos in 1..n {
            for j_pos in 0..i_pos {
                cplngs.data[tmp_i[i_pos]][tmp_i[j_pos]] += 1.0;
                cplngs.data[tmp_i[j_pos]][tmp_i[i_pos]] += 1.0;
            }
        }
    }
    cplngs.normalize(msa.len() as f32);

    return cplngs;
}

fn check_msa(msa: &Vec<Sequence>) -> bool {
    let first: &Sequence = &msa[0];
    for i in 1..msa.len() {
        if msa[i].len() != first.len() {
            eprintln!("\nSequence of incorrect length:\n {}\n",msa[i])
        }
    }
    true
}