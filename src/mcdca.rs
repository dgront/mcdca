use rand::Rng;
use std::collections::HashMap;

use bioshell_core::Sequence;


pub struct Couplings {
    pub n: usize,
    pub k: usize,
    cplngs: Vec<Vec<f32>>,
    index_to_aa: Vec<u8>,
    aa_to_index: HashMap<u8, usize>,
}

impl Couplings {
    /// Creates an empty Coupling instance i.e. none of amino acids are coupled
    pub fn new(seq_len: usize, aa_order: &str) -> Couplings {
        let m = vec![vec![0.0; seq_len * aa_order.len()]; seq_len * aa_order.len()];
        let index_to_aa = aa_order.as_bytes().to_vec();
        println!("amino acid order: {:?}",index_to_aa);
        let mut aa_to_index = HashMap::new();
        for i in 0..aa_order.len() {
            aa_to_index.insert(index_to_aa[i],i);
        }
        Couplings { n: seq_len, k: aa_order.len(), cplngs: m, index_to_aa, aa_to_index }
    }

    /// Initializes coupling diagonally.
    ///
    pub fn init_couplings_diagonaly(&mut self) {
        for i in 1..self.n {
            let ii = i * self.k;
            for j in 0..self.k {
                self.cplngs[ii + j][ii - self.k + j] = -1.0;
                self.cplngs[ii - self.k + j][ii + j] = -1.0;
            }
        }
        for j in 0..self.k {
            self.cplngs[j][self.k + j] = -1.0;
            self.cplngs[self.k + j][j] = -1.0;
        }
    }

    /// Initializes coupling from a given set of sequences.
    ///
    pub fn init_couplings_by_msa(&mut self, msa: &Vec<Sequence>) {
        for sequence in msa {
            for i_pos in 1..self.n {
                let i_aa = sequence.aa(i_pos);
                let ii = i_pos * self.k + self.aa_to_index[&i_aa];
                for j_pos in 0..i_pos {
                    let j_aa = self.aa_to_index[&sequence.aa(j_pos)];
                    self.cplngs[ii][j_pos * self.k + j_aa] += -1.0;
                }
            }
        }
        self.cplngs.iter_mut().for_each(|el| el.iter_mut().for_each(|iel| *iel /= msa.len() as f32));
    }

    /// Prints the large matrix of couplings on the screen
    pub fn show_matrix(n:usize, k:usize, m: &Vec<Vec<f32>>) {
        for (i, row) in m.iter().enumerate() {
            for (j, val) in row.iter().enumerate() {
                print!(" {:.3}",val);
                if j % k == (k - 1) { print!(" ") }
            }
            println!("");
            if i % k == (k - 1) { println!("#") }
        }
    }

    /// Prints the large matrix of couplings on the screen
    pub fn show(&self) { Couplings::show_matrix(self.n, self.k, &self.cplngs); }

    pub fn decode_sequence(&self, system: &Vec<u8>) -> String {
        let mut buffer: Vec<u8> = Vec::new();
        buffer.reserve(system.len());
        for i in 0..system.len() {
            buffer.push(self.index_to_aa[system[i] as usize]);
        }
        String::from_utf8_lossy(&buffer).to_string()
    }

    pub fn delta_energy(&self, system: &Vec<u8>, pos: usize, old: usize, new: usize) -> f32 {
        let mut en: f32 = 0.0;
        let mut pos_j: usize = 0;
        let pos_i: usize = pos * self.k;
        for aa_j in system.iter() {
            en += self.cplngs[pos_j + *aa_j as usize][pos_i + new] - self.cplngs[pos_j + *aa_j as usize][pos_i + old];
            // println!("{} {} {} {} {}", pos_j, aa_j, cplngs[pos_j + *aa_j as usize][pos_i + new], cplngs[pos_j + *aa_j as usize][pos_i + old], en);
            pos_j+=self.k;
        }

        return en;
    }

    pub fn total_energy(&self, system: &Vec<u8>) -> f32 {
        let mut en:f32 = 0.0;
        let mut pos_i :usize = 0;
        for aa_i in system.iter() {
            let mut pos_j :usize = 0;
            for aa_j in system.iter() {
                en += self.cplngs[pos_i + *aa_i as usize][pos_j + *aa_j as usize];
                pos_j += self.k;
            }
            pos_i += self.k;
        }
        return en/2.0;
    }
}


pub fn accumulate_counts(system: &Vec<u8>, n_aa:usize, counts: &mut Vec<Vec<f32>>) {
    let mut pos_i: usize = 0;
    for aa_i in system.iter() {
        let mut pos_j: usize = 0;
        for aa_j in system.iter() {
            counts[pos_i + *aa_i as usize][pos_j + *aa_j as usize] += 1.0;
            pos_j += n_aa;
        }
        pos_i += n_aa;
    }
}

pub fn isothermal_mc(system: &mut Vec<u8>, energy: &Couplings, inner_cycles: i32, outer_cycles: i32) -> Vec<Vec<f32>>{

    let mut counts: Vec<Vec<f32>> = vec![vec![0.0; energy.n * energy.k]; energy.n * energy.k];
    let mut rng = rand::thread_rng();
    let mut total_en:f32 = energy.total_energy(&system);
    let mut n_obs: f32 = 0.0;
    let mut n_succ: f32 = 0.0;
    for io in 0..outer_cycles {
        for ii in 0..inner_cycles {
            for is in 0..system.len() {
                let pos: usize = rng.gen_range(0..energy.n);
                let new_aa: usize = rng.gen_range(0..energy.k);
                let delta_en = energy.delta_energy(system, pos, system[pos] as usize, new_aa);
                if delta_en > 0.0 && (-delta_en).exp() < rng.gen_range(0.0..1.0) { continue; }
                system[pos] = new_aa as u8;
                total_en += delta_en;
                n_succ += 1.0;
            }
        }
        n_obs += 1.0;
        accumulate_counts(&system, energy.k, &mut counts);
        println!("{:4} {}",total_en, energy.decode_sequence(&system));
    }
    println!("#{}",n_succ / (outer_cycles*inner_cycles*system.len() as i32) as f32);
    counts.iter_mut().for_each(|el| el.iter_mut().for_each(|iel| *iel /= n_obs));
    return counts;
}

pub fn main() {
    let seq_len: usize = 100;
    let mut system: Vec<u8> = vec![0; seq_len];

    let mut en: Couplings = Couplings::new(seq_len, "ACDEFGHIKLMNPQRSTVWY-");
    // let mut en: Couplings = Couplings::new(seq_len, "ACD");

    en.init_couplings_diagonaly();

    en.total_energy(&system);
    let counts = isothermal_mc(&mut system, &en,10,100000);
    // en.show();
    Couplings::show_matrix(en.n, en.k, &counts);
}