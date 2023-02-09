use bioshell_sim::{ResizableSystem, System};
use bioshell_core::sequence::ResidueTypeOrder;
use crate::coupling_energy::CouplingEnergy;

#[derive(Clone)]
pub struct EvolvingSequence {
    /// Most recent total energy of the sequence system
    pub total_energy: f64,
    /// Protein or nucleic acid  sequence - encoded as `u8` integers
    /// Length of this vector should match the length of the evolving sequence
    pub sequence: Vec<u8>,
    pub res_mapping: ResidueTypeOrder,
    current_size: usize
}

impl EvolvingSequence {
    pub fn new(starting_sequence: &String, res_mapping: ResidueTypeOrder) -> EvolvingSequence {
        let n = starting_sequence.len();
        let sequence: Vec<u8> = Vec::new();
        let mut out = EvolvingSequence { total_energy: 0.0, sequence, res_mapping, current_size: n};
        out.sequence = out.encode_sequence(starting_sequence);

        return out;
    }

    pub fn seq_len(&self) -> usize { self.sequence.len() }

    pub fn decode_other(&self, system: &Vec<u8>) -> String {
        let mut buffer: Vec<char> = Vec::new();
        buffer.reserve(system.len());
        for aai in system {
            buffer.push(self.res_mapping.decode_letter(*aai));
        }
        buffer.iter().collect()
    }

    pub fn aa_cnt(&self) -> usize { self.res_mapping.size() }

    pub fn decode_sequence(&self) -> String {
        self.decode_other(&self.sequence)
    }

    pub fn encode_sequence(&self, seq: &String) -> Vec<u8> {
        let mut buffer: Vec<u8> = Vec::new();
        buffer.reserve(seq.len());
        for aai in seq.chars() {
            buffer.push(self.res_mapping.encode_letter(&aai));
        }
        return buffer;
    }

    pub fn delta_energy(&self, pos: usize, new_aa: u8, energy: &CouplingEnergy) -> f64 {

        let mut en: f32 = 0.0;
        let pos_i: usize = pos * self.aa_cnt();
        let pos_i_old = pos_i + self.sequence[pos] as usize;
        let pos_i_new = pos_i + new_aa as usize;
        let mut pos_j: usize = 0;
        for aa_j in self.sequence.iter() {
            en += energy.cplngs.data[pos_j + *aa_j as usize][pos_i_new];
            en -= energy.cplngs.data[pos_j + *aa_j as usize][pos_i_old];
            pos_j += self.aa_cnt();
        }

        return en as f64;
    }
}


impl System for EvolvingSequence {
    fn size(&self) -> usize { self.current_size }

    fn copy_from(&mut self, i: usize, rhs: &Self) { self.sequence[i] = rhs.sequence[i]; }
}

impl ResizableSystem for EvolvingSequence {
    /// TODO: throw excption if new_size is longer than the maximum size
    fn set_size(&mut self, new_size: usize) { self.current_size = new_size; }

    fn capacity(&self) -> usize { self.seq_len() }
}