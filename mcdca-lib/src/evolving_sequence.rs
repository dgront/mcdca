use bioshell_sim::{ResizableSystem, System};
use bioshell_core::sequence::ResidueTypeOrder;
use crate::coupling_energy::CouplingEnergy;


/// Stores a sequence and its energy
///
/// In the case of PERM sampling, the ``energy`` field is used to store Boltzmann weight rather than the energy itself
pub struct SequenceEntry {
    pub(crate) energy: f64,
    pub(crate) weight: f64,
    pub(crate) sequence: String
}

#[derive(Clone)]
pub struct EvolvingSequence {
    /// Most recent total energy of the sequence system
    pub total_energy: f64,
    /// Statistical weight of the current state of this sequence
    pub weight: f64,
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
        let mut out = EvolvingSequence { total_energy: 0.0, weight: 1.0, sequence, res_mapping, current_size: n};
        out.sequence = out.encode_sequence(starting_sequence);

        return out;
    }

    pub fn seq_len(&self) -> usize { self.sequence.len() }

    pub fn decode_other(&self, system: &Vec<u8>) -> String {
        let mut buffer: Vec<char> = Vec::new();
        buffer.reserve(system.len());
        for aai in system {
            buffer.push(self.res_mapping.index_to_letter(*aai));
        }
        buffer.iter().collect()
    }

    /// Returns the alphabet size for this sequence.
    ///
    /// The alphabet size is typically 20 for proteins and 4 for nucleic acids (21 and 5 respectively,
    /// when a gap symbol is also allowed in a sequence)
    pub fn aa_cnt(&self) -> usize { self.res_mapping.size() }

    pub fn decode_sequence(&self) -> String {
        self.decode_other(&self.sequence)
    }

    pub fn encode_sequence(&self, seq: &String) -> Vec<u8> {

        let types: Vec<u8> = seq.to_owned().into_bytes();
        let mut buffer: Vec<u8> = vec![];
        for aai in types {
            buffer.push(self.res_mapping.type_to_index(&aai) as u8);
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

    pub fn energy_by_letter(&self, pos: usize, energy: &CouplingEnergy, results: &mut Vec<f32>) {

        for i in 0..self.aa_cnt() { results[i] = 0.0; }

        let mut pos_j: usize = 0;
        let pos_i: usize = pos * self.aa_cnt();

        for aa_j in self.sequence.iter() {
            for i in 0..self.aa_cnt() {
                results[i] += energy.cplngs.data[pos_j + *aa_j as usize][pos_i+i];
            }
            pos_j += self.aa_cnt();
        }
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