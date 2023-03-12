use std::ops::Range;
use bioshell_sim::{Energy};

use crate::Couplings;
use crate::evolving_sequence::EvolvingSequence;

pub struct CouplingEnergy {
    pub cplngs: Couplings,
}

impl CouplingEnergy {

    /// Provides mutable access to the coupling matrix used by this energy function
    pub fn get_couplings_mut(&mut self) -> &mut Couplings { &mut self.cplngs }

    /// Prints each contribution to the total energy of this system
    pub fn explain(&self, system: &EvolvingSequence) {
        let mut en: f32 = 0.0;
        let mut e: f32 = 0.0;
        let mut pos_i: usize = 0;
        for aa_i in system.sequence.iter() {
            let mut pos_j: usize = 0;
            for aa_j in system.sequence.iter() {
                e = self.cplngs.data[pos_i + *aa_i as usize][pos_j + *aa_j as usize];
                en += e;
                println!("{} {} {} {}", pos_i, pos_j, e, en);
                pos_j += system.aa_cnt();
            }
            pos_i += system.aa_cnt();
        }
    }
}

impl Energy<EvolvingSequence> for CouplingEnergy {

    fn energy(&self, system: &EvolvingSequence) -> f64 {
        let mut en: f32 = 0.0;
        let mut pos_i: usize = 0;
        for aa_i in system.sequence.iter() {
            let mut pos_j: usize = 0;
            for aa_j in system.sequence.iter() {
                en += self.cplngs.data[pos_i + *aa_i as usize][pos_j + *aa_j as usize];
                pos_j += system.aa_cnt();
            }
            pos_i += system.aa_cnt();
        }
        return (en / 2.0) as f64;
    }

    fn energy_by_pos(&self, system: &EvolvingSequence, pos: usize) -> f64 {
        let mut en: f32 = 0.0;
        let pos_i: usize = pos * system.aa_cnt();
        let aa_i = system.sequence[pos];
        let mut pos_j: usize = 0;
        for aa_j in system.sequence.iter() {
            en += self.cplngs.data[pos_j + *aa_j as usize][pos_i + aa_i as usize];
            pos_j += system.aa_cnt();
        }

        return en as f64;
    }

    fn energy_by_range(&self, system: &EvolvingSequence, range: &Range<usize>) -> f64 {
        todo!()
    }

    fn name(&self) -> String { String::from("CouplingEnergy") }
}