use std::ops::Range;
use rand::{Rng, SeedableRng};

use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover};
use bioshell_sim::{Observer, ObserversSet, Energy};
use rand::rngs::SmallRng;
use crate::coupling_energy::CouplingEnergy;

use crate::evolving_sequence::EvolvingSequence;

pub struct FlipOnePos {
    succ_rate: AcceptanceStatistics,
    rng : SmallRng
}

impl FlipOnePos {
    pub fn new() -> FlipOnePos {
        FlipOnePos{ succ_rate: Default::default(), rng: SmallRng::from_entropy() }
    }
}

impl Mover<EvolvingSequence, CouplingEnergy> for FlipOnePos {

    fn perturb(&mut self, system: &mut EvolvingSequence, energy: &CouplingEnergy, acc: &mut dyn AcceptanceCriterion) -> Option<Range<usize>> {

        let pos: usize = self.rng.gen_range(0..system.seq_len());
        let new_aa: u8 = self.rng.gen_range(0..system.aa_cnt()) as u8;

        let delta_en: f64 = system.delta_energy(pos, new_aa, energy);
        if acc.check(0.0, delta_en) {
            system.sequence[pos] = new_aa;
            // ---------- test the energy consistency - in debug build only
            #[cfg(debug_assertions)]
                {
                    let new_total = energy.energy(&system);
                    if f64::abs((new_total - system.total_energy) - delta_en) > 0.01 {
                        let str = format!("Inconsistent energy! Total {} -> \
                        {new_total} with delta = {}, local delta: {delta_en} \
                        after flipping at pos {}", system.total_energy, (new_total - system.total_energy), pos);
                        panic!("{}", str);
                    }
                }

            self.succ_rate.n_succ += 1;
            system.total_energy += delta_en;
            return Some(pos..pos+1);
        }
        self.succ_rate.n_failed += 1;

        return None;
    }

    fn acceptance_statistics(&self) -> AcceptanceStatistics { self.succ_rate.clone() }

    fn max_range(&self) -> f64 { todo!() }

    fn set_max_range(&mut self, new_val: f64) { todo!() }
}
