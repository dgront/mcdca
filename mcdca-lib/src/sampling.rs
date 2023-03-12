use std::ops::Range;
use rand::{Rng, SeedableRng};

use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover, StepwiseMover};
use bioshell_sim::{Observer, ObserversSet, Energy, ResizableSystem, System};
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

/// Builds a protein / nucleic sequence
pub struct SequenceBuilder {
    /// temperature for the simulation
    pub temperature: f64,
    rng : SmallRng,
    energy_by_letter: Vec<f32>  // scratch memory to store energy evaluations
}

impl SequenceBuilder {
    /// Create a new builder that construct sequences from Boltzmann distribution at given temperature
    pub fn new(temperature: f64, alphabet_size: usize) -> SequenceBuilder {
        SequenceBuilder { temperature, rng: SmallRng::from_entropy(), energy_by_letter: vec![0.0; alphabet_size] }
    }
}

impl StepwiseMover<EvolvingSequence, CouplingEnergy> for SequenceBuilder {

    /// Starts a new sequence from a random letter
    fn start(&mut self, system: &mut EvolvingSequence, _energy: &CouplingEnergy) -> f64 {
        system.set_size(1);
        system.sequence[0] = self.rng.gen_range(0..system.aa_cnt()) as u8;

        return 1.0;
    }

    /// Grows the current chain by one bead
    fn grow_by_one(&mut self, system: &mut EvolvingSequence, energy: &CouplingEnergy) -> f64 {

        let i = system.size();
        system.set_size(i + 1);

        system.energy_by_letter(i, energy, &mut self.energy_by_letter);
        let mut total: f32 = 0.0;
        for i in 0..self.energy_by_letter.len() {
            self.energy_by_letter[i] = (-self.energy_by_letter[i]/self.temperature as f32).exp();
            total += self.energy_by_letter[i];
        }
        if total < 1e-20 || total.is_infinite() { return 0.0; }       // --- no suitable move generated
        // ---------- select one of the possible extension by importance sampling
        let mut rng = rand::thread_rng();
        let r = rng.gen_range(0.0..total);
        let mut which_v: usize = 0;
        let mut s = self.energy_by_letter[which_v];
        while s <= r {
            which_v += 1;
            s += self.energy_by_letter[which_v]
        }

        // ---------- set the selected letter
        system.sequence[i] = which_v as u8;

        // ---------- return the statistical weight
        return total as f64;
    }
}