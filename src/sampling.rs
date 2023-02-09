use std::ops::Range;
use rand::{Rng, SeedableRng};

use bioshell_montecarlo::{AcceptanceCriterion, AcceptanceStatistics, Mover};
use bioshell_sim::Energy;
use rand::rngs::SmallRng;
use crate::coupling_energy::CouplingEnergy;

use crate::observers::{Observer, ObserversSet};
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


/*

/// Can perform a single Monte Carlo sweep on a given system
/// A type that implement `MCSweep` must be able to propose a MC move and evaluate energy
pub trait MCSweep {
    type S;
    fn sweep(&mut self, system: &mut Self::S);
}


/// Runs one or more sweeps in sequence, making some observations
pub struct SimpleMCSampler<S: 'static> {
    pub sweeps: Vec<Box<dyn MCSweep<S = S>>>,
    pub observers: ObserversSet<S>,
}

impl<S> SimpleMCSampler<S> {

    /// creates new sampler
    pub fn new() -> SimpleMCSampler<S> {
        let obs: ObserversSet<S> = ObserversSet::new();
        SimpleMCSampler {sweeps: Vec::new(), observers: obs }
    }

    pub fn run(&mut self, system: &mut S, inner_cycles: u32, outer_cycles: u32) {
        for _io in 0..outer_cycles {
            for _ii in 0..inner_cycles {
                for i in 0..self.sweeps.len() {
                    self.sweeps[i].sweep(system);
                }
                self.observers.observe(system);
            }
        }
    }

    /// Add another sweep to this Monte Carlo sampler
    pub fn add_sweep(&mut self, sweep: Box<dyn MCSweep<S = S>>) { self.sweeps.push(sweep); }
}

// ---------- DCA-related stuff
pub struct SweepSingleAA;

impl MCSweep for SweepSingleAA {
    type S = EvolvingSequence2;

    fn sweep(&mut self, system: &mut Self::S) {
        let mut rng = rand::thread_rng();
        for _ in 0..system.seq_len() {
            let pos: usize = rng.gen_range(0..system.seq_len());
            let new_aa: u8 = rng.gen_range(0..system.aa_cnt()) as u8;
            let delta_en: f32 = system.delta_energy(pos, new_aa);
            if delta_en > 0.0 && (-delta_en).exp() < rng.gen_range(0.0..1.0) { continue }
            system.sequence[pos] = new_aa;
            system.total_energy += delta_en;
        }
    }
}

pub struct SweepAllAA {
    energies: Vec<f32>,
    w: Vec<f32>,
    sweep_order: Vec<usize>,
}

impl SweepAllAA {

    pub fn new(seq_len:usize, n_aa:usize) -> SweepAllAA {

        let energies: Vec<f32> = vec![0.0; n_aa];
        let w: Vec<f32> = vec![0.0; n_aa];
        let mut sweep_order: Vec<usize> =  vec![0; seq_len];
        for i in 0..seq_len { sweep_order[i] = i }
        SweepAllAA { energies, w, sweep_order }
    }
}

impl MCSweep for SweepAllAA {
    type S = EvolvingSequence2;

    fn sweep(&mut self, system: &mut Self::S) {
        let mut rng = rand::thread_rng();

        self.sweep_order.shuffle(&mut rng);

        let mut tmp_i: Vec<usize> = vec![0; system.seq_len()];
        for idx in 0..system.seq_len() {
            tmp_i[idx] = idx * system.aa_cnt() + system.sequence[idx] as usize;
        }

        let n = system.seq_len();
        let k = system.aa_cnt();
        for pos in &self.sweep_order {
            // ---------- energy value for each possible move, including self-aa
            system.energy_row(*pos, &mut self.energies);
            // ---------- Boltzmann weights
            for i in 0..self.energies.len() { self.w[i] = (-self.energies[i]).exp(); }
            let total_w: f32 = self.w.iter().sum();
            // ---------- Normalize the Boltzmann factors
            self.w.iter_mut().for_each(|iel| *iel /= total_w);

            // ---------- observe counts
            for i in 0..n {             // --- for each "other" position in an chain sequence (MSA column)
                let ii = tmp_i[i];      // --- get letter at that "other" pos
                for aa_j in 0..k {      // --- for index of each letter that could happen in the mutated position
                    system.observed.data[ii][pos * k + aa_j] += self.w[aa_j];
                    system.observed.data[pos * k + aa_j][ii] += self.w[aa_j];
                }
            }

            // ---------- Symmetric acceptance criterion
            let old_aa: u8 = system.sequence[*pos];
            let rnd: f32 = rng.gen_range(0.0..1.0);
            let mut new_aa: u8 = (system.aa_cnt() - 1) as u8;   // Just in case the loop below would not select the last letter
            let mut sum: f32 = 0.0;
            for i in 0..system.aa_cnt() {
                sum += self.w[i];
                if rnd < sum {
                    new_aa = i as u8;
                    break;
                }
            }

            // ---------- Metropolis criterion on Boltzmann weights
            #[cfg(debug_assertions)]
                {
                    let delta_w: f32 = self.w[new_aa as usize] / self.w[old_aa as usize];
                    // ---------- test whether delta-energy is computed correctly
                    let delta_en = self.energies[new_aa as usize] - self.energies[old_aa as usize];
                    let delta_en_ctrl: f32 = system.delta_energy(*pos, new_aa);
                    if (delta_en_ctrl-delta_en).abs() > 0.0001 {
                        panic!("Incorrect delta-energy: {} vs {}, err is: {}", delta_en, delta_en_ctrl, (delta_en - delta_en_ctrl).abs());
                    }
                    // ---------- test whether Boltzmann factor for acceptance criterion is OK
                    if (delta_w-(-delta_en_ctrl).exp()).abs() > 0.0001 {
                        panic!("Incorrect boltzmann factor!");
                    }
                }
            // --- It wouldn't be necessary to update stuff if we hadn't changed amino acid
            if new_aa == old_aa { continue; }

            system.sequence[*pos] = new_aa;
            tmp_i[*pos] = (*pos) * system.aa_cnt() + new_aa as usize;
            system.total_energy += self.energies[new_aa as usize] - self.energies[old_aa as usize];
        }
    }
}

 */