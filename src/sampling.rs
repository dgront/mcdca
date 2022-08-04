use rand::Rng;
use rand::seq::SliceRandom;

use crate::EvolvingSequence;
use crate::observers::{Observer, ObserversSet};


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
    type S = EvolvingSequence;

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
    pub observers: ObserversSet<S>,
    energies: Vec<f32>,
    sweep_order: Vec<usize>,
}

impl SweepAllAA {

    pub fn new(seq_len:usize, n_aa:usize) -> SweepAllAA {

        let observers: ObserversSet<S> = ObserversSet::new();
        let energies: Vec<f32> = vec![0.0; n_aa];
        let mut sweep_order: Vec<usize> =  vec![0; seq_len];
        for i in 0..seq_len { sweep_order[i] = i }
        SweepAllAA { observers, energies, sweep_order }
    }
}

impl MCSweep for SweepAllAA {
    type S = EvolvingSequence;

    fn sweep(&mut self, system: &mut Self::S) {
        let mut rng = rand::thread_rng();

        self.sweep_order.shuffle(&mut rng);

        for pos in &self.sweep_order {
            // ---------- energy value for each possible move, including self-aa
            system.energy_row(*pos, &mut self.energies);
            // ---------- Boltzmann weights
            for i in 0..self.energies.len() { self.energies[i] = (-self.energies[i]).exp(); }

            let old_aa: u8 = system.sequence[*pos];
            let new_aa: u8 = rng.gen_range(0..system.aa_cnt()) as u8;

            let delta_en: f32 = self.energies[new_aa as usize] / self.energies[old_aa as usize];
            if delta_en > 0.0 && (-delta_en).exp() < rng.gen_range(0.0..1.0) { continue }
            system.sequence[*pos] = new_aa;
            system.total_energy += delta_en;
        }
    }
}