use rand::Rng;

use crate::EvolvingSequence;
use crate::observers::Observer;


/// Can perform a single Monte Carlo sweep on a given system
/// A type that implement `MCSweep` must be able to propose a MC move and evaluate energy
pub trait MCSweep {
    type S;
    fn sweep(&mut self, system: &mut Self::S);
}

/// Runs one or more sweeps in sequence, making some observations
pub struct SimpleMCSampler<S> {
    pub sweeps: Vec<Box<dyn MCSweep<S = S>>>,
    pub inner_observers: Vec<Box<dyn Observer<O = S>>>,
    pub outer_observers: Vec<Box<dyn Observer<O = S>>>,
}

impl<S> SimpleMCSampler<S> {

    pub fn new() -> SimpleMCSampler<S> {
        SimpleMCSampler {sweeps: Vec::new(), inner_observers:Vec::new(), outer_observers:Vec::new() }
    }

    pub fn run(&mut self, system: &mut S, inner_cycles: i32, outer_cycles: i32) {
        for _io in 0..outer_cycles {
            for _ii in 0..inner_cycles {
                for i in 0..self.sweeps.len() {
                    self.sweeps[i].sweep(system);
                }
                for i in 0..self.inner_observers.len() {
                    self.inner_observers[i].observe(system);
                }
            }
            for i in 0..self.outer_observers.len() {
                self.outer_observers[i].observe(system);
            }
        }
    }

    pub fn get_observer<T: 'static>(&self, name: &str) -> Option<&T> {

        for o in self.inner_observers.iter() {
            if name == o.name() { return o.as_any().downcast_ref::<T>(); }
        }
        for o in self.outer_observers.iter() {
            if name == o.name() { return o.as_any().downcast_ref::<T>(); }
        }

        return None;
    }

    /// Add another sweep to this Monte Carlo sampler
    pub fn add_sweep(&mut self, sweep: Box<dyn MCSweep<S = S>>) { self.sweeps.push(sweep); }

    /// Call `close()` method for all observers this sampler posses
    /// This typically closes all opened files, computes statistics etc.
    pub fn close_observers(&mut self) {
        for o in self.inner_observers.iter_mut()  { o.close();}
        for o in self.outer_observers.iter_mut()  { o.close();}
    }

    /// Call `flush()` method for all observers this sampler posses
    /// This typically writes data to streams and clears buffers
    pub fn flush_observers(&mut self) {
        for o in self.inner_observers.iter_mut()  { o.flush();}
        for o in self.outer_observers.iter_mut()  { o.flush();}
    }
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
