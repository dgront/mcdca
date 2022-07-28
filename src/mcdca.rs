

use bioshell_core::Sequence;

mod couplings;
mod sampling;
mod observers;
use crate::couplings::{Couplings, EvolvingSequence};
use crate::sampling::{SweepSingleAA, SimpleMCSampler};
use crate::observers::{EnergyHistogram, ObservedCounts, PrintSequence};


pub fn main() {
    let seq: String = String::from("MTYKLILNGKTLKGETTTEAVDAATAEKVFKQYANDNGVDGEWTYDDATKTFTVTE");
    let mut system: EvolvingSequence = EvolvingSequence::new(&seq, "ACDEFGHIKLMNPQRSTVWY-");

    let sweep: SweepSingleAA = SweepSingleAA {};

    // ---------- Observers
    let print_seq = PrintSequence { };
    let energy_hist = EnergyHistogram::new(1.0,"en.dat");
    let coupled_pos = ObservedCounts::new(system.seq_len(), system.aa_cnt(), "couplings.dat");

    let mut sampler = SimpleMCSampler::new();
    sampler.add_sweep(Box::new(sweep));
    sampler.inner_observers.push( Box::new(energy_hist));
    sampler.inner_observers.push( Box::new(coupled_pos));
    sampler.outer_observers.push( Box::new(print_seq));

    sampler.run(&mut system, 10, 10);

    sampler.close_observers();
}