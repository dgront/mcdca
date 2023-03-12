mod args;

use std::env;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

use clap::Parser;

#[macro_use]
extern crate log;

use bioshell_core::chemical::ResidueType;
use bioshell_core::sequence::{from_fasta_reader, from_file, a3m_to_fasta, A3mConversionMode, Sequence,
                              SequenceProfile, ResidueTypeOrder};
use bioshell_montecarlo::{PERM, StepwiseBuilder};
use bioshell_sim::{Energy, Observer, ResizableSystem, System};

// use mcdca_lib::{counts_from_weighted_sequences};
use mcdca_lib::{counts_from_weighted_sequences, CouplingEnergy};
use mcdca_lib::{Couplings, counts_from_msa, update_couplings, EvolvingSequence, Pseudocounts};
use mcdca_lib::{EnergyHistogram, ObservedCounts, SequenceCollection, PrintSequence};
use mcdca_lib::{SequenceBuilder};

use args::Args;

// make it multi-core

pub fn main() {

    let args = Args::parse();
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    // ---------- Input sequence
    let seq: String = from_file(&args.fasta, from_fasta_reader)[0].to_string();

    // ---------- Read the target MSA and cleanup insertions
    let mut msa = from_file(&args.msa, from_fasta_reader);
    a3m_to_fasta(&mut msa, &A3mConversionMode::RemoveSmallCaps);

    // ---------- Other settings
    let n_cycles: u32 = args.optcycles;
    let newton_step: f64 = args.newton_step;
    let pseudo_fraction: f64 = args.pseudo_fraction;

    // ---------- Create sequence
    let alphabet: &str = if args.rna { "ACGU-" } else { "ACDEFGHIKLMNPQRSTVWY-" };
    let aa_order = ResidueTypeOrder::new(alphabet);
    let mut system: EvolvingSequence = EvolvingSequence::new(&seq, aa_order.clone());

    // ---------- Create couplings and sequence profile
    let target_counts = counts_from_msa(&system, &mut msa);
    let mut w = File::create("target_msa_counts.dat").unwrap();
    writeln!(&mut w, "{}", target_counts).unwrap();
    let prof: SequenceProfile = SequenceProfile::new(aa_order.clone(), &msa);
    // println!("{}", &prof);
    let pseudo_cnts = Pseudocounts::new(pseudo_fraction, prof);
    let mut energy = CouplingEnergy{cplngs: target_counts.clone() };
    system.total_energy = energy.energy(&system);

    let mut obs_seq = PrintSequence{};

    // ---------- Run the simulation!
    let start = Instant::now();
    for i_newt in 0..n_cycles {
        // ---------- Create a PERM sampler
        let builder = Box::new(SequenceBuilder::new(1.0,system.aa_cnt()));
        let mut sampler: PERM<EvolvingSequence, CouplingEnergy> = PERM::new(system.seq_len(), 0.1, 10.0,
                                                                            builder);
        sampler.set_w_scale(20.0);

        // ---------- Observers
        let mut collect_seq = SequenceCollection::new(format!("sequences-{}.dat", i_newt).as_str(), false);

        let mut i_chains: usize = 0;
        while i_chains < args.outer {
            // --- build a new chain
            system.weight = sampler.build(&mut system, &energy);
            // --- skip system if it was pruned
            if system.size() < system.capacity() { continue; }
            system.total_energy = energy.energy(&system);
            // println!("{}   {} {} {:?}", system.size(), system.total_energy, system.weight, sampler.weights(system.size() - 1));
            // obs_seq.observe(&system);
            collect_seq.observe(&system);
            i_chains += 1;
        }

        collect_seq.normalize_sequence_weights();
        let observed_counts = counts_from_weighted_sequences(&system, collect_seq.iter());
        collect_seq.flush();
        // debug!("Normalized observed:\n{}",observed_counts);
        // ---------- Newton optimisation step
        debug!("target counts:\n{}", &target_counts);
        debug!("observed counts:\n{}", &observed_counts);
        let err = update_couplings(&target_counts, &observed_counts,
                                   &pseudo_cnts, newton_step, &mut energy.cplngs);
        info!("ERROR: {}", err);
        // ---------- recalculate energy - the old value became obsolete when couplings were changed
        system.total_energy = energy.energy(&system);
    }
}
