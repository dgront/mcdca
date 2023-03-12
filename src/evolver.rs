mod args;

use std::env;
use std::fs::File;
use std::io::Write;
use std::time::Instant;

use clap::Parser;

#[macro_use]
extern crate log;

use bioshell_core::sequence::{from_file, from_fasta_reader, a3m_to_fasta, A3mConversionMode,
                              SequenceProfile, ResidueTypeOrder};
use bioshell_montecarlo::{AcceptanceStatistics, IsothermalMC, Sampler};
use bioshell_sim::{Energy, Observer};

use mcdca_lib::CouplingEnergy;
use mcdca_lib::{Couplings, counts_from_msa, update_couplings};
use mcdca_lib::Pseudocounts;
use mcdca_lib::EvolvingSequence;
use mcdca_lib::{FlipOnePos};
use mcdca_lib::{ObservedCounts, SequenceCollection, PrintSequence};

use args::Args;

// make it multi-core

pub fn main() {

    let args = Args::parse();
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "debug") }
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
    println!("{}", &prof);
    let pseudo_cnts = Pseudocounts::new(pseudo_fraction, prof);
    let mut energy = CouplingEnergy{cplngs: target_counts.clone() };
    system.total_energy = energy.energy(&system);

    // ---------- Create an MC sampler and plug a mover into it
    let mut sampler: IsothermalMC<EvolvingSequence, CouplingEnergy> = IsothermalMC::new(1.0);
    sampler.add_mover(Box::new(FlipOnePos::new()));

    let mut obs_seq = PrintSequence{};

    // ---------- Run the simulation!
    let start = Instant::now();
    let mut recent_acceptance = AcceptanceStatistics::default();
    for i_newt in 0..n_cycles {
        // ---------- Observers

        let mut collect_seq = SequenceCollection::new(format!("sequences-{}.dat", i_newt).as_str(), false);
        let mut obs_freq = ObservedCounts::new(&system, "observed_counts.dat");

        for i in 0..args.outer {
            let stats = sampler.get_mover(0).acceptance_statistics();
            sampler.make_sweeps(args.inner, &mut system, &energy);
            let f_succ = stats.recent_success_rate(&recent_acceptance);
            recent_acceptance = stats;
            // obs_seq.observe(&system);
            collect_seq.observe(&system);
            obs_freq.observe(&system);
            // println!("{:6} {:9.3} {:5.3} {:.2?}", i, energy.energy(&system) as f64, f_succ, start.elapsed());
        }
        collect_seq.flush();
        // ---------- Newton optimisation step
        debug!("target counts:\n{}", &target_counts);
        debug!("observed counts:\n{}", &obs_freq.get_counts());
        obs_freq.normalize();
        let observed_counts = obs_freq.get_counts();
        debug!("Normalized observed:\n{}",observed_counts);
        let err = update_couplings(&target_counts, &observed_counts,
                                   &pseudo_cnts, newton_step, &mut energy.cplngs);
        info!("ERROR: {}", err);
        // ---------- recalculate energy - the old value became obsolete when couplings were changed
        system.total_energy = energy.energy(&system);
    }
}
