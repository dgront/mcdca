use std::env;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::Write;
use std::time::Instant;
use clap::{Parser};

#[macro_use]
extern crate log;

use bioshell_core::chemical::ResidueType;
use bioshell_core::sequence::{from_fasta_file, a3m_to_fasta, A3mConversionMode, Sequence,
                              SequenceProfile, ResidueTypeOrder};
use bioshell_montecarlo::{AcceptanceStatistics, IsothermalMC, Sampler};
use bioshell_sim::{Energy, Observer};

use crate::coupling_energy::CouplingEnergy;

mod couplings;
mod sampling;
mod observers;
mod evolving_sequence;
mod coupling_energy;
mod pseudocounts;

use crate::couplings::{Couplings, counts_from_msa, update_couplings};
use crate::pseudocounts::Pseudocounts;
use crate::evolving_sequence::EvolvingSequence;
use crate::sampling::{FlipOnePos};
use crate::observers::{EnergyHistogram, ObservedCounts, SequenceCollection, PrintSequence};

// make it multi-core

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
/// Find direct coupling between amino acids in a given MSA
struct Args {
    /// input sequence in FASTA format
    #[clap(short, long, short='f')]
    fasta: String,
    /// input MSA in A3M format
    #[clap(short, long, short='m')]
    msa: String,
    /// number of inner MC cycles
    #[clap(short, long, default_value_t = 100)]
    inner: usize,
    /// number of outer MC cycles
    #[clap(short, long, default_value_t = 100)]
    outer: usize,
    /// number of optimization cycles
    #[clap(short, long, default_value_t = 10, short='c')]
    optcycles: u32,
    /// input is RNA rather than a protein
    #[clap( long)]
    rna: bool,
    /// number of optimization cycles
    #[clap(short, long, default_value_t = 0.01, short='n')]
    newton_step: f64,
    /// fraction of pseudocounts added to both observed and target statistics
    #[clap(short, long, default_value_t = 0.001, short='p')]
    pseudo_fraction: f64,
}


pub fn main() {

    let args = Args::parse();
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "debug") }
    env_logger::init();

    // ---------- Input sequence
    let seq: String = from_fasta_file(&args.fasta)[0].to_string();

    // ---------- Read the target MSA and cleanup insertions
    let mut msa = from_fasta_file(&args.msa);
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
