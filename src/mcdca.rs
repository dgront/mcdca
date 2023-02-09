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
use bioshell_sim::Energy;

use crate::coupling_energy::CouplingEnergy;

mod couplings;
mod sampling;
mod observers;
mod evolving_sequence;
mod coupling_energy;

use crate::couplings::{Couplings, counts_from_msa};
use crate::evolving_sequence::EvolvingSequence;
use crate::sampling::{FlipOnePos};
use crate::observers::{EnergyHistogram, ObservedCounts, Observer, SequenceCollection, PrintSequence};

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


/// Update the `couplings` matrix to lower the distance between expected and target observations
///
/// ```math
/// j^{k,l}_{A,B} = j^{k,l}_{A,B} - \text{c} * \frac{p^{k,l}_{A,B} - \hat{p}^{k,l}_{A,B}}{p^{k,l}_{A,B}}
/// ```
pub fn update_couplings(target_counts: &Couplings, observed_counts: &Couplings, profile: &SequenceProfile,
                        pseudo_fraction: f64, newton_step: f64, couplings: &mut Couplings) -> f64 {


    let mut error: f64 = 0.0;
    let k_aa = profile.residue_order().size();
    let n_res = profile.len();
    for i_pos in 0..n_res {
        for i_aa in 0..k_aa {
            let i_ind = i_pos * k_aa + i_aa;
            for j_pos in 0..n_res {
                if i_pos == j_pos { continue }
                for j_aa in 0..k_aa {
                    let mut pseudo: f64 = profile.fraction(i_pos, i_aa) as f64 * profile.fraction(j_pos, j_aa) as f64;
                    pseudo *= pseudo_fraction;
                    let j_ind = j_pos * k_aa + j_aa;
                    let target_val = target_counts.data[i_ind][j_ind] as f64;
                    let observed_val = observed_counts.data[i_ind][j_ind] as f64;
                    let mut delta = target_val - observed_val;
                    if observed_val.abs() < 1e-7 && target_val.abs() < 1e-7 {
                        couplings.data[i_ind][j_ind] = 0.0;
                        continue;
                    }
                    error += delta*delta;
                    println!("{} {} {} {} : {:6.3} {:6.3} {:6.3} {:6.3} {}",
                             i_pos, j_pos, i_aa, j_aa, target_val, observed_val, delta, couplings.data[i_ind][j_ind], error);
                    delta = delta / (2.0 * observed_val + pseudo);
                    delta *= newton_step;
                        // if i_pos == 0 && j_pos == 1 && i_aa == 0 && j_aa == 0 {
                        //     println!("{} {} {} {} : {} {} {} ", i_pos, j_pos, i_aa, j_aa, target_val, observed_val, delta);
                        // }
                    if delta.abs() > 10000.0 {
                        println!("{:3} {:2} {:3} {:2}  should be: {:5.3}  observed: {:5.3}, J: {:5.3} -> {:5.3} by {} err {}",
                                 i_pos, i_aa, j_pos, j_aa,
                                 target_val, observed_val,
                                 couplings.data[i_ind][j_ind], couplings.data[i_ind][j_ind] as f64 - delta, delta, delta * delta);
                    }
                    couplings.data[i_ind][j_ind] -= delta as f32;
                }
            }
        }
    }
    println!("");
    return error;
}


pub fn main() {

    let args = Args::parse();
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
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
    let mut energy = CouplingEnergy{cplngs: target_counts.clone() };
    system.total_energy = energy.energy(&system);

    // ---------- Create a MC mover a MC sampler
    let mut sampler: IsothermalMC<EvolvingSequence, CouplingEnergy> = IsothermalMC::new(1.0);
    sampler.add_mover(Box::new(FlipOnePos::new()));

    let mut obs_seq = PrintSequence{};

    // ---------- Run the simulation!
    let start = Instant::now();
    let mut recent_acceptance = AcceptanceStatistics::default();
    for i_newt in 0..n_cycles {
        // ---------- Observers
        let mut collect_seq = SequenceCollection::new("sequences", true);
        let mut obs_freq = ObservedCounts::new(&system, "observed_counts.dat");

        for i in 0..args.outer {
            let stats = sampler.get_mover(0).acceptance_statistics();
            sampler.make_sweeps(args.inner, &mut system, &energy);
            let f_succ = stats.recent_success_rate(&recent_acceptance);
            recent_acceptance = stats;
            // obs_seq.observe(&system);
            // collect_seq.observe(&system);
            obs_freq.observe(&system);
            // println!("{:6} {:9.3} {:5.3} {:.2?}", i, energy.energy(&system) as f64, f_succ, start.elapsed());
        }
        // ---------- Newton optimisation step
        debug!("target counts:\n{}", &target_counts);
        debug!("observed counts:\n{}", &obs_freq.get_counts());
        obs_freq.normalize();
        let observed_counts = obs_freq.get_counts();
        debug!("Normalized observed:\n{}",observed_counts);
        // ---------- Observed counts `obs_copy` should be used in Newton step with SweepSingleAA:
        // let mut obs_copy = counts.get_counts().clone();
        // obs_copy.normalize(counts.n_observed() * system.seq_len() as f32);
        // ---------- When SweepAllAA is used, `system.observed` works better; they have been normalised 10 lines above
        // let mut obs_copy = system.observed.clone();
        let err = update_couplings(&target_counts, &observed_counts,
                                   &prof, pseudo_fraction, newton_step, &mut energy.cplngs);
        // debug!("couplings after update:\n{}", &energy.cplngs);
        info!("ERROR: {}", err);
        // ---------- recalculate energy - the old value became obsolete when couplings were changed
        system.total_energy = energy.energy(&system);
    }

}

/*
pub fn main() {
    let args = Args::parse();
    if env::var("RUST_LOG").is_err() { env::set_var("RUST_LOG", "info") }
    env_logger::init();

    // ---------- Input sequence
    let seq: String = from_fasta_file(&args.fasta)[0].to_string();

    // ---------- Read the target MSA and cleanup insertions
    let mut msa = from_fasta_file(&args.msa);
    a3m_to_fasta(&mut msa, &A3mConversionMode::RemoveSmallCaps);

    // ---------- Other settings
    let n_inner: u32 = args.inner;
    let n_outer: u32 = args.outer;
    let n_cycles: u32 = args.optcycles;
    let newton_step: f32 = args.newton_step;
    let pseudo_fraction: f32 = args.pseudo_fraction;

    // ---------- Create sequence
    let alphabet: &str = if args.rna { "ACGU-" } else { "ACDEFGHIKLMNPQRSTVWY-" };
    let mut system: EvolvingSequence2 = EvolvingSequence2::new(&seq, alphabet);

    // ---------- Create couplings and sequence profile
    let target = counts_from_msa(&system, &mut msa);
    let mut w = File::create("target_msa_counts.dat").unwrap();
    writeln!(&mut w, "{}", target).unwrap();
    let prof: SeqProfile = SeqProfile::new(&system, &msa);
    println!("{}", &prof);

    // ---------- Observers
    let collect_seq = SequenceCollection::new("sequences", true);
    // let energy_hist = EnergyHistogram::new(1.0,"en.dat");
    let coupled_pos = ObservedCounts::new(system.seq_len(), system.aa_cnt(), "observed_counts.dat");

    // ---------- Create a MC sweep and a MC sampler
    let mut sampler = SimpleMCSampler::new();
    // let sweep: SweepSingleAA = SweepSingleAA {};
    let sweep: SweepAllAA= SweepAllAA::new(seq.len(), alphabet.len());
    sampler.add_sweep(Box::new(sweep));

    // sampler.inner_observers.push( Box::new(energy_hist));
    sampler.observers.add_observer( Box::new(coupled_pos), 1);
    sampler.observers.add_observer( Box::new(collect_seq), n_inner);


    // ---------- Run the simulation!
    for _ in 0..n_cycles {
        // ---------- Forward step - infer counts
        sampler.run(&mut system, n_inner, n_outer);
        system.observed.normalize((n_inner * n_outer * system.seq_len() as u32 * 2) as f32);
        println!("{}", system.observed);

        // ---------- Newton optimisation step
        debug!("target counts:\n{}", &target);
        let observed_counts : Option<&ObservedCounts> = sampler.observers.get_observer("ObservedCounts");
        match observed_counts {
            None => {}
            Some(counts) => {
                // ---------- Observed counts `obs_copy` should be used in Newton step with SweepSingleAA:
                // let mut obs_copy = counts.get_counts().clone();
                // obs_copy.normalize(counts.n_observed() * system.seq_len() as f32);
                // ---------- When SweepAllAA is used, `system.observed` works better; they have been normalised 10 lines above
                let mut obs_copy = system.observed.clone();
                let err = update_couplings(&target, &obs_copy,
                                           &prof, pseudo_fraction, newton_step, &mut system);
                debug!("observed counts:\n{}", &counts.get_counts());
                debug!("Normalized observed:\n{}",obs_copy);
                debug!("couplings after update:\n{}", &system.cplngs);
                info!("ERROR: {}", err);
            }
        };
        system.observed.clear();

        // ---------- Write observations
        sampler.observers.flush_observers();
    }
}
*/