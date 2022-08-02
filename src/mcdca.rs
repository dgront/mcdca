use std::fs::File;
use std::io::Write;
use clap::{Parser};

use bioshell_core::sequence::{from_fasta_file, a3m_to_fasta, A3mConversionMode};

mod couplings;
mod sampling;
mod observers;
use crate::couplings::{Couplings, EvolvingSequence, counts_from_msa};
use crate::sampling::{SweepSingleAA, SimpleMCSampler};
use crate::observers::{EnergyHistogram, ObservedCounts, Observer, SequenceCollection};


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
    inner: i32,
    /// number of outer MC cycles
    #[clap(short, long, default_value_t = 100)]
    outer: i32,
    /// number of optimization cycles
    #[clap(short, long, default_value_t = 10, short='c')]
    optcycles: i32,
    /// input is RNA rather than a protein
    #[clap( long)]
    rna: bool,
}

pub fn update_couplings(target_counts: &Couplings, observed_counts: &Couplings, system: &mut EvolvingSequence) -> f32 {


    let dumping = 0.1;
    let nk = system.seq_len() * system.aa_cnt();
    let mut error: f32 = 0.0;
    for i_pos in 0..system.seq_len() {
        for i_aa in 0..system.aa_cnt() {
            let i_ind = i_pos* system.aa_cnt() + i_aa;
            for j_pos in 0..system.seq_len() {
                if i_pos == j_pos { continue }
                for j_aa in 0..system.aa_cnt() {
                    let j_ind = j_pos* system.aa_cnt() + j_aa;
                    let target_val = target_counts.data[i_ind][j_ind];
                    let observed_val = observed_counts.data[i_ind][j_ind];
                    let mut delta = target_val - observed_val;
                    if observed_val == 0.0 && target_val == 0.0 {
                        system.cplngs.data[i_ind][j_ind] = 0.0;
                        continue;
                    }
                    error += delta*delta;
                    delta /= 2.0 * observed_val;
                    delta *= dumping;
                    println!("{:3} {:2} {:3} {:2}  should be: {:5.3}  observed: {:5.3}, J: {:5.3} -> {:5.3} by {} err {}",
                             i_pos, i_aa, j_pos, j_aa,
                             target_val, observed_val,
                             system.cplngs.data[i_ind][j_ind], system.cplngs.data[i_ind][j_ind] - delta, delta, delta*delta);
                    system.cplngs.data[i_ind][j_ind] -= delta;
                }
            }
        }
    }

    return error;
}

pub fn main() {
    let args = Args::parse();

    // ---------- Input sequence
    let seq: String = from_fasta_file(&args.fasta)[0].to_string();

    // ---------- Read the target MSA and cleanup insertions
    let mut msa = from_fasta_file(&args.msa);
    a3m_to_fasta(&mut msa, &A3mConversionMode::RemoveSmallCaps);

    // ---------- Create sequence
    let alphabet: &str = if args.rna { "ACGU-" } else { "ACDEFGHIKLMNPQRSTVWY-" };
    let mut system: EvolvingSequence = EvolvingSequence::new(&seq, alphabet);

    let target = counts_from_msa(&system, &msa);
    let mut w = File::create("target_msa_counts.dat").unwrap();
    writeln!(&mut w, "{}", target).unwrap();

    // ---------- Create a MC sweep and a MC sampler
    let sweep: SweepSingleAA = SweepSingleAA {};
    let mut sampler = SimpleMCSampler::new();
    sampler.add_sweep(Box::new(sweep));

    // ---------- Observers
    let collect_seq = SequenceCollection::new("sequences", true);
    let energy_hist = EnergyHistogram::new(1.0,"en.dat");
    let coupled_pos = ObservedCounts::new(system.seq_len(), system.aa_cnt(), "observed_counts.dat");

    // sampler.inner_observers.push( Box::new(energy_hist));
    sampler.inner_observers.push( Box::new(coupled_pos));
    sampler.outer_observers.push( Box::new(collect_seq));

    // ---------- Other settings
    let n_inner: i32 = args.inner;
    let n_outer: i32 = args.outer;
    let n_cycles: i32 = args.optcycles;

    // ---------- Run the simulation!
    for _ in 0..n_cycles {
        // ---------- Forward step - infer counts
        sampler.run(&mut system, n_inner, n_outer);

        println!("target counts:\n{}", &target);
        let observed_counts : Option<&ObservedCounts> = sampler.get_observer("ObservedCounts");
        match observed_counts {
            None => {}
            Some(counts) => {
                let mut obs_copy = counts.get_counts().clone();
                obs_copy.normalize(counts.n_observed());
                let err = update_couplings(&target, &obs_copy, &mut system);
                println!("observed counts:\n{}", &counts.get_counts());
                println!("Normalized observed:\n{}",obs_copy);
                println!("couplings after update:\n{}", &system.cplngs);
                println!("ERROR: {}", err);
            }
        };

        // ---------- Write observations
        sampler.flush_observers();
    }
}