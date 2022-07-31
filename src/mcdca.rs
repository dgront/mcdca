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

pub fn update_couplings(target_counts: &Couplings, observed_counts: &Couplings, n_observed: f32, system: &mut EvolvingSequence) -> f32 {

    let nk = system.seq_len()*system.aa_cnt();
    let mut error: f32 = 0.0;
    for i in 0..nk {
        let i_pos = i / system.aa_cnt();
        for j in 0..nk {
            let j_pos = j / system.aa_cnt();
            if i_pos == j_pos { continue }
            if observed_counts.data[i][j] == 0.0 {
                system.cplngs.data[i][j] = 0.0;
            } else {
                let obs_counts_normalised = observed_counts.data[i][j] / n_observed;
                let mut delta = target_counts.data[i][j] - obs_counts_normalised;
                error += delta;
                delta /= 2.0 * obs_counts_normalised;
                system.cplngs.data[i][j] -= delta;
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
    let collect_seq = SequenceCollection::new("sequences");
    let energy_hist = EnergyHistogram::new(1.0,"en.dat");
    let coupled_pos = ObservedCounts::new(system.seq_len(), system.aa_cnt(), "observed_counts.dat");

    sampler.inner_observers.push( Box::new(energy_hist));
    sampler.outer_observers.push( Box::new(coupled_pos));
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
                // counts.normalize();
                let err = update_couplings(&target, counts.get_counts(),counts.n_observed(), &mut system);
                println!("observed counts:\n{}", &counts.get_counts());
                println!("couplings after update:\n{}", &system.cplngs);
                println!("ERROR: {}", err);
            }
        };

        // ---------- Write observations
        sampler.flush_observers();
    }
}