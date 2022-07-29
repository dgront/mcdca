use std::fs::File;
use std::io::Write;
use clap::{Parser};

use bioshell_core::sequence::{from_fasta_file, a3m_to_fasta, A3mConversionMode};

mod couplings;
mod sampling;
mod observers;
use crate::couplings::{Couplings, EvolvingSequence, counts_from_msa};
use crate::sampling::{SweepSingleAA, SimpleMCSampler};
use crate::observers::{EnergyHistogram, ObservedCounts, SequenceCollection};


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
    /// Inner cycles
    #[clap(short, long, default_value_t = 100)]
    inner: i32,
    /// outer cycles
    #[clap(short, long, default_value_t = 100)]
    outer: i32,
    /// input is RNA rather than a protein
    #[clap( long)]
    rna: bool,
}

pub fn update_couplings(target_counts: &Couplings, observed_counts: &Couplings, system: &mut EvolvingSequence) {

    let nk = system.seq_len()*system.aa_cnt();
    for i in 0..nk {
        for j in 0..nk {
            let mut delta = target_counts.data[i][j] - observed_counts.data[i][j];
            delta /= (2.0 * observed_counts.data[i][j]);
            system.cplngs.data[i][j] -= delta;
        }
    }
}

pub fn main() {
    let args = Args::parse();

    // ---------- Input sequence
    let seq: String = from_fasta_file(&args.fasta)[0].to_string();

    // ---------- Read the target MSA and cleanup insertions
    let mut msa = from_fasta_file(&args.msa);
    a3m_to_fasta(&mut msa, &A3mConversionMode::RemoveSmallCaps);

    // ---------- Create sequence
    let alphabet: &str = if args.rna { "ACGU" } else { "ACDEFGHIKLMNPQRSTVWY-" };
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

    // ---------- Run the simulation!
    for _ in 0..5 {
        // ---------- Forward step - infer counts
        sampler.run(&mut system, n_inner, n_outer);

        // ---------- Write observations
        sampler.close_observers();

        // update_couplings(&target, coupled_pos.get_counts(), &mut system);
    }
}