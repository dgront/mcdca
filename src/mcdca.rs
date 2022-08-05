use std::env;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::Write;
use clap::{Parser};

#[macro_use]
extern crate log;

use bioshell_core::sequence::{from_fasta_file, a3m_to_fasta, A3mConversionMode, Sequence};

mod couplings;
mod sampling;
mod observers;
use crate::couplings::{Couplings, EvolvingSequence, counts_from_msa};
use crate::sampling::{SweepSingleAA, SimpleMCSampler};
use crate::observers::{EnergyHistogram, ObservedCounts, Observer, SequenceCollection};

// Run 2GB1 case or FDX
// make it multi-core
// implement recycked sampling

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
    inner: u32,
    /// number of outer MC cycles
    #[clap(short, long, default_value_t = 100)]
    outer: u32,
    /// number of optimization cycles
    #[clap(short, long, default_value_t = 10, short='c')]
    optcycles: u32,
    /// input is RNA rather than a protein
    #[clap( long)]
    rna: bool,
    /// number of optimization cycles
    #[clap(short, long, default_value_t = 0.01, short='n')]
    newton_step: f32,
}

pub fn update_couplings(target_counts: &Couplings, observed_counts: &Couplings, profile: &SeqProfile,
                        pseudo_fraction: f32, newton_step: f32, system: &mut EvolvingSequence) -> f32 {


    let mut error: f32 = 0.0;
    for i_pos in 0..system.seq_len() {
        for i_aa in 0..system.aa_cnt() {
            let i_ind = i_pos* system.aa_cnt() + i_aa;
            for j_pos in 0..system.seq_len() {
                if i_pos == j_pos { continue }
                for j_aa in 0..system.aa_cnt() {
                    let pseudo = profile.fraction(i_pos, i_aa) * profile.fraction(j_pos, j_aa) * pseudo_fraction;
                    let j_ind = j_pos* system.aa_cnt() + j_aa;
                    let target_val = target_counts.data[i_ind][j_ind];
                    let observed_val = observed_counts.data[i_ind][j_ind];
                    let mut delta = target_val - observed_val;
                    if observed_val.abs() < 1e-7 && target_val.abs() < 1e-7 {
                        system.cplngs.data[i_ind][j_ind] = 0.0;
                        continue;
                    }
                    error += delta*delta;
                    delta = delta / (2.0 * observed_val + pseudo);
                    delta *= newton_step;
                    if delta.abs() > 10000.0 {
                        println!("{:3} {:2} {:3} {:2}  should be: {:5.3}  observed: {:5.3}, J: {:5.3} -> {:5.3} by {} err {}",
                                 i_pos, i_aa, j_pos, j_aa,
                                 target_val, observed_val,
                                 system.cplngs.data[i_ind][j_ind], system.cplngs.data[i_ind][j_ind] - delta, delta, delta * delta);
                    }
                    system.cplngs.data[i_ind][j_ind] -= delta;
                }
            }
        }
    }

    return error;
}

pub struct SeqProfile {
    n: usize,
    data: Vec<Vec<f32>>,
}

impl SeqProfile {
    pub fn new(system: &EvolvingSequence, msa: &Vec<Sequence>) -> SeqProfile {

        let n = system.seq_len();
        let k = system.aa_cnt();
        let mut data: Vec<Vec<f32>> = vec![vec![0.0; k]; n];

        for sequence in msa {
            if sequence.len() != n {
                eprintln!("\nSequence of incorrect length! Is: {}, should be: {}. The sequence skipped::\n {}\n",
                          sequence.len(), n, sequence);
                continue;
            }            for idx in 0..n {
                let aa_idx = system.aa_to_index(&sequence.char(idx)) as usize;
                data[idx][aa_idx] += 1.0;
            }
        }

        for idx in 0..n {
            let sum: f32 = data[idx].iter().sum();
            for aa_idx in 0..k {
                data[idx][aa_idx] /= sum;
            }
        }

        SeqProfile{n, data}
    }

    pub fn fraction(&self, pos: usize, aa: usize) -> f32 { self.data[pos][aa] }

}

impl Display for SeqProfile {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (i, row) in self.data.iter().enumerate() {
            write!(f, "{i:4} ");
            for val in row {
                write!(f, "{val:5.3} ");
            }
            writeln!(f, "");
        }
        Ok(())
    }
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
    let n_inner: u32 = args.inner;
    let n_outer: u32 = args.outer;
    let n_cycles: u32 = args.optcycles;
    let newton_step: f32 = args.newton_step;

    // ---------- Create sequence
    let alphabet: &str = if args.rna { "ACGU-" } else { "ACDEFGHIKLMNPQRSTVWY-" };
    let mut system: EvolvingSequence = EvolvingSequence::new(&seq, alphabet);

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
    let sweep: SweepSingleAA = SweepSingleAA {};
    // let sweep: SweepAllAA= SweepAllAA::new(seq.len(), alphabet.len());
    sampler.add_sweep(Box::new(sweep));

    // sampler.inner_observers.push( Box::new(energy_hist));
    sampler.observers.add_observer( Box::new(coupled_pos), 1);
    sampler.observers.add_observer( Box::new(collect_seq), n_inner);


    // ---------- Run the simulation!
    for _ in 0..n_cycles {
        // ---------- Forward step - infer counts
        sampler.run(&mut system, n_inner, n_outer);
        system.observed.normalize((n_inner*n_outer) as f32 );
        println!("{}", system.observed);
        system.observed.clear();


        debug!("target counts:\n{}", &target);
        let observed_counts : Option<&ObservedCounts> = sampler.observers.get_observer("ObservedCounts");
        match observed_counts {
            None => {}
            Some(counts) => {
                let mut obs_copy = counts.get_counts().clone();
                obs_copy.normalize(counts.n_observed());
                let err = update_couplings(&target, &obs_copy,
                                           &prof, 0.001, newton_step, &mut system);
                debug!("observed counts:\n{}", &counts.get_counts());
                debug!("Normalized observed:\n{}",obs_copy);
                debug!("couplings after update:\n{}", &system.cplngs);
                debug!("ERROR: {}", err);
            }
        };

        // ---------- Write observations
        sampler.observers.flush_observers();
    }
}