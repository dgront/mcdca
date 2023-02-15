use clap::{Parser};


/// Command line arguments for the `evolver` and `builder` applications
///
/// Since both apps share most of their flags, this struct has been moved into a separate file
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
pub(crate) struct Args {
    /// input sequence in FASTA format
    #[clap(short, long, short='f')]
    pub fasta: String,
    /// input MSA in A3M format
    #[clap(short, long, short='m')]
    pub msa: String,
    /// number of inner MC cycles
    #[clap(short, long, default_value_t = 100)]
    pub inner: usize,
    /// number of outer MC cycles
    #[clap(short, long, default_value_t = 100)]
    pub outer: usize,
    /// number of optimization cycles
    #[clap(short, long, default_value_t = 10, short='c')]
    pub optcycles: u32,
    /// input is RNA rather than a protein
    #[clap( long)]
    pub rna: bool,
    /// number of optimization cycles
    #[clap(short, long, default_value_t = 0.01, short='n')]
    pub newton_step: f64,
    /// fraction of pseudocounts added to both observed and target statistics
    #[clap(short, long, default_value_t = 0.001, short='p')]
    pub pseudo_fraction: f64,
}