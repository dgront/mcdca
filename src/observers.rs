use std::os::macos::raw::stat;

use bioshell_numerical::statistics::Histogram;
use bioshell_core::utils::out_writer;

// ---------- generic traits and types that should be relocated to bioshell

use crate::{Couplings, EvolvingSequence};

pub trait Observer {
    type O;
    fn observe(&mut self, object: &Self::O);
    fn close(&mut self);
}

// ---------- DCA-related stuff
pub struct PrintSequence;

impl Observer for PrintSequence {
    type O = EvolvingSequence;

    fn observe(&mut self, obj: &Self::O) {
        println!("{:4} {}", obj.energy(), obj.decode_sequence());
    }

    fn close(&mut self) {}
}


pub struct ObservedCounts {
    counts: Couplings,
    output: String
}

impl ObservedCounts {

    pub fn new(seq_len: usize, n_aa_types: usize, out_name:&str) -> ObservedCounts {
        ObservedCounts{counts: Couplings::new(seq_len, n_aa_types), output: out_name.parse().unwrap()}
    }
}

impl Observer for ObservedCounts {
    type O = EvolvingSequence;

    fn observe(&mut self, obj: &Self::O) {
        let mut pos_i: usize = 0;
        for aa_i in obj.sequence.iter() {
            let mut pos_j: usize = 0;
            for aa_j in obj.sequence.iter() {
                self.counts.data[pos_i + *aa_i as usize][pos_j + *aa_j as usize] += 1.0;
                pos_j += self.counts.k;
            }
            pos_i += self.counts.k;
        }
    }

    fn close(&mut self) {
        let mut out_writer = out_writer(&self.output.as_str());
        out_writer.write(format!("{}", self.counts).as_bytes()).ok() ;
    }
}

/// Observer that collects energy observations into a histogram
pub struct EnergyHistogram {
    stats: Histogram,
    output: String
}

impl EnergyHistogram {

    /// Create a new observer and initialise its histogram with a given bin width.
    ///
    /// # Arguments
    ///
    /// * `energy_bin` - energy bin width
    /// * `out_name` - name of the output file; use `stdout` or an empty string to print on the screen
    ///
    /// # Examples
    ///
    /// ```
    /// use observers::EnergyHistogram;
    /// let en_obs = EnergyHistogram::new(1.0, "");
    /// ```
    pub fn new(energy_bin:f64, out_name:&str) -> EnergyHistogram {
        EnergyHistogram{stats: Histogram::by_bin_width(energy_bin), output: out_name.parse().unwrap() }
    }
}

impl Observer for EnergyHistogram {
    type O = EvolvingSequence;

    fn observe(&mut self, obj: &Self::O) {
        self.stats.insert(obj.total_energy as f64);
    }

    fn close(&mut self) {
        let mut out_writer = out_writer(&self.output.as_str());
        out_writer.write(format!("{}", self.stats).as_bytes()).ok() ;
    }
}