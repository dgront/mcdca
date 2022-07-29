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
    output: String,
    n_observ: f32,
    tmp_i: Vec<usize>,
}

impl ObservedCounts {

    pub fn new(seq_len: usize, n_aa_types: usize, out_name:&str) -> ObservedCounts {
        let mut tmp_i: Vec<usize> = vec![0; seq_len];
        ObservedCounts {
            counts: Couplings::new(seq_len, n_aa_types),
            output: out_name.parse().unwrap(),
            n_observ: 0.0,
            tmp_i
        }
    }
}

impl Observer for ObservedCounts {
    type O = EvolvingSequence;

    fn observe(&mut self, obj: &Self::O) {

        self.n_observ += 1.0;
        for idx in 0..self.counts.n {
            self.tmp_i[idx] = idx * self.counts.k + obj.sequence[idx] as usize;
        }
        for i in 0..self.counts.n {
            let ii = self.tmp_i[i];
            for j in 0..self.counts.n {
                self.counts.data[ii][self.tmp_i[j]] += 1.0;
            }
        }
    }

    fn close(&mut self) {
        let mut out_writer = out_writer(&self.output.as_str());
        if self.n_observ > 0.0 { self.counts.normalize(self.n_observ); }
        out_writer.write(format!("{}", self.counts).as_bytes()).ok();
    }
}

/// Helper struct to store a single sequence data
struct SequenceEntry {
    energy: f32,
    sequence: String
}

/// Collects sequences together with their energy
///
/// The collection can estimate the observed counts from sequences it keeps assuming canonical ensemble
pub struct SequenceCollection {
    output: String,
    sequences: Vec<SequenceEntry>
}

impl SequenceCollection {

    /// Creates a new empty collection of sequence entries
    pub fn new(out_name:&str) -> SequenceCollection {
        SequenceCollection{output: out_name.parse().unwrap(), sequences: Vec::new()}
    }
}

impl Observer for SequenceCollection {
    type O = EvolvingSequence;

    /// Copy the current sequence and energy from an observed `EvolvingSequence` instance
    fn observe(&mut self, obj: &Self::O) {
        self.sequences.push(SequenceEntry{energy:obj.energy(), sequence:obj.decode_sequence()})
    }

    /// Stores observations in a file, or prints on the screen
    fn close(&mut self) {
        let mut out_writer = out_writer(&self.output.as_str());
        for s in self.sequences.iter() {
            out_writer.write(format!("{:3} {}\n", s.energy, s.sequence).as_bytes()).ok();
        }
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