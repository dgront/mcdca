use std::any::Any;
use std::slice::Iter;

use bioshell_sim::{Observer, ObserversSet};
use bioshell_statistics::Histogram;
use bioshell_core::utils::{out_writer, writes_to_screen};


use crate::{Couplings, EvolvingSequence, SequenceEntry};

/// Prints a sequence on a screen.
///
/// This simple observer just print energy of a sequence and the sequence itself in a single row
pub struct PrintSequence;

impl Observer for PrintSequence {
    type S = EvolvingSequence;

    fn observe(&mut self, obj: &Self::S) { println!("{} {}", obj.total_energy, obj.decode_sequence()); }

    fn flush(&mut self) {}

    fn name(&self) -> &str { "PrintSequence" }

    fn as_any(&self) -> &dyn Any { self }
}


pub struct ObservedCounts {
    counts: Couplings,
    output: String,
    n_observ: f32,
    tmp_i: Vec<usize>,
}

impl ObservedCounts {

    pub fn new(system: &EvolvingSequence, out_name:&str) -> ObservedCounts {
        let mut tmp_i: Vec<usize> = vec![0; system.seq_len()];
        ObservedCounts {
            counts: Couplings::new(system.seq_len(), system.aa_cnt()),
            output: out_name.parse().unwrap(),
            n_observ: 0.0,
            tmp_i
        }
    }

    /// Provide read-only access to the current state of observed counts
    pub fn get_counts(&self) -> &Couplings { &self.counts }

    /// Normalize counts
    pub fn normalize(&mut self) { self.counts.normalize(self.n_observ); }

    /// The number of sequences seen by this observer.
    /// This is the constant used to normalize the counts.
    pub fn n_observed(& self) -> f32 { self.n_observ }

}

impl Observer for ObservedCounts {
    type S = EvolvingSequence;

    fn observe(&mut self, obj: &Self::S) {

        self.n_observ += 1.0;
        for idx in 0..self.counts.seq_length {
            self.tmp_i[idx] = idx * self.counts.k_aa_types + obj.sequence[idx] as usize;
        }

        for i_pos in 1..self.counts.seq_length {
            for j_pos in 0..i_pos {
                self.counts.data[self.tmp_i[i_pos]][self.tmp_i[j_pos]] += 0.5;
                self.counts.data[self.tmp_i[j_pos]][self.tmp_i[i_pos]] += 0.5;
            }
        }
    }

    fn flush(&mut self) {
        let mut out_writer = out_writer(&self.output.as_str(),true);
        if self.n_observ > 0.0 { self.counts.normalize(self.n_observ); }
        out_writer.write(format!("{}", self.counts).as_bytes()).ok();
        self.counts.clear();
        self.n_observ = 0.0;
    }

    fn name(&self) -> &str { "ObservedCounts" }
    fn as_any(&self) -> &dyn Any { self }
}


/// Collects sequences.
///
///
pub struct SequenceCollection {
    flush_separately: bool,
    flush_id:i16,
    output: String,
    sequences: Vec<SequenceEntry>
}

impl SequenceCollection {

    /// Creates a new empty collection of sequence entries
    ///
    /// # Arguments
    ///
    /// * `out_name` -  name of the output file; use `stdout` or an empty string to print on the screen
    /// * `flush_separately` - each flush goes to a separate file, which will be named `out_name + i`
    ///     where `i` is the flush index and `out_name` is the original name given here. This flag
    ///     is not relevant when the sequences are flushed into a screen
    pub fn new(out_name:&str, flush_separately:bool) -> SequenceCollection {
        SequenceCollection{flush_separately, flush_id:0, output: out_name.parse().unwrap(), sequences: Vec::new()}
    }

    /// Provides iterator over [`SequenceEntry`] objects contained in this collection.
    ///
    /// Note that a [`SequenceCollection`] is cleared after each [`SequenceCollection::flush()`] call,
    ///so the returned iterator will provide only sequences gathered since the most recent  flush.
    pub fn iter(&self) -> Iter<'_, SequenceEntry> { self.sequences.iter() }

    /// Normalizes weights of sequences stored in this collection
    ///
    /// After this call, weight of all the sequences will sum up to 1.0
    pub fn normalize_sequence_weights(&mut self) {
        let mut total = 0.0;
        for sequence in &self.sequences {
            total += sequence.weight;
        }
        for sequence in &mut self.sequences {
            sequence.weight /= total;
        }
    }
}

impl Observer for SequenceCollection {
    type S = EvolvingSequence;

    /// Copy the current sequence and energy from an observed `EvolvingSequence` instance
    fn observe(&mut self, obj: &Self::S) {
        self.sequences.push(SequenceEntry{energy:obj.total_energy, weight: obj.weight, sequence:obj.decode_sequence()})
    }

    /// Stores observations in a file or prints on the screen.
    ///
    /// When the `flush_separately` flag was set to `true`, it will create a separate file for the current pool
    /// of sequences. Otherwise they will be appended to the existing one. The collection will be cleared after this call
    fn flush(&mut self) {

        let fname: String = if writes_to_screen(&self.output.as_str()) {
            self.output.to_string()         // --- need to_string() makes a copy as we can't move it
        } else {
            self.flush_id += 1;
            format!("{}-{}",&self.output, self.flush_id)
        };
        let mut out_writer = out_writer(fname.as_str(), true);
        for s in self.sequences.iter() {
            out_writer.write(format!("{:11.6} {:8.6} {}\n", s.energy, s.weight, s.sequence).as_bytes()).ok();
        }
        self.sequences.clear();
    }

    fn name(&self) -> &str { "SequenceCollection" }

    fn as_any(&self) -> &dyn Any { self }
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
    type S = EvolvingSequence;

    fn observe(&mut self, obj: &Self::S) {
        self.stats.insert(obj.total_energy as f64);
    }

    fn flush(&mut self) {
        let mut out_writer = out_writer(&self.output.as_str(), true);
        out_writer.write(format!("{}", self.stats).as_bytes()).ok() ;
    }

    fn name(&self) -> &str { "EnergyHistogram" }

    fn as_any(&self) -> &dyn Any { self }
}