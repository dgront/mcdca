use std::any::Any;

use bioshell_numerical::statistics::Histogram;
use bioshell_core::utils::{out_writer, writes_to_screen};

// ---------- generic traits and types that should be relocated to bioshell

use crate::{Couplings, EvolvingSequence};

/// Observer takes observations of a system of a generic type `O`
pub trait Observer {
    /// The type of objects being observed byt his observer
    type O;
    /// Takes observations
    fn observe(&mut self, object: &Self::O);
    fn flush(&mut self);
    fn close(&mut self);
    fn name(&self) -> &str;
    fn as_any(&self) -> &dyn Any;
}

/// A set of observers, that observe a system of a generic type `S`
///
///
pub struct  ObserversSet<S: 'static> {
    n_called: u32,
    observers: Vec<Box<dyn Observer<O = S>>>,
    lag_times : Vec<u32>
}

impl<S> ObserversSet<S> {

    pub fn new() -> ObserversSet<S> {
        ObserversSet {n_called: 0, observers:Vec::new(), lag_times:Vec::new() }
    }

    pub fn add_observer(&mut self, o: Box<dyn Observer<O=S>>, lag_time: u32) {
        self.observers.push(o);
        self.lag_times.push(lag_time);
    }

    pub fn get_observers(&self) -> &Vec<Box<dyn Observer<O=S>>> { &self.observers }

    pub fn observe(&mut self, object: &S) {
        for i in 0..self.observers.len() {
            if self.n_called % self.lag_times[i] == 0 {
                self.observers[i].observe(object);
            }
        }
        self.n_called += 1;
    }

    /// Call `flush()` method for all observers this sampler posses
    /// This typically writes data to streams and clears buffers
    pub fn flush_observers(&mut self) { for o in self.observers.iter_mut() { o.flush(); } }

    /// Call `close()` method for all observers this sampler posses
    /// This typically closes all opened files, computes statistics etc.
    pub fn close_observers(&mut self) { for o in self.observers.iter_mut() { o.close(); } }

    pub fn get_observer<T: 'static>(&self, name: &str) -> Option<&T> {
        for o in self.get_observers() {
            if name == o.name() { return o.as_any().downcast_ref::<T>(); }
        }

        return None;
    }
}

// ---------- DCA-related stuff
pub struct PrintSequence;

impl Observer for PrintSequence {
    type O = EvolvingSequence;

    fn observe(&mut self, obj: &Self::O) {
        println!("{:4} {}", obj.energy(), obj.decode_sequence());
    }

    fn flush(&mut self) {}

    fn close(&mut self) {}

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

    pub fn new(seq_len: usize, n_aa_types: usize, out_name:&str) -> ObservedCounts {
        let mut tmp_i: Vec<usize> = vec![0; seq_len];
        ObservedCounts {
            counts: Couplings::new(seq_len, n_aa_types),
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

    fn flush(&mut self) {
        let mut out_writer = out_writer(&self.output.as_str());
        if self.n_observ > 0.0 { self.counts.normalize(self.n_observ); }
        out_writer.write(format!("{}", self.counts).as_bytes()).ok();
        self.counts.clear();
        self.n_observ = 0.0;
    }

    fn close(&mut self) {}

    fn name(&self) -> &str { "ObservedCounts" }
    fn as_any(&self) -> &dyn Any { self }
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
}

impl Observer for SequenceCollection {
    type O = EvolvingSequence;

    /// Copy the current sequence and energy from an observed `EvolvingSequence` instance
    fn observe(&mut self, obj: &Self::O) {
        self.sequences.push(SequenceEntry{energy:obj.energy(), sequence:obj.decode_sequence()})
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
        let mut out_writer = out_writer(fname.as_str());
        for s in self.sequences.iter() {
            out_writer.write(format!("{:3} {}\n", s.energy, s.sequence).as_bytes()).ok();
        }
        self.sequences.clear();
    }

    fn close(&mut self) {}

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
    type O = EvolvingSequence;

    fn observe(&mut self, obj: &Self::O) {
        self.stats.insert(obj.total_energy as f64);
    }

    fn flush(&mut self) {
        let mut out_writer = out_writer(&self.output.as_str());
        out_writer.write(format!("{}", self.stats).as_bytes()).ok() ;
    }

    fn close(&mut self) {}

    fn name(&self) -> &str { "EnergyHistogram" }

    fn as_any(&self) -> &dyn Any { self }
}