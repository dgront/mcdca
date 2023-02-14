use bioshell_core::sequence::{SequenceProfile};

pub struct Pseudocounts {
    pseudo_fraction: f64,
    profile:SequenceProfile
}

impl Pseudocounts {
    pub fn new(pseudo_fraction: f64, seq_profile:SequenceProfile) -> Pseudocounts {
        Pseudocounts{ pseudo_fraction, profile: seq_profile }
    }

    pub fn pseudo_fraction(&self,i_pos: usize, i_aa: usize, j_pos: usize, j_aa: usize) -> f64 {
        self.pseudo_fraction * self.profile.fraction(i_pos, i_aa) as f64 * self.profile.fraction(j_pos, j_aa) as f64
    }
}