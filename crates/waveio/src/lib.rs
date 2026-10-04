//! `waveio`: the file formats around `wavecore`.
//!
//! A *bed file* is what the fetch tool writes and the solver reads; a *run* is what
//! the solver writes and a viewer reads. Both are a small JSON header next to raw
//! little-endian `f32` data, specified in `docs/data-format.md`.

mod bed;
mod error;
mod run;

pub use bed::{BED_FORMAT, Bed, BedHeader, Frame};
pub use error::Error;
pub use run::{RUN_FORMAT, Run, RunHeader, RunWriter, Stats};
