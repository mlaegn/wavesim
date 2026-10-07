//! `waveslice`: one line across a reef, seen side-on, with the water's full motion.
//!
//! The main model averages the water over its depth, so it knows where a wave breaks and how
//! tall it is, but not its shape: a depth-averaged surface cannot fold over itself. Here the flow
//! in a vertical slice is solved as potential flow (inviscid, irrotational) with the exact
//! conditions at the free surface, which follows the surface as it steepens, overturns and
//! throws a lip, up to the moment the lip lands, after which the water is foam and splash and
//! this way of computing it ends (Longuet-Higgins & Cokelet 1976, Grilli et al. 1989).
//!
//! The flow is found by a boundary element method: only the edge of the water is discretised,
//! the free surface on top, the bed below and walls at the ends, and the potential inside follows
//! from the potential and its normal derivative on that edge. See [`bem`].
//!
//! Like `wavecore` it is pure numerics with no dependencies, and it is a separate step that a
//! normal run never uses.

pub mod bem;
mod linalg;

pub use bem::{Kind, Side, Solution, solve};
