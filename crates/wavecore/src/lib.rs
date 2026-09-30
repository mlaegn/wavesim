//! `wavecore`: pure-numerics nearshore wave solver.
//!
//! No file I/O, no graphics, no GDAL, no threads. That keeps it compilable to WASM
//! and portable to a GPU. Reading files and writing output belong in other crates.
//!
//! ```
//! use wavecore::{Bathymetry, Grid, Solver, State};
//!
//! let grid = Grid::new(64, 48, 3.0, 3.0);
//! let bed = Bathymetry::plane_beach(&grid, 0.02, 8.0);
//! let mut state = State::lake_at_rest(&grid, &bed, 0.0);
//! let solver = Solver::new(grid, bed);
//!
//! let dt = solver.stable_dt(&state);
//! solver.step(&mut state, dt);
//! ```

pub mod bathymetry;
pub mod grid;
pub mod solver;
pub mod state;

pub use bathymetry::Bathymetry;
pub use grid::{GHOST, Grid};
pub use solver::{G, H_DRY, Order, Solver};
pub use state::State;
