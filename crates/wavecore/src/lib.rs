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
//!
//! Generating waves: a wave-maker radiates a monochromatic wave into a channel, a sponge
//! layer absorbs what reaches the far end, and Manning friction is optional.
//!
//! ```
//! use wavecore::{Bathymetry, Grid, Solver, Sponge, State, WaveMaker};
//!
//! let grid = Grid::new(400, 4, 0.5, 0.5);
//! let bed = Bathymetry::flat(&grid, 2.0);
//! let maker = WaveMaker {
//!     x: 100.0,           // source line, metres
//!     amplitude: 0.05,    // metres
//!     period: 12.0,       // seconds
//!     angle: 0.0,         // radians from +x
//!     sigma: 2.0,         // source width, metres
//!     ramp_periods: 2.0,
//! };
//! let solver = Solver::new(grid, bed)
//!     .with_wavemaker(maker)
//!     .with_sponge(Sponge::new(60, 1.0).edges(true, true, false, false))
//!     .with_manning(0.02);
//! let mut state = State::lake_at_rest(&grid, &solver.bed, 0.0);
//!
//! let dt = solver.stable_dt(&state);
//! solver.step(&mut state, dt);
//! assert!(state.time > 0.0);
//! ```

pub mod bathymetry;
pub mod dispersion;
pub mod forcing;
pub mod grid;
mod par;
pub mod solver;
pub mod state;

pub use bathymetry::Bathymetry;
pub use dispersion::{Breaking, Dispersion, DispersionStats};
pub use forcing::{LinearWave, Relaxation, Sponge, WaveMaker, linear_wave, manning_factor};
pub use grid::{GHOST, Grid};
pub use solver::{G, H_DRY, Order, Solver};
pub use state::State;
