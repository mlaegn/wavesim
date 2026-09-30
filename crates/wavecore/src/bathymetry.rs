use crate::grid::{Grid, mirror_ghosts};

/// Bed elevation `b` in metres relative to still-water level (negative is underwater).
/// Stored padded; ghost cells mirror the interior.
#[derive(Clone, Debug)]
pub struct Bathymetry {
    pub b: Vec<f64>,
}

impl Bathymetry {
    /// Build from a function of cell-centre coordinates `(x, y)` in metres.
    pub fn from_fn(grid: &Grid, f: impl Fn(f64, f64) -> f64) -> Self {
        let mut b = vec![0.0; grid.cells()];
        for (i, j) in grid.interior() {
            let (x, y) = grid.centre(i, j);
            b[grid.idx(i, j)] = f(x, y);
        }
        mirror_ghosts(grid, &mut b);
        Self { b }
    }

    /// Build from interior values, row-major with `j` (y) outer and `i` (x) inner, as
    /// stored in a bed file. Panics if the length does not match the grid.
    pub fn from_interior(grid: &Grid, values: &[f64]) -> Self {
        assert_eq!(
            values.len(),
            grid.nx * grid.ny,
            "bed has {} values, grid needs {}",
            values.len(),
            grid.nx * grid.ny
        );
        let mut b = vec![0.0; grid.cells()];
        for (n, (i, j)) in grid.interior().enumerate() {
            b[grid.idx(i, j)] = values[n];
        }
        mirror_ghosts(grid, &mut b);
        Self { b }
    }

    /// Constant depth: a flat bed `depth` metres below still water.
    pub fn flat(grid: &Grid, depth: f64) -> Self {
        Self::from_fn(grid, |_, _| -depth)
    }

    /// Plane beach rising towards +x with the given slope (rise over run).
    /// The bed starts `offshore_depth` below still water at x = 0.
    pub fn plane_beach(grid: &Grid, slope: f64, offshore_depth: f64) -> Self {
        Self::from_fn(grid, |x, _| -offshore_depth + slope * x)
    }

    /// Deterministic rough bed, used to test that still water stays still.
    /// Oscillates between `-mean_depth - amplitude` and `-mean_depth + amplitude`.
    pub fn rough(grid: &Grid, mean_depth: f64, amplitude: f64) -> Self {
        let lx = grid.nx as f64 * grid.dx;
        let ly = grid.ny as f64 * grid.dy;
        let tau = std::f64::consts::TAU;
        Self::from_fn(grid, |x, y| {
            let s = 0.6 * (tau * 3.0 * x / lx).sin() * (tau * 2.0 * y / ly).cos()
                + 0.4 * (tau * 7.0 * x / lx + 1.3).sin();
            -mean_depth + amplitude * s
        })
    }
}
