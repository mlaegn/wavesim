use crate::grid::Grid;

/// Bed elevation `b` in metres relative to still-water level (negative is underwater).
/// Stored padded; ghost cells copy the nearest interior cell.
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
        fill_ghosts_copy(grid, &mut b);
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

/// Copy the nearest interior value into every ghost cell.
pub(crate) fn fill_ghosts_copy(grid: &Grid, a: &mut [f64]) {
    let (nx, ny) = (grid.nx, grid.ny);
    for j in 1..=ny {
        a[grid.idx(0, j)] = a[grid.idx(1, j)];
        a[grid.idx(nx + 1, j)] = a[grid.idx(nx, j)];
    }
    for i in 0..nx + 2 {
        a[grid.idx(i, 0)] = a[grid.idx(i, 1)];
        a[grid.idx(i, ny + 1)] = a[grid.idx(i, ny)];
    }
}
