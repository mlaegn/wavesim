use crate::bathymetry::Bathymetry;
use crate::grid::{GHOST, Grid, mirror_ghosts};

/// Conserved variables: water depth `h` and depth-integrated momentum `hu`, `hv`.
/// All arrays are padded, see [`Grid`].
#[derive(Clone, Debug)]
pub struct State {
    pub h: Vec<f64>,
    pub hu: Vec<f64>,
    pub hv: Vec<f64>,
}

impl State {
    pub fn zeros(grid: &Grid) -> Self {
        let n = grid.cells();
        Self {
            h: vec![0.0; n],
            hu: vec![0.0; n],
            hv: vec![0.0; n],
        }
    }

    /// Still water at surface elevation `eta`. Cells whose bed is above `eta` are dry.
    pub fn lake_at_rest(grid: &Grid, bed: &Bathymetry, eta: f64) -> Self {
        let mut s = Self::zeros(grid);
        for (i, j) in grid.interior() {
            let k = grid.idx(i, j);
            s.h[k] = (eta - bed.b[k]).max(0.0);
        }
        s.fill_walls(grid);
        s
    }

    /// Reflective walls: ghost cells mirror depth and tangential momentum,
    /// and flip the momentum component normal to the wall.
    pub fn fill_walls(&mut self, grid: &Grid) {
        mirror_ghosts(grid, &mut self.h);
        mirror_ghosts(grid, &mut self.hu);
        mirror_ghosts(grid, &mut self.hv);
        let (nx, ny, w) = (grid.nx, grid.ny, grid.width());
        for j in 0..ny + 2 * GHOST {
            for k in 0..GHOST {
                let (l, r) = (grid.idx(GHOST - 1 - k, j), grid.idx(GHOST + nx + k, j));
                self.hu[l] = -self.hu[l];
                self.hu[r] = -self.hu[r];
            }
        }
        for i in 0..w {
            for k in 0..GHOST {
                let (lo, hi) = (grid.idx(i, GHOST - 1 - k), grid.idx(i, GHOST + ny + k));
                self.hv[lo] = -self.hv[lo];
                self.hv[hi] = -self.hv[hi];
            }
        }
    }

    /// Total water volume in cubic metres.
    pub fn volume(&self, grid: &Grid) -> f64 {
        grid.interior()
            .map(|(i, j)| self.h[grid.idx(i, j)])
            .sum::<f64>()
            * grid.dx
            * grid.dy
    }

    /// Largest absolute depth-integrated momentum component over the interior.
    pub fn max_abs_momentum(&self, grid: &Grid) -> f64 {
        grid.interior()
            .map(|(i, j)| {
                let k = grid.idx(i, j);
                self.hu[k].abs().max(self.hv[k].abs())
            })
            .fold(0.0, f64::max)
    }

    /// Smallest depth over the interior; negative means the scheme lost positivity.
    pub fn min_depth(&self, grid: &Grid) -> f64 {
        grid.interior()
            .map(|(i, j)| self.h[grid.idx(i, j)])
            .fold(f64::INFINITY, f64::min)
    }
}
