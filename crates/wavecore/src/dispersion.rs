//! Weak dispersion: the Serre-Green-Naghdi (SGN) correction to the shallow-water
//! equations, switched off where waves break.
//!
//! Shallow water has no dispersion, so a wave that is tall for its depth steepens into
//! a shock wherever it is. Real water spreads the steep front out, because short
//! components travel slower than long ones, and only breaks when the wave is nearly as
//! tall as the water is deep. SGN puts that back.
//!
//! For a flat bed the depth-averaged momentum equation becomes
//!
//! ```text
//! d(hu)/dt + div(h u u) + g h grad(eta) = (1/3) grad( h^3 Gamma )
//! Gamma = d(div u)/dt + u . grad(div u) - (div u)^2
//! ```
//!
//! (from a vertical velocity that grows linearly from the bed). The unknown `w = du/dt`
//! appears inside `Gamma`, so each evaluation of the time derivative solves
//!
//! ```text
//! h w - grad( c div w ) = r,    c = phi h^3 / 3,
//! r = (hyperbolic momentum rate) - u (mass rate) + grad( c N ),   N = u.grad(div u) - (div u)^2
//! ```
//!
//! This operator is symmetric and positive definite, so conjugate gradients with a
//! Jacobi preconditioner solves it; the previous solution is a good first guess. `phi`
//! is the hybrid breaking switch: 1 where the wave is smooth, 0 where it is steep or
//! tall for its depth, which leaves plain shallow water there and lets the shock
//! dissipate the wave. Keeping `phi` inside the gradient keeps the operator symmetric.
//!
//! Bed-slope terms of the dispersive operator are not included; they are of order
//! slope x (kh)^2.

use crate::bathymetry::Bathymetry;
use crate::grid::{GHOST, Grid, fill_vector_ghosts};
use crate::solver::{H_DRY, velocity};
use crate::state::State;

/// Depths below this get no dispersion; a few centimetres of water has no wave to disperse.
const H_DISPERSION_MIN: f64 = 0.05;

/// When to switch the dispersive terms off. Each pair is `(starts fading, fully off)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Breaking {
    /// Surface height above still water relative to the local depth, `eta / h`.
    pub eta_ratio: (f64, f64),
    /// Magnitude of the surface slope, `|grad eta|`.
    pub slope: (f64, f64),
}

impl Default for Breaking {
    fn default() -> Self {
        Self {
            eta_ratio: (0.30, 0.55),
            slope: (0.20, 0.45),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Dispersion {
    pub breaking: Breaking,
    /// Keep the nonlinear dispersive terms `u.grad(div u) - (div u)^2`. Switching them off
    /// leaves the linear (Peregrine-type) dispersion.
    pub nonlinear: bool,
    /// Stop the linear solve when the residual is this fraction of the right-hand side.
    /// 1e-4 and 1e-6 give the same waves to three digits; 1e-4 takes half the iterations.
    pub tolerance: f64,
    pub max_iterations: usize,
}

impl Default for Dispersion {
    fn default() -> Self {
        Self {
            breaking: Breaking::default(),
            nonlinear: true,
            tolerance: 1e-4,
            max_iterations: 200,
        }
    }
}

/// How hard the linear solves have been working.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct DispersionStats {
    pub solves: u64,
    pub iterations: u64,
    /// Solves that stopped at `max_iterations` without reaching the tolerance.
    pub unconverged: u64,
}

fn smoothstep(lo: f64, hi: f64, x: f64) -> f64 {
    let t = ((x - lo) / (hi - lo)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// The smoothness switch `phi` for every cell: 1 where the surface is smooth, falling to 0
/// where it is steep or tall for its depth (a front about to break, a bore, the swash).
/// Both the dispersive terms and the third-order reconstruction are scaled by it.
///
/// It is computed everywhere the stencils reach; the outermost ring copies its neighbour.
pub(crate) fn switch(
    grid: &Grid,
    bed: &Bathymetry,
    s: &State,
    breaking: Breaking,
    phi: &mut [f64],
) {
    let w = grid.width();
    let rows = grid.ny + 2 * GHOST;
    let (ix, iy) = (0.5 / grid.dx, 0.5 / grid.dy);
    let eta = |m: usize| s.h[m] + bed.b[m];
    for j in 1..rows - 1 {
        for i in 1..w - 1 {
            let k = j * w + i;
            let h = s.h[k];
            if h < H_DISPERSION_MIN {
                phi[k] = 0.0;
                continue;
            }
            let ex = (eta(k + 1) - eta(k - 1)) * ix;
            let ey = (eta(k + w) - eta(k - w)) * iy;
            let ratio = (eta(k) / h).max(0.0);
            let off = smoothstep(breaking.eta_ratio.0, breaking.eta_ratio.1, ratio)
                .max(smoothstep(breaking.slope.0, breaking.slope.1, ex.hypot(ey)));
            phi[k] = 1.0 - off;
        }
    }
    for j in 0..rows {
        phi[j * w] = phi[j * w + 1];
        phi[j * w + w - 1] = phi[j * w + w - 2];
    }
    for i in 0..w {
        phi[i] = phi[w + i];
        phi[(rows - 1) * w + i] = phi[(rows - 2) * w + i];
    }
}

/// Scratch arrays, all padded like the grid. Allocated once.
pub(crate) struct Work {
    u: Vec<f64>,
    v: Vec<f64>,
    div: Vec<f64>,
    /// `c = phi h^3 / 3`.
    c: Vec<f64>,
    /// `c * (something)`: first `c N`, then reused as `c div p` in the matrix product.
    f: Vec<f64>,
    heff: Vec<f64>,
    diag_x: Vec<f64>,
    diag_y: Vec<f64>,
    /// Solution `w = du/dt`, kept between calls as the next first guess.
    wx: Vec<f64>,
    wy: Vec<f64>,
    bx: Vec<f64>,
    by: Vec<f64>,
    rx: Vec<f64>,
    ry: Vec<f64>,
    zx: Vec<f64>,
    zy: Vec<f64>,
    px: Vec<f64>,
    py: Vec<f64>,
    ax: Vec<f64>,
    ay: Vec<f64>,
    pub(crate) stats: DispersionStats,
}

impl Work {
    pub(crate) fn new(grid: &Grid) -> Self {
        let n = grid.cells();
        let z = || vec![0.0; n];
        Self {
            u: z(),
            v: z(),
            div: z(),
            c: z(),
            f: z(),
            heff: z(),
            diag_x: z(),
            diag_y: z(),
            wx: z(),
            wy: z(),
            bx: z(),
            by: z(),
            rx: z(),
            ry: z(),
            zx: z(),
            zy: z(),
            px: z(),
            py: z(),
            ax: z(),
            ay: z(),
            stats: DispersionStats::default(),
        }
    }
}

/// `A p` on the interior, with `A p = h p - grad(c div p)`. Fills the ghost cells of `p`
/// first, and uses `f` to hold `c div p`.
#[allow(clippy::too_many_arguments)]
fn matvec(
    grid: &Grid,
    periodic_y: bool,
    heff: &[f64],
    c: &[f64],
    f: &mut [f64],
    px: &mut [f64],
    py: &mut [f64],
    ax: &mut [f64],
    ay: &mut [f64],
) {
    fill_vector_ghosts(grid, px, py, periodic_y);
    let w = grid.width();
    let (ix, iy) = (0.5 / grid.dx, 0.5 / grid.dy);
    // c div p over the interior and one ring of ghost cells.
    for j in GHOST - 1..grid.ny + GHOST + 1 {
        for i in GHOST - 1..grid.nx + GHOST + 1 {
            let k = j * w + i;
            f[k] = c[k] * ((px[k + 1] - px[k - 1]) * ix + (py[k + w] - py[k - w]) * iy);
        }
    }
    for (i, j) in grid.interior() {
        let k = j * w + i;
        ax[k] = heff[k] * px[k] - (f[k + 1] - f[k - 1]) * ix;
        ay[k] = heff[k] * py[k] - (f[k + w] - f[k - w]) * iy;
    }
}

fn dot(grid: &Grid, ax: &[f64], ay: &[f64], bx: &[f64], by: &[f64]) -> f64 {
    let w = grid.width();
    let mut s = 0.0;
    for j in GHOST..grid.ny + GHOST {
        let row = j * w;
        for k in row + GHOST..row + GHOST + grid.nx {
            s += ax[k] * bx[k] + ay[k] * by[k];
        }
    }
    s
}

impl Dispersion {
    /// Replace the shallow-water momentum rates `(lhu, lhv)` with the SGN ones. `s` must
    /// have its ghost cells filled; `lh` is the mass rate.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn apply(
        &self,
        grid: &Grid,
        periodic_y: bool,
        s: &State,
        lh: &[f64],
        lhu: &mut [f64],
        lhv: &mut [f64],
        phi: &[f64],
        work: &mut Work,
    ) {
        let w = grid.width();
        let rows = grid.ny + 2 * GHOST;
        let (dx, dy) = (grid.dx, grid.dy);
        let (ix, iy) = (0.5 / dx, 0.5 / dy);
        let wk = work;

        for k in 0..grid.cells() {
            wk.u[k] = velocity(s.h[k], s.hu[k]);
            wk.v[k] = velocity(s.h[k], s.hv[k]);
        }
        for j in 1..rows - 1 {
            for i in 1..w - 1 {
                let k = j * w + i;
                wk.div[k] = (wk.u[k + 1] - wk.u[k - 1]) * ix + (wk.v[k + w] - wk.v[k - w]) * iy;
            }
        }

        // The coefficient c = phi h^3 / 3 and the nonlinear forcing c N, over the interior
        // and one ring of ghost cells.
        for j in GHOST - 1..grid.ny + GHOST + 1 {
            for i in GHOST - 1..grid.nx + GHOST + 1 {
                let k = j * w + i;
                let h = s.h[k];
                let c = if h < H_DISPERSION_MIN {
                    0.0
                } else {
                    phi[k] * h * h * h / 3.0
                };
                wk.c[k] = c;
                wk.f[k] = if self.nonlinear && c > 0.0 {
                    let n = wk.u[k] * (wk.div[k + 1] - wk.div[k - 1]) * ix
                        + wk.v[k] * (wk.div[k + w] - wk.div[k - w]) * iy
                        - wk.div[k] * wk.div[k];
                    c * n
                } else {
                    0.0
                };
            }
        }

        // Right-hand side, the effective depth, and the Jacobi preconditioner.
        for (i, j) in grid.interior() {
            let k = j * w + i;
            let h = s.h[k];
            wk.heff[k] = h.max(H_DRY);
            if h <= H_DRY {
                wk.bx[k] = 0.0;
                wk.by[k] = 0.0;
                wk.wx[k] = 0.0;
                wk.wy[k] = 0.0;
            } else {
                wk.bx[k] = lhu[k] - wk.u[k] * lh[k] + (wk.f[k + 1] - wk.f[k - 1]) * ix;
                wk.by[k] = lhv[k] - wk.v[k] * lh[k] + (wk.f[k + w] - wk.f[k - w]) * iy;
            }
            wk.diag_x[k] = wk.heff[k] + (wk.c[k + 1] + wk.c[k - 1]) * 0.25 / (dx * dx);
            wk.diag_y[k] = wk.heff[k] + (wk.c[k + w] + wk.c[k - w]) * 0.25 / (dy * dy);
        }

        // Preconditioned conjugate gradients, starting from the previous solution.
        matvec(
            grid, periodic_y, &wk.heff, &wk.c, &mut wk.f, &mut wk.wx, &mut wk.wy, &mut wk.ax,
            &mut wk.ay,
        );
        for (i, j) in grid.interior() {
            let k = j * w + i;
            wk.rx[k] = wk.bx[k] - wk.ax[k];
            wk.ry[k] = wk.by[k] - wk.ay[k];
            wk.zx[k] = wk.rx[k] / wk.diag_x[k];
            wk.zy[k] = wk.ry[k] / wk.diag_y[k];
            wk.px[k] = wk.zx[k];
            wk.py[k] = wk.zy[k];
        }
        let bnorm = dot(grid, &wk.bx, &wk.by, &wk.bx, &wk.by).sqrt();
        let target = self.tolerance * bnorm + 1e-10;
        let mut rz = dot(grid, &wk.rx, &wk.ry, &wk.zx, &wk.zy);
        let mut iterations = 0;
        let mut converged = dot(grid, &wk.rx, &wk.ry, &wk.rx, &wk.ry).sqrt() <= target;
        while !converged && iterations < self.max_iterations {
            matvec(
                grid, periodic_y, &wk.heff, &wk.c, &mut wk.f, &mut wk.px, &mut wk.py, &mut wk.ax,
                &mut wk.ay,
            );
            let pap = dot(grid, &wk.px, &wk.py, &wk.ax, &wk.ay);
            let alpha = rz / pap;
            for (i, j) in grid.interior() {
                let k = j * w + i;
                wk.wx[k] += alpha * wk.px[k];
                wk.wy[k] += alpha * wk.py[k];
                wk.rx[k] -= alpha * wk.ax[k];
                wk.ry[k] -= alpha * wk.ay[k];
                wk.zx[k] = wk.rx[k] / wk.diag_x[k];
                wk.zy[k] = wk.ry[k] / wk.diag_y[k];
            }
            iterations += 1;
            converged = dot(grid, &wk.rx, &wk.ry, &wk.rx, &wk.ry).sqrt() <= target;
            let rz_next = dot(grid, &wk.rx, &wk.ry, &wk.zx, &wk.zy);
            let beta = rz_next / rz;
            rz = rz_next;
            for (i, j) in grid.interior() {
                let k = j * w + i;
                wk.px[k] = wk.zx[k] + beta * wk.px[k];
                wk.py[k] = wk.zy[k] + beta * wk.py[k];
            }
        }
        wk.stats.solves += 1;
        wk.stats.iterations += iterations as u64;
        if !converged {
            wk.stats.unconverged += 1;
        }

        // d(hu)/dt = h w + u dh/dt.
        for (i, j) in grid.interior() {
            let k = j * w + i;
            lhu[k] = s.h[k] * wk.wx[k] + wk.u[k] * lh[k];
            lhv[k] = s.h[k] * wk.wy[k] + wk.v[k] * lh[k];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random numbers in [-1, 1).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        }
    }

    fn setup(periodic: bool) -> (Grid, Work, Lcg) {
        let grid = Grid::new(13, 9, 2.0, 3.0);
        let mut work = Work::new(&grid);
        let mut rng = Lcg(7);
        for k in 0..grid.cells() {
            work.heff[k] = 0.5 + 4.0 * (rng.next() + 1.0);
            work.c[k] = 3.0 * (rng.next() + 1.0);
        }
        // c and heff are only read at the interior and one ring around it, but filling
        // everything keeps the test simple.
        let _ = periodic;
        (grid, work, rng)
    }

    fn random_field(grid: &Grid, rng: &mut Lcg) -> (Vec<f64>, Vec<f64>) {
        let mut x = vec![0.0; grid.cells()];
        let mut y = vec![0.0; grid.cells()];
        for (i, j) in grid.interior() {
            x[grid.idx(i, j)] = rng.next();
            y[grid.idx(i, j)] = rng.next();
        }
        (x, y)
    }

    fn apply_a(
        grid: &Grid,
        work: &mut Work,
        periodic: bool,
        x: &[f64],
        y: &[f64],
    ) -> (Vec<f64>, Vec<f64>) {
        let (mut px, mut py) = (x.to_vec(), y.to_vec());
        matvec(
            grid,
            periodic,
            &work.heff,
            &work.c,
            &mut work.f,
            &mut px,
            &mut py,
            &mut work.ax,
            &mut work.ay,
        );
        (work.ax.clone(), work.ay.clone())
    }

    #[test]
    fn the_operator_is_symmetric_and_positive_definite() {
        for periodic in [false, true] {
            let (grid, mut work, mut rng) = setup(periodic);
            // Make c and heff consistent with the wrap when periodic: copy the interior
            // rows into the ghost rows so the operator sees one periodic medium.
            if periodic {
                crate::grid::wrap_ghosts_y(&grid, &mut work.c);
                crate::grid::wrap_ghosts_y(&grid, &mut work.heff);
            }
            // Mirror c in the x ghosts, as the solver's mirrored state would give.
            crate::grid::mirror_ghosts(&grid, &mut work.c);
            if periodic {
                crate::grid::wrap_ghosts_y(&grid, &mut work.c);
            }
            for _ in 0..5 {
                let (x1, y1) = random_field(&grid, &mut rng);
                let (x2, y2) = random_field(&grid, &mut rng);
                let (ax1, ay1) = apply_a(&grid, &mut work, periodic, &x1, &y1);
                let (ax2, ay2) = apply_a(&grid, &mut work, periodic, &x2, &y2);
                let left = dot(&grid, &x2, &y2, &ax1, &ay1);
                let right = dot(&grid, &x1, &y1, &ax2, &ay2);
                assert!(
                    (left - right).abs() < 1e-10 * left.abs().max(1.0),
                    "periodic={periodic}: <y, A x> = {left}, <x, A y> = {right}"
                );
                assert!(
                    dot(&grid, &x1, &y1, &ax1, &ay1) > 0.0,
                    "periodic={periodic}: not positive"
                );
            }
        }
    }

    #[test]
    fn smoothstep_is_a_clamped_cubic() {
        assert_eq!(smoothstep(1.0, 3.0, 0.0), 0.0);
        assert_eq!(smoothstep(1.0, 3.0, 4.0), 1.0);
        assert!((smoothstep(1.0, 3.0, 2.0) - 0.5).abs() < 1e-15);
    }
}
