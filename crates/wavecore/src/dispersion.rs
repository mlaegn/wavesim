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
use crate::par::{rows1, rows2, rows3, sum_rows};
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
    let h = &s.h[..];
    let b = &bed.b[..];
    let eta = |m: usize| h[m] + b[m];
    rows1(w, phi, |j, row| {
        if j == 0 || j == rows - 1 {
            return;
        }
        for (i, out) in row.iter_mut().enumerate().take(w - 1).skip(1) {
            let k = j * w + i;
            if h[k] < H_DISPERSION_MIN {
                *out = 0.0;
                continue;
            }
            let ex = (eta(k + 1) - eta(k - 1)) * ix;
            let ey = (eta(k + w) - eta(k - w)) * iy;
            let ratio = (eta(k) / h[k]).max(0.0);
            let off = smoothstep(breaking.eta_ratio.0, breaking.eta_ratio.1, ratio)
                .max(smoothstep(breaking.slope.0, breaking.slope.1, ex.hypot(ey)));
            *out = 1.0 - off;
        }
    });
    for j in 0..rows {
        phi[j * w] = phi[j * w + 1];
        phi[j * w + w - 1] = phi[j * w + w - 2];
    }
    for i in 0..w {
        phi[i] = phi[w + i];
        phi[(rows - 1) * w + i] = phi[(rows - 2) * w + i];
    }
}

/// How close each interior cell is to breaking: 0 for smooth water up to 1 where the surface
/// is steep or tall for its depth (the dispersive terms are off there and the shock
/// dissipates the wave). Row-major over the interior, `y` outer. It is 0 in water too thin
/// to carry a wave, so dry land and the film at the waterline never read as breaking.
pub(crate) fn breaking_indicator(
    grid: &Grid,
    bed: &Bathymetry,
    s: &State,
    breaking: Breaking,
    periodic_y: bool,
) -> Vec<f64> {
    let mut filled = s.clone();
    filled.fill_boundaries(grid, periodic_y);
    let mut phi = vec![0.0; grid.cells()];
    switch(grid, bed, &filled, breaking, &mut phi);
    grid.interior()
        .map(|(i, j)| {
            let k = grid.idx(i, j);
            if s.h[k] < H_DISPERSION_MIN {
                0.0
            } else {
                1.0 - phi[k]
            }
        })
        .collect()
}

/// Scratch arrays, all padded like the grid. Allocated once.
pub(crate) struct Work {
    u: Vec<f64>,
    v: Vec<f64>,
    /// 1 in cells deep enough to take part in the dispersive system, 0 elsewhere.
    mask: Vec<f64>,
    /// The velocities with thin cells zeroed, for the divergence and the nonlinear term.
    um: Vec<f64>,
    vm: Vec<f64>,
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
            mask: z(),
            um: z(),
            vm: z(),
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

/// `A p` on the interior, with `A p = h p - grad(c div p)`, restricted to the cells in
/// `mask`: `p` must be zero outside it and the result is zeroed there, which keeps the
/// operator symmetric and positive definite. Fills the ghost cells of `p` first, and uses
/// `f` to hold `c div p`.
#[allow(clippy::too_many_arguments)]
fn matvec(
    grid: &Grid,
    periodic_y: bool,
    heff: &[f64],
    mask: &[f64],
    c: &[f64],
    f: &mut [f64],
    px: &mut [f64],
    py: &mut [f64],
    ax: &mut [f64],
    ay: &mut [f64],
) {
    fill_vector_ghosts(grid, px, py, periodic_y);
    let (px, py): (&[f64], &[f64]) = (px, py);
    let w = grid.width();
    let (ix, iy) = (0.5 / grid.dx, 0.5 / grid.dy);
    let (nx, ny) = (grid.nx, grid.ny);
    // c div p over the interior and one ring of ghost cells.
    rows1(w, f, |j, row| {
        if j < GHOST - 1 || j > ny + GHOST {
            return;
        }
        for (i, out) in row
            .iter_mut()
            .enumerate()
            .take(nx + GHOST + 1)
            .skip(GHOST - 1)
        {
            let k = j * w + i;
            *out = c[k] * ((px[k + 1] - px[k - 1]) * ix + (py[k + w] - py[k - w]) * iy);
        }
    });
    let f: &[f64] = f;
    rows2(w, ax, ay, |j, rx, ry| {
        if j < GHOST || j >= ny + GHOST {
            return;
        }
        for i in GHOST..nx + GHOST {
            let k = j * w + i;
            rx[i] = mask[k] * (heff[k] * px[k] - (f[k + 1] - f[k - 1]) * ix);
            ry[i] = mask[k] * (heff[k] * py[k] - (f[k + w] - f[k - w]) * iy);
        }
    });
}

fn dot(grid: &Grid, ax: &[f64], ay: &[f64], bx: &[f64], by: &[f64]) -> f64 {
    let (w, nx) = (grid.width(), grid.nx);
    sum_rows(GHOST, GHOST + grid.ny, w, |j| {
        let row = j * w + GHOST;
        let mut s = 0.0;
        for k in row..row + nx {
            s += ax[k] * bx[k] + ay[k] * by[k];
        }
        s
    })
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
        let (nx, ny) = (grid.nx, grid.ny);
        let rows = ny + 2 * GHOST;
        let (dx, dy) = (grid.dx, grid.dy);
        let (ix, iy) = (0.5 / dx, 0.5 / dy);
        let nonlinear = self.nonlinear;
        let Work {
            u,
            v,
            mask,
            um,
            vm,
            div,
            c,
            f,
            heff,
            diag_x,
            diag_y,
            wx,
            wy,
            bx,
            by,
            rx,
            ry,
            zx,
            zy,
            px,
            py,
            ax,
            ay,
            stats,
        } = work;
        let (h, hu, hv) = (&s.h[..], &s.hu[..], &s.hv[..]);

        // Water thinner than H_DISPERSION_MIN follows plain shallow water: it takes no part
        // in the dispersive system, neither as unknown nor as a source of divergence. Left
        // in, a film a millimetre deep has an enormous acceleration w = r / h that leaks into
        // its neighbours, and their dispersive force feeds back into it divided by that
        // depth: a loop that grows without bound at the edge of a run-up.
        rows2(w, u, v, |j, ru, rv| {
            for i in 0..w {
                let k = j * w + i;
                ru[i] = velocity(h[k], hu[k]);
                rv[i] = velocity(h[k], hv[k]);
            }
        });
        let (u, v): (&[f64], &[f64]) = (u, v);
        rows3(w, mask, um, vm, |j, rm, rum, rvm| {
            for i in 0..w {
                let k = j * w + i;
                let m = if h[k] >= H_DISPERSION_MIN { 1.0 } else { 0.0 };
                rm[i] = m;
                rum[i] = m * u[k];
                rvm[i] = m * v[k];
            }
        });
        let (mask, um, vm): (&[f64], &[f64], &[f64]) = (mask, um, vm);
        rows1(w, div, |j, row| {
            if j == 0 || j == rows - 1 {
                return;
            }
            for (i, out) in row.iter_mut().enumerate().take(w - 1).skip(1) {
                let k = j * w + i;
                *out = (um[k + 1] - um[k - 1]) * ix + (vm[k + w] - vm[k - w]) * iy;
            }
        });
        let div: &[f64] = div;

        // The coefficient c = phi h^3 / 3 and the nonlinear forcing c N, over the interior
        // and one ring of ghost cells.
        rows2(w, c, f, |j, rc, rf| {
            if j < GHOST - 1 || j > ny + GHOST {
                return;
            }
            for i in GHOST - 1..nx + GHOST + 1 {
                let k = j * w + i;
                let depth = h[k];
                let coeff = if depth < H_DISPERSION_MIN {
                    0.0
                } else {
                    phi[k] * depth * depth * depth / 3.0
                };
                rc[i] = coeff;
                rf[i] = if nonlinear && coeff > 0.0 {
                    let n = um[k] * (div[k + 1] - div[k - 1]) * ix
                        + vm[k] * (div[k + w] - div[k - w]) * iy
                        - div[k] * div[k];
                    coeff * n
                } else {
                    0.0
                };
            }
        });
        let c: &[f64] = c;

        // The effective depth and the Jacobi preconditioner.
        rows3(w, heff, diag_x, diag_y, |j, rh, rdx, rdy| {
            if j < GHOST || j >= ny + GHOST {
                return;
            }
            for i in GHOST..nx + GHOST {
                let k = j * w + i;
                rh[i] = h[k].max(H_DRY);
                if mask[k] == 0.0 {
                    rdx[i] = 1.0;
                    rdy[i] = 1.0;
                } else {
                    rdx[i] = rh[i] + (c[k + 1] + c[k - 1]) * 0.25 / (dx * dx);
                    rdy[i] = rh[i] + (c[k + w] + c[k - w]) * 0.25 / (dy * dy);
                }
            }
        });
        let (heff, diag_x, diag_y): (&[f64], &[f64], &[f64]) = (heff, diag_x, diag_y);

        // The right-hand side, and a zero first guess where a cell takes no part. `f` holds
        // c N here and is scratch for the matrix product afterwards.
        {
            let (lhu_r, lhv_r, f_r): (&[f64], &[f64], &[f64]) = (lhu, lhv, f);
            rows2(w, bx, by, |j, rbx, rby| {
                if j < GHOST || j >= ny + GHOST {
                    return;
                }
                for i in GHOST..nx + GHOST {
                    let k = j * w + i;
                    if mask[k] == 0.0 {
                        rbx[i] = 0.0;
                        rby[i] = 0.0;
                    } else {
                        rbx[i] = lhu_r[k] - u[k] * lh[k] + (f_r[k + 1] - f_r[k - 1]) * ix;
                        rby[i] = lhv_r[k] - v[k] * lh[k] + (f_r[k + w] - f_r[k - w]) * iy;
                    }
                }
            });
        }
        rows2(w, wx, wy, |j, rwx, rwy| {
            if j < GHOST || j >= ny + GHOST {
                return;
            }
            for i in GHOST..nx + GHOST {
                if mask[j * w + i] == 0.0 {
                    rwx[i] = 0.0;
                    rwy[i] = 0.0;
                }
            }
        });
        let (bx, by): (&[f64], &[f64]) = (bx, by);

        // Preconditioned conjugate gradients, starting from the previous solution.
        matvec(grid, periodic_y, heff, mask, c, f, wx, wy, ax, ay);
        {
            let (ax_r, ay_r): (&[f64], &[f64]) = (ax, ay);
            rows2(w, rx, ry, |j, rrx, rry| {
                if j < GHOST || j >= ny + GHOST {
                    return;
                }
                for i in GHOST..nx + GHOST {
                    let k = j * w + i;
                    rrx[i] = bx[k] - ax_r[k];
                    rry[i] = by[k] - ay_r[k];
                }
            });
        }
        precondition(w, nx, ny, rx, ry, diag_x, diag_y, zx, zy);
        copy_interior(w, nx, ny, zx, zy, px, py);
        let bnorm = dot(grid, bx, by, bx, by).sqrt();
        let target = self.tolerance * bnorm + 1e-10;
        let mut rz = dot(grid, rx, ry, zx, zy);
        let mut iterations = 0;
        let mut converged = dot(grid, rx, ry, rx, ry).sqrt() <= target;
        while !converged && iterations < self.max_iterations {
            matvec(grid, periodic_y, heff, mask, c, f, px, py, ax, ay);
            let pap = dot(grid, px, py, ax, ay);
            let alpha = rz / pap;
            axpy(w, nx, ny, alpha, px, py, wx, wy);
            axpy(w, nx, ny, -alpha, ax, ay, rx, ry);
            precondition(w, nx, ny, rx, ry, diag_x, diag_y, zx, zy);
            iterations += 1;
            converged = dot(grid, rx, ry, rx, ry).sqrt() <= target;
            let rz_next = dot(grid, rx, ry, zx, zy);
            let beta = rz_next / rz;
            rz = rz_next;
            {
                let (zx_r, zy_r): (&[f64], &[f64]) = (zx, zy);
                rows2(w, px, py, |j, rpx, rpy| {
                    if j < GHOST || j >= ny + GHOST {
                        return;
                    }
                    for i in GHOST..nx + GHOST {
                        let k = j * w + i;
                        rpx[i] = zx_r[k] + beta * rpx[i];
                        rpy[i] = zy_r[k] + beta * rpy[i];
                    }
                });
            }
        }
        stats.solves += 1;
        stats.iterations += iterations as u64;
        if !converged {
            stats.unconverged += 1;
        }

        // d(hu)/dt = h w + u dh/dt.
        {
            let (wx_r, wy_r): (&[f64], &[f64]) = (wx, wy);
            rows2(w, lhu, lhv, |j, rlu, rlv| {
                if j < GHOST || j >= ny + GHOST {
                    return;
                }
                for i in GHOST..nx + GHOST {
                    let k = j * w + i;
                    if mask[k] > 0.0 {
                        rlu[i] = h[k] * wx_r[k] + u[k] * lh[k];
                        rlv[i] = h[k] * wy_r[k] + v[k] * lh[k];
                    }
                }
            });
        }
    }
}

/// `z = r / diag` over the interior, for both components.
#[allow(clippy::too_many_arguments)]
fn precondition(
    w: usize,
    nx: usize,
    ny: usize,
    rx: &[f64],
    ry: &[f64],
    diag_x: &[f64],
    diag_y: &[f64],
    zx: &mut [f64],
    zy: &mut [f64],
) {
    rows2(w, zx, zy, |j, rzx, rzy| {
        if j < GHOST || j >= ny + GHOST {
            return;
        }
        for i in GHOST..nx + GHOST {
            let k = j * w + i;
            rzx[i] = rx[k] / diag_x[k];
            rzy[i] = ry[k] / diag_y[k];
        }
    });
}

/// `dst = src` over the interior, for both components.
fn copy_interior(
    w: usize,
    nx: usize,
    ny: usize,
    sx: &[f64],
    sy: &[f64],
    dx_: &mut [f64],
    dy_: &mut [f64],
) {
    rows2(w, dx_, dy_, |j, rdx, rdy| {
        if j < GHOST || j >= ny + GHOST {
            return;
        }
        for i in GHOST..nx + GHOST {
            let k = j * w + i;
            rdx[i] = sx[k];
            rdy[i] = sy[k];
        }
    });
}

/// `y += a * x` over the interior, for both components.
#[allow(clippy::too_many_arguments)]
fn axpy(
    w: usize,
    nx: usize,
    ny: usize,
    a: f64,
    xx: &[f64],
    xy: &[f64],
    yx: &mut [f64],
    yy: &mut [f64],
) {
    rows2(w, yx, yy, |j, ryx, ryy| {
        if j < GHOST || j >= ny + GHOST {
            return;
        }
        for i in GHOST..nx + GHOST {
            let k = j * w + i;
            ryx[i] += a * xx[k];
            ryy[i] += a * xy[k];
        }
    });
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
            work.mask[k] = 1.0;
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
            &work.mask,
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
    fn with_thin_cells_masked_out_the_operator_is_still_symmetric_and_positive_definite() {
        // The mask removes shallow cells from the system. A patchy mask must leave `M A M`
        // symmetric and positive definite, and a masked cell must stay exactly zero.
        let (grid, mut work, mut rng) = setup(false);
        crate::grid::mirror_ghosts(&grid, &mut work.c);
        for k in 0..grid.cells() {
            work.mask[k] = if rng.next() > -0.3 { 1.0 } else { 0.0 };
        }
        crate::grid::mirror_ghosts(&grid, &mut work.mask);
        let masked = |grid: &Grid, mask: &[f64], mut v: (Vec<f64>, Vec<f64>)| {
            for (i, j) in grid.interior() {
                let k = grid.idx(i, j);
                v.0[k] *= mask[k];
                v.1[k] *= mask[k];
            }
            v
        };
        assert!(
            grid.interior()
                .any(|(i, j)| work.mask[grid.idx(i, j)] == 0.0),
            "the test mask removes nothing"
        );
        for _ in 0..5 {
            let (x1, y1) = masked(&grid, &work.mask, random_field(&grid, &mut rng));
            let (x2, y2) = masked(&grid, &work.mask, random_field(&grid, &mut rng));
            let (ax1, ay1) = apply_a(&grid, &mut work, false, &x1, &y1);
            let (ax2, ay2) = apply_a(&grid, &mut work, false, &x2, &y2);
            let left = dot(&grid, &x2, &y2, &ax1, &ay1);
            let right = dot(&grid, &x1, &y1, &ax2, &ay2);
            assert!(
                (left - right).abs() < 1e-10 * left.abs().max(1.0),
                "<y, A x> = {left}, <x, A y> = {right}"
            );
            assert!(dot(&grid, &x1, &y1, &ax1, &ay1) > 0.0);
            for (i, j) in grid.interior() {
                let k = grid.idx(i, j);
                if work.mask[k] == 0.0 {
                    assert_eq!(
                        (ax1[k], ay1[k]),
                        (0.0, 0.0),
                        "a masked cell was not left at zero"
                    );
                }
            }
        }
    }

    #[test]
    fn a_thin_film_neither_feels_nor_feeds_the_dispersive_system() {
        // Regression: a film a millimetre deep at the edge of a run-up once had an enormous
        // acceleration that leaked into its neighbours and fed back divided by its depth,
        // growing without bound. Its velocity must not change what its neighbours get, and
        // its own momentum rate must pass through untouched.
        let grid = Grid::new(12, 8, 2.0, 2.0);
        let film = grid.idx(GHOST + 6, GHOST + 4);
        let run = |film_speed: f64| {
            let mut s = State::zeros(&grid);
            for (i, j) in grid.interior() {
                let k = grid.idx(i, j);
                s.h[k] = 2.0;
                s.hu[k] = 2.0 * 0.3 * (0.7 * i as f64).sin();
                s.hv[k] = 2.0 * 0.2 * (0.5 * j as f64).cos();
            }
            s.h[film] = 0.001;
            s.hu[film] = 0.001 * film_speed;
            s.hv[film] = 0.0;
            s.fill_boundaries(&grid, false);
            let bed = Bathymetry::flat(&grid, 2.0);
            let mut phi = vec![1.0; grid.cells()];
            phi[film] = 0.0;
            let lh: Vec<f64> = (0..grid.cells()).map(|k| 0.01 * (k % 7) as f64).collect();
            let mut lhu: Vec<f64> = (0..grid.cells()).map(|k| 0.05 * (k % 5) as f64).collect();
            let mut lhv: Vec<f64> = (0..grid.cells()).map(|k| 0.03 * (k % 3) as f64).collect();
            let before = (lhu[film], lhv[film]);
            let mut work = Work::new(&grid);
            let _ = bed;
            Dispersion::default().apply(&grid, false, &s, &lh, &mut lhu, &mut lhv, &phi, &mut work);
            (lhu, lhv, before)
        };
        let (a_u, a_v, before) = run(1.0);
        let (b_u, b_v, _) = run(1000.0);
        for (i, j) in grid.interior() {
            let k = grid.idx(i, j);
            if k == film {
                assert_eq!((a_u[k], a_v[k]), before, "the film's own rate was altered");
                assert_eq!((b_u[k], b_v[k]), before, "the film's own rate was altered");
            } else {
                assert_eq!(a_u[k], b_u[k], "cell ({i}, {j}) felt the film's speed");
                assert_eq!(a_v[k], b_v[k], "cell ({i}, {j}) felt the film's speed");
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
