//! First-order, well-balanced finite-volume solver for the nonlinear shallow-water
//! equations. Hydrostatic reconstruction (Audusse et al. 2004) keeps still water
//! still over any bed and keeps depth non-negative at moving shorelines. HLL flux.
//!
//! Every face computation is a pure function of the two cells that share it, so the
//! same kernels can be ported one-to-one to a compute shader.

use crate::bathymetry::Bathymetry;
use crate::grid::Grid;
use crate::state::State;

pub const G: f64 = 9.81;
/// Below this depth a cell counts as dry and carries no velocity.
pub const H_DRY: f64 = 1e-8;

/// What one cell contributes to a face: depth, momentum normal to the face,
/// momentum along it, and bed elevation.
#[derive(Clone, Copy)]
struct Side {
    h: f64,
    hn: f64,
    ht: f64,
    b: f64,
}

/// Numerical flux across a face plus the hydrostatic pressure corrections that
/// the cell on each side adds to its normal-momentum flux.
struct Face {
    flux: [f64; 3],
    corr_left: f64,
    corr_right: f64,
}

fn velocity(h: f64, hm: f64) -> f64 {
    if h > H_DRY { hm / h } else { 0.0 }
}

fn face(l: Side, r: Side) -> Face {
    // Hydrostatic reconstruction: both sides see the same bed step, so a flat
    // free surface produces exactly balanced fluxes.
    let b_star = l.b.max(r.b);
    let hl = (l.h + l.b - b_star).max(0.0);
    let hr = (r.h + r.b - b_star).max(0.0);
    let ul = velocity(l.h, l.hn);
    let ur = velocity(r.h, r.hn);
    let vl = velocity(l.h, l.ht);
    let vr = velocity(r.h, r.ht);

    let cl = (G * hl).sqrt();
    let cr = (G * hr).sqrt();
    let s_l = (ul - cl).min(ur - cr);
    let s_r = (ul + cl).max(ur + cr);

    let f_l = [hl * ul, hl * ul * ul + 0.5 * G * hl * hl, hl * ul * vl];
    let f_r = [hr * ur, hr * ur * ur + 0.5 * G * hr * hr, hr * ur * vr];
    let u_l = [hl, hl * ul, hl * vl];
    let u_r = [hr, hr * ur, hr * vr];

    let flux = if hl <= H_DRY && hr <= H_DRY {
        [0.0; 3]
    } else if s_l >= 0.0 {
        f_l
    } else if s_r <= 0.0 {
        f_r
    } else {
        let inv = 1.0 / (s_r - s_l);
        std::array::from_fn(|k| (s_r * f_l[k] - s_l * f_r[k] + s_l * s_r * (u_r[k] - u_l[k])) * inv)
    };

    Face {
        flux,
        corr_left: 0.5 * G * (l.h * l.h - hl * hl),
        corr_right: 0.5 * G * (r.h * r.h - hr * hr),
    }
}

pub struct Solver {
    pub grid: Grid,
    pub bed: Bathymetry,
    pub cfl: f64,
}

impl Solver {
    pub fn new(grid: Grid, bed: Bathymetry) -> Self {
        Self {
            grid,
            bed,
            cfl: 0.4,
        }
    }

    /// Largest stable time step for the current state. Infinite if nothing can move.
    pub fn stable_dt(&self, s: &State) -> f64 {
        let g = &self.grid;
        let mut rate: f64 = 0.0;
        for (i, j) in g.interior() {
            let k = g.idx(i, j);
            let c = (G * s.h[k]).sqrt();
            let u = velocity(s.h[k], s.hu[k]).abs();
            let v = velocity(s.h[k], s.hv[k]).abs();
            rate = rate.max((u + c) / g.dx + (v + c) / g.dy);
        }
        if rate > 0.0 {
            self.cfl / rate
        } else {
            f64::INFINITY
        }
    }

    /// Advance one forward-Euler step of size `dt` with reflective walls.
    pub fn step(&self, s: &mut State, dt: f64) {
        let g = &self.grid;
        s.fill_walls(g);
        let mut next = s.clone();
        let (rx, ry) = (dt / g.dx, dt / g.dy);
        let side = |k: usize, normal_is_x: bool| Side {
            h: s.h[k],
            hn: if normal_is_x { s.hu[k] } else { s.hv[k] },
            ht: if normal_is_x { s.hv[k] } else { s.hu[k] },
            b: self.bed.b[k],
        };

        for (i, j) in g.interior() {
            let k = g.idx(i, j);
            let c = side(k, true);
            let west = face(side(g.idx(i - 1, j), true), c);
            let east = face(c, side(g.idx(i + 1, j), true));
            let cy = side(k, false);
            let south = face(side(g.idx(i, j - 1), false), cy);
            let north = face(cy, side(g.idx(i, j + 1), false));

            next.h[k] =
                s.h[k] - rx * (east.flux[0] - west.flux[0]) - ry * (north.flux[0] - south.flux[0]);
            next.hu[k] = s.hu[k]
                - rx * ((east.flux[1] + east.corr_left) - (west.flux[1] + west.corr_right))
                - ry * (north.flux[2] - south.flux[2]);
            next.hv[k] = s.hv[k]
                - rx * (east.flux[2] - west.flux[2])
                - ry * ((north.flux[1] + north.corr_left) - (south.flux[1] + south.corr_right));

            if next.h[k] < H_DRY {
                // Only rounding-level negatives are expected under the CFL limit.
                next.h[k] = next.h[k].max(0.0);
                next.hu[k] = 0.0;
                next.hv[k] = 0.0;
            }
        }
        *s = next;
    }
}
