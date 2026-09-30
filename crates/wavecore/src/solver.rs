//! Well-balanced finite-volume solver for the nonlinear shallow-water equations.
//!
//! Hydrostatic reconstruction (Audusse et al. 2004) keeps still water still over any
//! bed and keeps depth non-negative at moving shorelines. The flux is HLL. Second order
//! comes from MUSCL reconstruction of the free surface and velocities with an MC
//! limiter, advanced by SSP-RK2 (Heun). Wherever a cell or a neighbour is nearly dry,
//! or reconstruction would produce a face without water, the cell drops to first order.
//!
//! Every face computation is a pure function of the states on its two sides, so the
//! same kernels can be ported one-to-one to a compute shader.

use crate::bathymetry::Bathymetry;
use crate::grid::Grid;
use crate::state::State;

pub const G: f64 = 9.81;
/// Below this depth a cell counts as dry and carries no velocity.
pub const H_DRY: f64 = 1e-8;
/// Second-order reconstruction needs the cell and both neighbours at least this deep.
const H_SECOND: f64 = 1e-3;
/// ...and every reconstructed face depth at least this deep.
const H_FACE_MIN: f64 = 1e-4;

/// Spatial and temporal order of the scheme.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Order {
    /// Piecewise-constant states, forward Euler.
    First,
    /// MUSCL reconstruction with an MC limiter, SSP-RK2.
    Second,
}

/// Cell-averaged primitive variables along one axis: depth, velocity normal to the
/// axis' faces, velocity along them, and bed elevation.
#[derive(Clone, Copy)]
struct Cell {
    h: f64,
    un: f64,
    ut: f64,
    b: f64,
}

/// State on one side of a face: depth, normal and tangential velocity, and the bed
/// elevation there.
#[derive(Clone, Copy)]
struct Side {
    h: f64,
    un: f64,
    ut: f64,
    z: f64,
}

/// Numerical flux across a face plus the hydrostatic pressure corrections that the
/// cell on each side adds to its normal-momentum flux.
struct Face {
    flux: [f64; 3],
    corr_left: f64,
    corr_right: f64,
}

/// Everything one axis contributes to a cell: the fluxes through its two faces on
/// that axis and the centred bed-slope source for the normal momentum.
struct Axis {
    low: Face,
    high: Face,
    source: f64,
}

fn velocity(h: f64, hm: f64) -> f64 {
    if h > H_DRY { hm / h } else { 0.0 }
}

fn minmod3(a: f64, b: f64, c: f64) -> f64 {
    if a > 0.0 && b > 0.0 && c > 0.0 {
        a.min(b).min(c)
    } else if a < 0.0 && b < 0.0 && c < 0.0 {
        a.max(b).max(c)
    } else {
        0.0
    }
}

/// Monotonised-central limited slope (change across one full cell) from the
/// differences to the lower and upper neighbour.
fn limited_slope(d_low: f64, d_high: f64) -> f64 {
    minmod3(2.0 * d_low, 0.5 * (d_low + d_high), 2.0 * d_high)
}

/// Face states of cell `c` on its low and high side, given its neighbours.
///
/// The free surface `h + b`, not the depth, is reconstructed: over still water its
/// slope is exactly zero, which is what keeps a lake at rest at rest. The bed is
/// reconstructed with the central slope, so the source term in [`axis_terms`] cancels
/// the face pressures exactly for a flat surface.
fn reconstruct(low: Cell, c: Cell, high: Cell, order: Order) -> (Side, Side) {
    let first = Side {
        h: c.h,
        un: c.un,
        ut: c.ut,
        z: c.b,
    };
    if order == Order::First || low.h < H_SECOND || c.h < H_SECOND || high.h < H_SECOND {
        return (first, first);
    }
    let eta = |x: Cell| x.h + x.b;
    let s_eta = limited_slope(eta(c) - eta(low), eta(high) - eta(c));
    let s_un = limited_slope(c.un - low.un, high.un - c.un);
    let s_ut = limited_slope(c.ut - low.ut, high.ut - c.ut);
    let s_z = 0.25 * (high.b - low.b);

    let side = |sign: f64| Side {
        h: eta(c) + sign * 0.5 * s_eta - (c.b + sign * s_z),
        un: c.un + sign * 0.5 * s_un,
        ut: c.ut + sign * 0.5 * s_ut,
        z: c.b + sign * s_z,
    };
    let (lo, hi) = (side(-1.0), side(1.0));
    if lo.h < H_FACE_MIN || hi.h < H_FACE_MIN {
        return (first, first);
    }
    (lo, hi)
}

fn face(l: Side, r: Side) -> Face {
    // Hydrostatic reconstruction: both sides see the same bed step, so a flat
    // free surface produces exactly balanced fluxes.
    let z_star = l.z.max(r.z);
    let hl = (l.h + l.z - z_star).max(0.0);
    let hr = (r.h + r.z - z_star).max(0.0);

    let cl = (G * hl).sqrt();
    let cr = (G * hr).sqrt();
    let s_l = (l.un - cl).min(r.un - cr);
    let s_r = (l.un + cl).max(r.un + cr);

    let f_l = [
        hl * l.un,
        hl * l.un * l.un + 0.5 * G * hl * hl,
        hl * l.un * l.ut,
    ];
    let f_r = [
        hr * r.un,
        hr * r.un * r.un + 0.5 * G * hr * hr,
        hr * r.un * r.ut,
    ];
    let u_l = [hl, hl * l.un, hl * l.ut];
    let u_r = [hr, hr * r.un, hr * r.ut];

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

/// Fluxes and source for the cell in the middle of five cells along one axis.
fn axis_terms(stencil: [Cell; 5], spacing: f64, order: Order) -> Axis {
    let [m2, m1, c, p1, p2] = stencil;
    let (_, high_of_m1) = reconstruct(m2, m1, c, order);
    let (low, high) = reconstruct(m1, c, p1, order);
    let (low_of_p1, _) = reconstruct(c, p1, p2, order);
    Axis {
        low: face(high_of_m1, low),
        high: face(high, low_of_p1),
        source: -G * 0.5 * (low.h + high.h) * (high.z - low.z) / spacing,
    }
}

/// Time derivative of the interior cells; padded arrays, ghost entries stay zero.
struct Rhs {
    h: Vec<f64>,
    hu: Vec<f64>,
    hv: Vec<f64>,
}

pub struct Solver {
    pub grid: Grid,
    pub bed: Bathymetry,
    pub cfl: f64,
    pub order: Order,
}

impl Solver {
    pub fn new(grid: Grid, bed: Bathymetry) -> Self {
        Self {
            grid,
            bed,
            cfl: 0.4,
            order: Order::Second,
        }
    }

    /// Choose the scheme order; the default from [`Solver::new`] is second order.
    pub fn with_order(mut self, order: Order) -> Self {
        self.order = order;
        self
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

    /// Advance one step of size `dt` with reflective walls.
    pub fn step(&self, s: &mut State, dt: f64) {
        s.fill_walls(&self.grid);
        match self.order {
            Order::First => {
                *s = self.advance(s, &self.rhs(s), dt);
            }
            Order::Second => {
                let mut s1 = self.advance(s, &self.rhs(s), dt);
                s1.fill_walls(&self.grid);
                let s2 = self.advance(&s1, &self.rhs(&s1), dt);
                for (a, b) in [(&mut s.h, &s2.h), (&mut s.hu, &s2.hu), (&mut s.hv, &s2.hv)] {
                    for (x, y) in a.iter_mut().zip(b) {
                        *x = 0.5 * (*x + y);
                    }
                }
                self.clean_dry(s);
            }
        }
    }

    /// Zero the momentum of dry cells and remove rounding-level negative depth.
    fn clean_dry(&self, s: &mut State) {
        for (i, j) in self.grid.interior() {
            let k = self.grid.idx(i, j);
            if s.h[k] < H_DRY {
                s.h[k] = s.h[k].max(0.0);
                s.hu[k] = 0.0;
                s.hv[k] = 0.0;
            }
        }
    }

    /// `s + dt * rhs`, interior only; ghost cells are refilled by the caller.
    fn advance(&self, s: &State, rhs: &Rhs, dt: f64) -> State {
        let mut next = s.clone();
        for (i, j) in self.grid.interior() {
            let k = self.grid.idx(i, j);
            next.h[k] += dt * rhs.h[k];
            next.hu[k] += dt * rhs.hu[k];
            next.hv[k] += dt * rhs.hv[k];
        }
        self.clean_dry(&mut next);
        next
    }

    fn rhs(&self, s: &State) -> Rhs {
        let g = &self.grid;
        let n = g.cells();
        let mut out = Rhs {
            h: vec![0.0; n],
            hu: vec![0.0; n],
            hv: vec![0.0; n],
        };
        let cell = |i: usize, j: usize, normal_is_x: bool| {
            let k = g.idx(i, j);
            let (un, ut) = if normal_is_x {
                (s.hu[k], s.hv[k])
            } else {
                (s.hv[k], s.hu[k])
            };
            Cell {
                h: s.h[k],
                un: velocity(s.h[k], un),
                ut: velocity(s.h[k], ut),
                b: self.bed.b[k],
            }
        };

        for (i, j) in g.interior() {
            let k = g.idx(i, j);
            let x = axis_terms(
                std::array::from_fn(|m| cell(i + m - 2, j, true)),
                g.dx,
                self.order,
            );
            let y = axis_terms(
                std::array::from_fn(|m| cell(i, j + m - 2, false)),
                g.dy,
                self.order,
            );

            out.h[k] =
                -(x.high.flux[0] - x.low.flux[0]) / g.dx - (y.high.flux[0] - y.low.flux[0]) / g.dy;
            out.hu[k] = -((x.high.flux[1] + x.high.corr_left) - (x.low.flux[1] + x.low.corr_right))
                / g.dx
                + x.source
                - (y.high.flux[2] - y.low.flux[2]) / g.dy;
            out.hv[k] = -((y.high.flux[1] + y.high.corr_left) - (y.low.flux[1] + y.low.corr_right))
                / g.dy
                + y.source
                - (x.high.flux[2] - x.low.flux[2]) / g.dx;
        }
        out
    }
}
