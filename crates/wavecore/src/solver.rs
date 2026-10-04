//! Well-balanced finite-volume solver for the nonlinear shallow-water equations.
//!
//! Hydrostatic reconstruction (Audusse et al. 2004) keeps still water still over any
//! bed and keeps depth non-negative at moving shorelines. The flux is HLL. Second order
//! comes from MUSCL reconstruction of the free surface and velocities with an MC
//! limiter, advanced by SSP-RK2 (Heun); [`Order::Fifth`] replaces the limiter with a
//! fifth-order reconstruction wherever the wave is smooth. Wherever a cell or a neighbour is
//! nearly dry, or reconstruction would produce a face without water, the cell drops to first
//! order.
//!
//! Every face computation is a pure function of the states on its two sides, so the
//! same kernels can be ported one-to-one to a compute shader.

use std::cell::RefCell;

use crate::bathymetry::Bathymetry;
use crate::dispersion::{Breaking, Dispersion, DispersionStats, Work, breaking_indicator, switch};
use crate::forcing::{Drive, H_FRICTION_MIN, Sponge, WaveMaker, manning_factor, relax_cell};
use crate::grid::{GHOST, Grid, wrap_ghosts_y};
use crate::par::{max_rows, rows3};
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
    /// MUSCL reconstruction with an MC limiter, SSP-RK2. Robust at shocks, but the limiter
    /// clips every smooth crest and trough to first order, which damps waves noticeably at
    /// fewer than about 40 cells per wavelength.
    Second,
    /// Unlimited fifth-order upwind-biased reconstruction where the solution is smooth (the
    /// switch `phi` from [`crate::dispersion::switch`] near 1), blended back to the MC limiter
    /// as waves steepen towards breaking. SSP-RK2 in time (SSP-RK3 was tried and made no
    /// measurable difference). Keeps the height of smooth waves on a coarse grid, which a TVD
    /// limiter cannot: a 14 s swell crossing a kilometre of shelf on 6 m cells arrives within
    /// 1.5% of the height it has on 3 m cells, where the third-order scheme it replaced lost 14%.
    Fifth,
}

/// Cell-averaged primitive variables along one axis: depth, velocity normal to the
/// axis' faces, velocity along them, and bed elevation.
#[derive(Clone, Copy)]
struct Cell {
    h: f64,
    un: f64,
    ut: f64,
    b: f64,
    /// Smoothness switch, 1 where the wave is smooth and 0 where it is steep or breaking.
    phi: f64,
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

pub(crate) fn velocity(h: f64, hm: f64) -> f64 {
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

/// Values at the low and high face of the cell `v[2]`, from the cell averages `v` of it and
/// two neighbours on each side.
///
/// `blend` is 0 for the MC limiter, which uses only the nearest neighbours, and 1 for the
/// unlimited fifth-order upwind-biased scheme
/// `(2 v0 - 13 v1 + 47 v2 + 27 v3 - 3 v4) / 60` at the high face (mirrored at the low face),
/// which is exact for quartics.
pub(crate) fn face_values(v: [f64; 5], blend: f64) -> (f64, f64) {
    let [low2, low, c, high, high2] = v;
    let s = limited_slope(c - low, high - c);
    let (mc_lo, mc_hi) = (c - 0.5 * s, c + 0.5 * s);
    if blend <= 0.0 {
        return (mc_lo, mc_hi);
    }
    let fifth_hi = (2.0 * low2 - 13.0 * low + 47.0 * c + 27.0 * high - 3.0 * high2) / 60.0;
    let fifth_lo = (2.0 * high2 - 13.0 * high + 47.0 * c + 27.0 * low - 3.0 * low2) / 60.0;
    (
        mc_lo + blend * (fifth_lo - mc_lo),
        mc_hi + blend * (fifth_hi - mc_hi),
    )
}

/// Face states of the middle cell of five on its low and high side.
///
/// The free surface `h + b`, not the depth, is reconstructed: over still water its
/// slope is exactly zero, which is what keeps a lake at rest at rest. The bed is
/// reconstructed with the central slope, so the source term in [`axis_terms`] cancels
/// the face pressures exactly for a flat surface.
fn reconstruct(cells: [Cell; 5], order: Order) -> (Side, Side) {
    let [low2, low, c, high, high2] = cells;
    let first = Side {
        h: c.h,
        un: c.un,
        ut: c.ut,
        z: c.b,
    };
    if order == Order::First || low.h < H_SECOND || c.h < H_SECOND || high.h < H_SECOND {
        return (first, first);
    }
    // Smooth only if the whole stencil is wet and smooth: a front or a shoreline within two
    // cells must keep the limiter.
    let blend = if order == Order::Fifth && low2.h >= H_SECOND && high2.h >= H_SECOND {
        cells.iter().map(|x| x.phi).fold(1.0, f64::min)
    } else {
        0.0
    };
    let (eta_lo, eta_hi) = face_values(cells.map(|x| x.h + x.b), blend);
    let (un_lo, un_hi) = face_values(cells.map(|x| x.un), blend);
    let (ut_lo, ut_hi) = face_values(cells.map(|x| x.ut), blend);
    let s_z = 0.25 * (high.b - low.b);

    let side = |sign: f64, eta_f: f64, un: f64, ut: f64| Side {
        h: eta_f - (c.b + sign * s_z),
        un,
        ut,
        z: c.b + sign * s_z,
    };
    let (lo, hi) = (
        side(-1.0, eta_lo, un_lo, ut_lo),
        side(1.0, eta_hi, un_hi, ut_hi),
    );
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

/// Fluxes and source for the cell in the middle of seven cells along one axis.
fn axis_terms(stencil: [Cell; 7], spacing: f64, order: Order) -> Axis {
    let five = |first: usize| std::array::from_fn(|m| stencil[first + m]);
    let (_, high_of_m1) = reconstruct(five(0), order);
    let (low, high) = reconstruct(five(1), order);
    let (low_of_p1, _) = reconstruct(five(2), order);
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
    /// Manning roughness `n` in s/m^(1/3); zero switches friction off.
    pub manning: f64,
    /// The south and north edges wrap around instead of reflecting.
    pub periodic_y: bool,
    wavemaker_spec: Option<WaveMaker>,
    wavemaker: Option<Drive>,
    sponges: Vec<Sponge>,
    dispersion: Option<(Dispersion, RefCell<Work>)>,
    /// The smoothness switch, recomputed on every evaluation.
    phi: RefCell<Vec<f64>>,
}

impl Solver {
    pub fn new(grid: Grid, bed: Bathymetry) -> Self {
        Self {
            grid,
            bed,
            cfl: 0.4,
            order: Order::Second,
            manning: 0.0,
            periodic_y: false,
            wavemaker_spec: None,
            wavemaker: None,
            sponges: Vec::new(),
            dispersion: None,
            phi: RefCell::new(vec![0.0; grid.cells()]),
        }
    }

    /// Choose the scheme order; the default from [`Solver::new`] is second order.
    pub fn with_order(mut self, order: Order) -> Self {
        self.order = order;
        self
    }

    /// Add Manning bottom friction with roughness `n` in s/m^(1/3).
    pub fn with_manning(mut self, n: f64) -> Self {
        self.manning = n;
        self
    }

    /// Add an absorbing sponge layer. Call again for layers of different widths.
    pub fn with_sponge(mut self, sponge: Sponge) -> Self {
        self.sponges.push(sponge);
        self
    }

    /// Add an internal wave-maker. The still-water depth at the source line is read
    /// from the bed, so the bed must be set first and lie below `y = 0` there.
    pub fn with_wavemaker(mut self, spec: WaveMaker) -> Self {
        self.wavemaker_spec = Some(spec);
        self.rebuild_drive();
        self
    }

    /// Replace shallow water with the Serre-Green-Naghdi equations: waves of different
    /// lengths travel at different speeds, so big swells no longer steepen into shocks in
    /// deep water, and dispersion fades out where waves break. See [`crate::dispersion`].
    pub fn with_dispersion(mut self, dispersion: Dispersion) -> Self {
        self.dispersion = Some((dispersion, RefCell::new(Work::new(&self.grid))));
        self.rebuild_drive();
        self
    }

    /// How close each interior cell is to breaking, from 0 (smooth) to 1 (steep or tall for
    /// its depth), row-major with `y` outer. These are the criteria that switch the dispersive
    /// terms off, so with dispersion this is where the model treats the wave as breaking; for
    /// shallow water it marks where the same criteria hold. It is 0 where there is no water
    /// deep enough to carry a wave.
    pub fn breaking_indicator(&self, s: &State) -> Vec<f64> {
        breaking_indicator(
            &self.grid,
            &self.bed,
            s,
            self.breaking_config(),
            self.periodic_y,
        )
    }

    /// How hard the dispersive solves have worked so far, or `None` for shallow water.
    pub fn dispersion_stats(&self) -> Option<DispersionStats> {
        self.dispersion
            .as_ref()
            .map(|(_, work)| work.borrow().stats)
    }

    /// The wave-maker's strength depends on the linear wave relation, so it is worked out
    /// again whenever the bed, the wave-maker or the dispersion changes.
    fn rebuild_drive(&mut self) {
        let Some(spec) = self.wavemaker_spec else {
            return;
        };
        let g = &self.grid;
        let i = ((spec.x / g.dx - 0.5).round().max(0.0) as usize).min(g.nx - 1) + GHOST;
        let depth = -self.bed.b[g.idx(i, GHOST + g.ny / 2)];
        self.wavemaker = Some(Drive::new(spec, g, depth, self.dispersion.is_some()));
        self.check_periodic();
    }

    /// Make the domain periodic in `y`. Set the bed first; its ghost rows are wrapped.
    /// An oblique wave-maker needs `k_y * Ly` to be a whole multiple of 2 pi so the
    /// wave fits the period; this panics otherwise.
    pub fn with_periodic_y(mut self) -> Self {
        self.periodic_y = true;
        wrap_ghosts_y(&self.grid, &mut self.bed.b);
        self.check_periodic();
        self
    }

    fn check_periodic(&self) {
        if let (true, Some(drive)) = (self.periodic_y, &self.wavemaker) {
            let turns = drive.ky() * self.grid.ny as f64 * self.grid.dy / std::f64::consts::TAU;
            assert!(
                (turns - turns.round()).abs() < 1e-6,
                "wave-maker does not fit the periodic domain: k_y * Ly / 2pi = {turns}, must be a whole number"
            );
        }
    }

    /// Largest stable time step for the current state. Infinite if nothing can move.
    pub fn stable_dt(&self, s: &State) -> f64 {
        let g = &self.grid;
        let w = g.width();
        let rate = max_rows(GHOST, GHOST + g.ny, w, |j| {
            let mut rate: f64 = 0.0;
            for k in j * w + GHOST..j * w + GHOST + g.nx {
                let c = (G * s.h[k]).sqrt();
                let u = velocity(s.h[k], s.hu[k]).abs();
                let v = velocity(s.h[k], s.hv[k]).abs();
                rate = rate.max((u + c) / g.dx + (v + c) / g.dy);
            }
            rate
        });
        if rate > 0.0 {
            self.cfl / rate
        } else {
            f64::INFINITY
        }
    }

    /// Advance one step of size `dt` with reflective walls, then apply sponge damping
    /// and friction. Advances `s.time` by `dt`.
    pub fn step(&self, s: &mut State, dt: f64) {
        let t = s.time;
        s.fill_boundaries(&self.grid, self.periodic_y);
        match self.order {
            Order::First => {
                *s = self.advance(s, &self.rhs(s, t), dt);
            }
            Order::Second | Order::Fifth => {
                let mut s1 = self.advance(s, &self.rhs(s, t), dt);
                s1.fill_boundaries(&self.grid, self.periodic_y);
                let s2 = self.advance(&s1, &self.rhs(&s1, t + dt), dt);
                let w = self.grid.width();
                rows3(w, &mut s.h, &mut s.hu, &mut s.hv, |j, rh, rhu, rhv| {
                    for i in 0..w {
                        let k = j * w + i;
                        rh[i] = 0.5 * (rh[i] + s2.h[k]);
                        rhu[i] = 0.5 * (rhu[i] + s2.hu[k]);
                        rhv[i] = 0.5 * (rhv[i] + s2.hv[k]);
                    }
                });
                self.clean_dry(s);
            }
        }
        self.damp(s, dt);
        s.time = t + dt;
    }

    /// Sponge layer and bottom friction, applied after the flux update.
    fn damp(&self, s: &mut State, dt: f64) {
        if self.sponges.is_empty() && self.manning == 0.0 {
            return;
        }
        let g = self.grid;
        let (sponges, manning, bed) = (&self.sponges, self.manning, &self.bed.b);
        rows3(
            g.width(),
            &mut s.h,
            &mut s.hu,
            &mut s.hv,
            |j, rh, rhu, rhv| {
                if j < GHOST || j >= g.ny + GHOST {
                    return;
                }
                for i in GHOST..g.nx + GHOST {
                    let k = j * g.width() + i;
                    for sp in sponges {
                        let rate = sp.rate(&g, i, j);
                        if rate > 0.0 {
                            // Relax towards still water, depth and momentum at the same rate.
                            // Relaxing the surface towards the level instead drained swash on
                            // land to a film microns deep that kept its momentum and moved at
                            // tens of metres per second, which set the time step for the grid.
                            let f = (-rate * dt).exp();
                            let cell = (&mut rh[i], &mut rhu[i], &mut rhv[i]);
                            relax_cell(cell, bed[k], (sp.level, 0.0, 0.0), f);
                        }
                    }
                    if manning > 0.0 && rh[i] > H_FRICTION_MIN {
                        let speed = rhu[i].hypot(rhv[i]) / rh[i];
                        let f = manning_factor(rh[i], speed, dt, manning);
                        rhu[i] *= f;
                        rhv[i] *= f;
                    }
                }
            },
        );
        self.clean_dry(s);
    }

    /// Zero the momentum of dry cells and remove rounding-level negative depth.
    fn clean_dry(&self, s: &mut State) {
        let (w, nx, ny) = (self.grid.width(), self.grid.nx, self.grid.ny);
        rows3(w, &mut s.h, &mut s.hu, &mut s.hv, |j, rh, rhu, rhv| {
            if j < GHOST || j >= ny + GHOST {
                return;
            }
            for i in GHOST..nx + GHOST {
                if rh[i] < H_DRY {
                    rh[i] = rh[i].max(0.0);
                    rhu[i] = 0.0;
                    rhv[i] = 0.0;
                }
            }
        });
    }

    /// `s + dt * rhs`, interior only; ghost cells are refilled by the caller.
    fn advance(&self, s: &State, rhs: &Rhs, dt: f64) -> State {
        let mut next = s.clone();
        let (w, nx, ny) = (self.grid.width(), self.grid.nx, self.grid.ny);
        rows3(
            w,
            &mut next.h,
            &mut next.hu,
            &mut next.hv,
            |j, rh, rhu, rhv| {
                if j < GHOST || j >= ny + GHOST {
                    return;
                }
                for i in GHOST..nx + GHOST {
                    let k = j * w + i;
                    rh[i] += dt * rhs.h[k];
                    rhu[i] += dt * rhs.hu[k];
                    rhv[i] += dt * rhs.hv[k];
                }
            },
        );
        self.clean_dry(&mut next);
        next
    }

    /// The breaking thresholds in force: the dispersive model's, or the defaults.
    fn breaking_config(&self) -> Breaking {
        self.dispersion
            .as_ref()
            .map_or_else(Breaking::default, |(d, _)| d.breaking)
    }

    fn rhs(&self, s: &State, t: f64) -> Rhs {
        let g = &self.grid;
        let n = g.cells();
        // How smooth the wave is, everywhere. The fifth-order reconstruction and the
        // dispersive terms both fade out where it is not.
        let mut phi_guard = self.phi.borrow_mut();
        if self.order == Order::Fifth || self.dispersion.is_some() {
            switch(g, &self.bed, s, self.breaking_config(), &mut phi_guard);
        }
        let phi: &[f64] = &phi_guard;
        let mut out = Rhs {
            h: vec![0.0; n],
            hu: vec![0.0; n],
            hv: vec![0.0; n],
        };
        let bed = &self.bed.b[..];
        let order = self.order;
        let drive = self.wavemaker.as_ref();
        let g = *g;
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
                b: bed[k],
                phi: phi[k],
            }
        };

        rows3(
            g.width(),
            &mut out.h,
            &mut out.hu,
            &mut out.hv,
            |j, rh, rhu, rhv| {
                if j < GHOST || j >= g.ny + GHOST {
                    return;
                }
                for i in GHOST..g.nx + GHOST {
                    let x = axis_terms(
                        std::array::from_fn(|m| cell(i + m - 3, j, true)),
                        g.dx,
                        order,
                    );
                    let y = axis_terms(
                        std::array::from_fn(|m| cell(i, j + m - 3, false)),
                        g.dy,
                        order,
                    );

                    rh[i] = -(x.high.flux[0] - x.low.flux[0]) / g.dx
                        - (y.high.flux[0] - y.low.flux[0]) / g.dy;
                    rhu[i] = -((x.high.flux[1] + x.high.corr_left)
                        - (x.low.flux[1] + x.low.corr_right))
                        / g.dx
                        + x.source
                        - (y.high.flux[2] - y.low.flux[2]) / g.dy;
                    rhv[i] = -((y.high.flux[1] + y.high.corr_left)
                        - (y.low.flux[1] + y.low.corr_right))
                        / g.dy
                        + y.source
                        - (x.high.flux[2] - x.low.flux[2]) / g.dx;
                    if let Some(drive) = drive {
                        let (px, py) = g.centre(i, j);
                        rh[i] += drive.source(px, py, t);
                    }
                }
            },
        );
        let g = &self.grid;
        if let Some((dispersion, work)) = &self.dispersion {
            dispersion.apply(
                g,
                self.periodic_y,
                s,
                &out.h,
                &mut out.hu,
                &mut out.hv,
                phi,
                &mut work.borrow_mut(),
            );
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fifth_order_face_values_are_exact_for_quartics() {
        // Cell averages of a quartic over unit cells centred on -2..2 give its exact value at
        // the faces of the middle cell, x = -1/2 and 1/2.
        let p = |x: f64| 1.0 + 2.0 * x - x * x + 0.5 * x.powi(3) - 0.25 * x.powi(4);
        let integral = |x: f64| x + x * x - x.powi(3) / 3.0 + 0.125 * x.powi(4) - 0.05 * x.powi(5);
        let v: [f64; 5] =
            std::array::from_fn(|m| integral(m as f64 - 1.5) - integral(m as f64 - 2.5));
        let (lo, hi) = face_values(v, 1.0);
        assert!((lo - p(-0.5)).abs() < 1e-12, "{lo} vs {}", p(-0.5));
        assert!((hi - p(0.5)).abs() < 1e-12, "{hi} vs {}", p(0.5));
        // With the limiter fully on, the faces stay between the neighbouring averages.
        let (lo, hi) = face_values([0.0, 0.0, 1.0, 1.0, 1.0], 0.0);
        assert!((0.0..=1.0).contains(&lo) && (0.0..=1.0).contains(&hi));
    }
}
