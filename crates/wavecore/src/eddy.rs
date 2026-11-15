//! Breaking fronts, followed from step to step: where they are, the eddy viscosity of their
//! turbulence (Kennedy et al. 2000, the closure FUNWAVE uses), and the broken water behind them
//! that the switch keeps without dispersion. Opt-in, through [`crate::Breaking::viscosity`].
//!
//! The switch alone reads the surface afresh at every step. It turns the dispersive terms off
//! at the steep front of a wave about to break, and the shock then dissipates it; but once the
//! front has broken and eased, the switch turns them back on and the wave goes on as a smooth,
//! steep bore that loses little. On Ting & Kirby's spilling beach the inner surf zone then
//! stays at H/d 1.0, where the laboratory measured 0.6. Plain shallow water over the inner surf
//! zone gives 0.65 to 0.70: the shock dissipates the right amount once it is left alone.
//!
//! So a breaking front is followed. It starts where the switch, read afresh, says the wave
//! breaks and the surface is rising, and goes on breaking, cell to cell, by Kennedy et al.'s
//! criterion: while the surface rises faster than a threshold that falls from
//! `initial * sqrt(g h)` to `settled * sqrt(g h)` over `duration * sqrt(h / g)`, the rate of
//! rise being `-div(hu)` from the mass equation. So a bore keeps breaking in the inner surf zone
//! although its face is far less steep than at the break point. Where a front breaks the switch
//! reads breaking, and so it does for `linger * sqrt(h / g)` after the front has passed, over
//! the roller behind it. The front also gets Kennedy's eddy viscosity in the momentum
//! equations:
//!
//! ```text
//! d(hu)/dt += d/dx( nu d(hu)/dx ) + 1/2 d/dy( nu (d(hu)/dy + d(hv)/dx) )
//! d(hv)/dt += d/dy( nu d(hv)/dy ) + 1/2 d/dx( nu (d(hv)/dx + d(hu)/dy) )
//! nu = B delta^2 h eta_t
//! ```
//!
//! where `B` ramps from 0 to 1 as the rate of rise `eta_t` goes from the threshold to twice it.
//! It is applied after each step, split from the rest, in as many explicit sub-steps as it needs
//! to stay stable. Depths here are the total depth. Only with dispersion: in plain shallow
//! water the shock alone is the breaking.
//!
//! On Ting & Kirby's beach (`tests/lab.rs`) the inner surf zone is 1.66 times the measured
//! height without this, and 1.04 with it. Kennedy's viscosity with the switch left as it was
//! gives only 1.53. With fronts started by Kennedy's own criterion (the surface first rising
//! faster than `initial * sqrt(g h)`): the fronts with the switch off behind them and no
//! viscosity gave 1.13, with the viscosity 1.06; lingering 5 or 10 `sqrt(h / g)` gave the same,
//! 20 gave 1.10, and none 1.47.
//!
//! **Why it is not on by default.** With it, the plunging breaker is tallest 0.86 m before the
//! laboratory's break point, against 0.5 m allowed. On that steep 5 s wave the switch reads
//! breaking at x = 6.93 m however a front is started (at 0.8, fully off, or by Kennedy's rate of
//! rise), and the memory holds on to that reading where the switch alone lets the crest recover
//! and grow on to 7.3 m. Leaving a front to the switch until it has broken for `duration` moved
//! the plunging break only to 7.08 m, and the spilling breaker then broke 0.8 m late and its
//! inner surf zone came out 1.26 times the laboratory's.

use crate::dispersion::H_DISPERSION_MIN;
use crate::grid::{GHOST, Grid};
use crate::par::{max_rows, rows1, rows2};
use crate::solver::G;
use crate::state::State;

/// Water thinner than this has no breaking front and gets no eddy viscosity, and none passes
/// into it: diffusing momentum into a film a few millimetres deep gives it an absurd velocity.
const H_BREAKING_MIN: f64 = 0.01;

/// Where the switch, `1 - phi`, counts a cell as breaking: the reading a run reports.
const BREAKING_ONSET: f64 = 0.8;

/// How breaking fronts are followed, Kennedy et al.'s (2000) coefficients by default.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EddyViscosity {
    /// A front that has just broken keeps breaking while the surface rises faster than this
    /// times `sqrt(g h)`.
    pub initial: f64,
    /// ...and a front that has broken for [`EddyViscosity::duration`] keeps breaking while it
    /// rises faster than this times `sqrt(g h)`.
    pub settled: f64,
    /// How long the threshold takes to fall from `initial` to `settled`, in `sqrt(h / g)`.
    pub duration: f64,
    /// The mixing length as a multiple of the depth, `delta`. 0 follows the fronts without
    /// any eddy viscosity.
    pub mixing: f64,
    /// How long a cell a breaking front has passed stays without dispersion, in `sqrt(h / g)`.
    pub linger: f64,
}

impl Default for EddyViscosity {
    fn default() -> Self {
        Self {
            initial: 0.65,
            settled: 0.15,
            duration: 5.0,
            mixing: 1.2,
            linger: 5.0,
        }
    }
}

/// Scratch arrays, padded like the grid. Allocated once.
pub(crate) struct EddyWork {
    rise: Vec<f64>,
    nu: Vec<f64>,
    memory: Vec<f64>,
    hu: Vec<f64>,
    hv: Vec<f64>,
}

impl EddyWork {
    pub(crate) fn new(grid: &Grid) -> Self {
        let n = grid.cells();
        Self {
            rise: vec![0.0; n],
            nu: vec![0.0; n],
            memory: vec![0.0; n],
            hu: vec![0.0; n],
            hv: vec![0.0; n],
        }
    }
}

impl EddyViscosity {
    /// Kennedy et al.'s transition, `duration * sqrt(h / g)` seconds in water `h` deep.
    fn span(&self, h: f64) -> f64 {
        self.duration * (h / G).sqrt()
    }

    /// Whether the switch stays off in a cell with breaking memory `memory` (see
    /// [`State::breaking`]): on a breaking front, and in the water a front has just left.
    pub(crate) fn keeps_off(memory: f64) -> bool {
        memory != 0.0
    }

    /// The rate of rise, as a multiple of `sqrt(g h)`, above which a front that has been
    /// breaking for `age` seconds in water `h` deep goes on breaking.
    fn threshold(&self, age: f64, h: f64) -> f64 {
        let span = self.span(h);
        if age >= span {
            self.settled
        } else {
            self.initial + (self.settled - self.initial) * age / span
        }
    }

    /// Follow the breaking fronts in `s.breaking` on to the state just reached, then apply
    /// their eddy viscosity over `dt`. `phi` is the switch read afresh from `s`, without the
    /// fronts. `s` must have its ghost cells filled; `fill` refills them between the viscous
    /// sub-steps.
    pub(crate) fn step(
        &self,
        grid: &Grid,
        s: &mut State,
        phi: &[f64],
        dt: f64,
        work: &mut EddyWork,
        fill: impl FnMut(&mut State),
    ) {
        self.follow(grid, s, phi, dt, work);
        if self.mixing > 0.0 {
            diffuse(grid, s, dt, work, fill);
        }
    }

    /// Where the fronts are now, into `s.breaking`, and their eddy viscosity into `work.nu`.
    ///
    /// A cell breaks if the switch `phi` reads breaking there while the surface rises, or if the
    /// surface rises faster than the threshold for the oldest breaking front among it and its
    /// four neighbours. It then takes that front's age plus `dt` (the front moves less than a
    /// cell per step, the time step being CFL-limited, so the neighbours carry it along), or
    /// starts a front of its own. A cell a bore has just left starts to count down `linger`,
    /// and a cell counting down goes on until it reaches 0.
    fn follow(&self, grid: &Grid, s: &mut State, phi: &[f64], dt: f64, work: &mut EddyWork) {
        let w = grid.width();
        let (nx, ny) = (grid.nx, grid.ny);
        let (ix, iy) = (0.5 / grid.dx, 0.5 / grid.dy);
        let (h, hu, hv) = (&s.h[..], &s.hu[..], &s.hv[..]);
        rows1(w, &mut work.rise, |j, row| {
            row.fill(0.0);
            if j < GHOST || j >= ny + GHOST {
                return;
            }
            for (i, out) in row.iter_mut().enumerate().take(nx + GHOST).skip(GHOST) {
                let k = j * w + i;
                *out = -((hu[k + 1] - hu[k - 1]) * ix + (hv[k + w] - hv[k - w]) * iy);
            }
        });
        let rise: &[f64] = &work.rise;
        let old = &s.breaking[..];
        let this = *self;
        rows2(w, &mut work.memory, &mut work.nu, |j, rm, rn| {
            rm.fill(0.0);
            rn.fill(0.0);
            if j < GHOST || j >= ny + GHOST {
                return;
            }
            for i in GHOST..nx + GHOST {
                let k = j * w + i;
                let depth = h[k];
                if depth < H_BREAKING_MIN {
                    continue;
                }
                let c = (G * depth).sqrt();
                let front = [k, k - 1, k + 1, k - w, k + w]
                    .iter()
                    .map(|&m| old[m])
                    .fold(0.0, f64::max);
                let starts =
                    depth >= H_DISPERSION_MIN && 1.0 - phi[k] >= BREAKING_ONSET && rise[k] > 0.0;
                let goes_on = front > 0.0 && rise[k] > this.threshold(front, depth) * c;
                let age = if front > 0.0 && (starts || goes_on) {
                    front + dt
                } else if starts {
                    dt
                } else {
                    0.0
                };
                if age > 0.0 {
                    let b = (rise[k] / (this.threshold(age, depth) * c) - 1.0).clamp(0.0, 1.0);
                    rm[i] = age;
                    rn[i] = b * this.mixing * this.mixing * depth * rise[k];
                } else if old[k] > 0.0 {
                    rm[i] = -this.linger * (depth / G).sqrt();
                } else if old[k] < 0.0 {
                    rm[i] = (old[k] + dt).min(0.0);
                }
            }
        });
        std::mem::swap(&mut s.breaking, &mut work.memory);
    }
}

/// Apply the eddy viscosity `work.nu` (zero outside the interior) to the momentum of `s` over
/// `dt`, in explicit sub-steps short enough to be stable. Returns how many it took.
fn diffuse(
    grid: &Grid,
    s: &mut State,
    dt: f64,
    work: &mut EddyWork,
    mut fill: impl FnMut(&mut State),
) -> usize {
    let w = grid.width();
    let (nx, ny) = (grid.nx, grid.ny);
    let nu: &[f64] = &work.nu;
    let top = max_rows(GHOST, GHOST + ny, w, |j| {
        nu[j * w + GHOST..j * w + GHOST + nx]
            .iter()
            .fold(0.0, |a: f64, &b| a.max(b))
    });
    if top <= 0.0 {
        return 0;
    }
    let (idx2, idy2) = (1.0 / (grid.dx * grid.dx), 1.0 / (grid.dy * grid.dy));
    // Forward Euler on d/dx(nu d/dx) + 1/2 d/dy(nu d/dy) is stable for
    // dt nu (2/dx^2 + 1/dy^2) <= 1; half that leaves room for the cross terms.
    let limit = 0.5 / (top * (2.0 * idx2 + idy2).max(idx2 + 2.0 * idy2));
    let steps = (dt / limit).ceil().max(1.0) as usize;
    let sub = dt / steps as f64;
    let cross = 0.25 / (grid.dx * grid.dy);
    for step in 0..steps {
        if step > 0 {
            fill(s);
        }
        let (h, hu, hv) = (&s.h[..], &s.hu[..], &s.hv[..]);
        // Viscosity between two cells: none if either is too thin to take part.
        let link = |a: usize, b: usize| {
            if h[a] < H_BREAKING_MIN || h[b] < H_BREAKING_MIN {
                0.0
            } else {
                0.5 * (nu[a] + nu[b])
            }
        };
        let wet = |a: usize| if h[a] < H_BREAKING_MIN { 0.0 } else { nu[a] };
        rows2(w, &mut work.hu, &mut work.hv, |j, ru, rv| {
            if j < GHOST || j >= ny + GHOST {
                return;
            }
            for i in GHOST..nx + GHOST {
                let k = j * w + i;
                if h[k] < H_BREAKING_MIN {
                    ru[i] = hu[k];
                    rv[i] = hv[k];
                    continue;
                }
                let (e, west, north, south) = (k + 1, k - 1, k + w, k - w);
                // d/dx(nu d/dx q) and d/dy(nu d/dy q), conservative between cells.
                let xx = |q: &[f64]| {
                    (link(k, e) * (q[e] - q[k]) - link(k, west) * (q[k] - q[west])) * idx2
                };
                let yy = |q: &[f64]| {
                    (link(k, north) * (q[north] - q[k]) - link(k, south) * (q[k] - q[south])) * idy2
                };
                // d/dy(nu d/dx q) and d/dx(nu d/dy q), centred.
                let yx = |q: &[f64]| {
                    (wet(north) * (q[north + 1] - q[north - 1])
                        - wet(south) * (q[south + 1] - q[south - 1]))
                        * cross
                };
                let xy = |q: &[f64]| {
                    (wet(e) * (q[e + w] - q[e - w]) - wet(west) * (q[west + w] - q[west - w]))
                        * cross
                };
                ru[i] = hu[k] + sub * (xx(hu) + 0.5 * (yy(hu) + yx(hv)));
                rv[i] = hv[k] + sub * (yy(hv) + 0.5 * (xx(hv) + xy(hu)));
            }
        });
        // Only the interior was written; the ghost cells are refilled before they are read.
        std::mem::swap(&mut s.hu, &mut work.hu);
        std::mem::swap(&mut s.hv, &mut work.hv);
    }
    steps
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_threshold_falls_from_onset_to_settled_over_the_duration() {
        let e = EddyViscosity::default();
        let h = 1.0;
        let span = 5.0 * (h / G).sqrt();
        assert_eq!(e.threshold(0.0, h), 0.65);
        assert!((e.threshold(0.5 * span, h) - 0.40).abs() < 1e-12);
        assert_eq!(e.threshold(2.0 * span, h), 0.15);
    }

    /// Still water 1 m deep on cells of 0.1 m, with `hu` set by `momentum(i)` along x.
    fn channel(grid: &Grid, breaking: &[f64], momentum: impl Fn(usize) -> f64) -> State {
        let mut s = State::zeros(grid);
        for (i, j) in grid.interior() {
            let k = grid.idx(i, j);
            s.h[k] = 1.0;
            s.hu[k] = momentum(i);
        }
        s.breaking.copy_from_slice(breaking);
        s.fill_walls(grid);
        s
    }

    #[test]
    fn diffusion_conserves_momentum_and_smooths_a_jump() {
        // A jump in hu along x spreads out and keeps its total; the viscosity stays clear of
        // the walls, which take momentum out.
        let grid = Grid::new(40, 3, 0.1, 0.1);
        let none = vec![0.0; grid.cells()];
        let mut s = channel(&grid, &none, |i| if i < GHOST + 20 { 0.5 } else { 0.0 });
        let mut work = EddyWork::new(&grid);
        for (i, j) in grid.interior() {
            if (GHOST + 8..GHOST + 32).contains(&i) {
                work.nu[grid.idx(i, j)] = 0.05;
            }
        }
        let total = |s: &State| {
            grid.interior()
                .map(|(i, j)| s.hu[grid.idx(i, j)])
                .sum::<f64>()
        };
        let before = total(&s);
        let steps = diffuse(&grid, &mut s, 0.5, &mut work, |s| s.fill_walls(&grid));
        assert!(steps > 1, "the sub-steps were not needed: {steps}");
        assert!(
            (total(&s) - before).abs() < 1e-9 * before,
            "{} vs {before}",
            total(&s)
        );
        let j = GHOST + 1;
        let (a, b) = (s.hu[grid.idx(GHOST + 19, j)], s.hu[grid.idx(GHOST + 20, j)]);
        assert!(a < 0.5 && b > 0.0 && a > b, "not smoothed: {a} {b}");
        assert!(
            grid.interior()
                .all(|(i, j)| (-1e-12..=0.5 + 1e-12).contains(&s.hu[grid.idx(i, j)])),
            "diffusion overshot"
        );
    }

    #[test]
    fn a_front_starts_where_the_switch_says_moves_on_with_its_age_and_lingers_behind() {
        let grid = Grid::new(20, 3, 0.1, 0.1);
        let e = EddyViscosity::default();
        let (dt, c) = (0.01, G.sqrt());
        let cell = |i: usize| grid.idx(GHOST + i, GHOST + 1);
        // Momentum `a` towards cell `at` from both sides: it rises at a / dx, its neighbours
        // at half that.
        let towards = |at: usize, rate: f64| {
            move |i: usize| {
                let a = rate * c * 0.1;
                match i.cmp(&(GHOST + at)) {
                    std::cmp::Ordering::Less => a,
                    std::cmp::Ordering::Equal => 0.0,
                    std::cmp::Ordering::Greater => -a,
                }
            }
        };
        let smooth = vec![1.0; grid.cells()];
        let mut work = EddyWork::new(&grid);
        let none = vec![0.0; grid.cells()];
        // Rising fast but smooth by the switch: no front.
        let mut s = channel(&grid, &none, towards(5, 0.8));
        e.follow(&grid, &mut s, &smooth, dt, &mut work);
        assert!(
            s.breaking.iter().all(|&m| m == 0.0),
            "a front without the switch"
        );
        // The switch reads breaking at cells 4 to 6, but only 4 to 6 rise: fronts there.
        let mut broken = smooth.clone();
        for i in 3..=7 {
            broken[cell(i)] = 0.0;
        }
        let mut s = channel(&grid, &none, towards(5, 0.8));
        e.follow(&grid, &mut s, &broken, dt, &mut work);
        assert_eq!(s.breaking[cell(5)], dt);
        assert_eq!(s.breaking[cell(4)], dt);
        assert_eq!(
            s.breaking[cell(3)],
            0.0,
            "a front where the surface does not rise"
        );
        assert!(work.nu[cell(5)] > 0.0);
        assert_eq!(work.nu[cell(4)], 0.0, "rising at 0.4 has a viscosity");
        // Kennedy's threshold falls with a front's age: a front broken for longer than the
        // transition goes on where the surface rises at 0.35 sqrt(g h). A step on, smooth by
        // the switch, the rise has moved a cell: the front goes with it and keeps its age, and
        // the cell it has left behind starts to linger, with the switch off.
        let old = 2.0 * e.span(1.0);
        let mut bore = none.clone();
        for i in 4..=6 {
            bore[cell(i)] = old;
        }
        let mut s = channel(&grid, &bore, towards(6, 0.7));
        e.follow(&grid, &mut s, &smooth, dt, &mut work);
        assert_eq!(s.breaking[cell(6)], old + dt);
        assert_eq!(
            s.breaking[cell(5)],
            old + dt,
            "a bore rising at 0.35 stopped"
        );
        let linger = e.linger * (1.0 / G).sqrt();
        assert_eq!(s.breaking[cell(4)], -linger);
        assert!(EddyViscosity::keeps_off(old) && EddyViscosity::keeps_off(-linger));
        assert!(!EddyViscosity::keeps_off(0.0));
        // In still water the fronts stop and linger too, the lingering counts down, and
        // nothing has a viscosity.
        let mut s = channel(&grid, &s.breaking, |_| 0.0);
        e.follow(&grid, &mut s, &smooth, dt, &mut work);
        assert_eq!(s.breaking[cell(6)], -linger);
        assert!((s.breaking[cell(4)] + linger - dt).abs() < 1e-15);
        assert!(
            work.nu.iter().all(|&v| v == 0.0),
            "still water has a viscosity"
        );
    }
}
