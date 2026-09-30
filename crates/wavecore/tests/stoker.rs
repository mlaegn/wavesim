//! Verification ladder, rung 2: dam breaks against their exact solutions.
//!
//! Wet bed: the Stoker solution (rarefaction, constant middle state, shock).
//! Dry bed: the Ritter solution (rarefaction into a dry channel).

use wavecore::{Bathymetry, G, GHOST, Grid, Order, Solver, State};

const LENGTH: f64 = 100.0;
const DAM: f64 = 50.0;

/// Exact solution of a 1D dam break as a function of `xi = (x - DAM) / t`.
struct Exact {
    h_l: f64,
    h_r: f64,
    h_m: f64,
    u_m: f64,
    shock: f64,
}

impl Exact {
    fn new(h_l: f64, h_r: f64) -> Self {
        let c_l = (G * h_l).sqrt();
        if h_r <= 0.0 {
            return Self {
                h_l,
                h_r,
                h_m: 0.0,
                u_m: 0.0,
                shock: 2.0 * c_l,
            };
        }
        // Middle state: rarefaction on the left meets a shock on the right.
        let f = |h_m: f64| {
            2.0 * (c_l - (G * h_m).sqrt())
                - (h_m - h_r) * (0.5 * G * (h_m + h_r) / (h_m * h_r)).sqrt()
        };
        let (mut lo, mut hi) = (h_r, h_l);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if f(mid) > 0.0 { lo = mid } else { hi = mid }
        }
        let h_m = 0.5 * (lo + hi);
        let u_m = 2.0 * (c_l - (G * h_m).sqrt());
        Self {
            h_l,
            h_r,
            h_m,
            u_m,
            shock: h_m * u_m / (h_m - h_r),
        }
    }

    fn depth(&self, xi: f64) -> f64 {
        let c_l = (G * self.h_l).sqrt();
        if xi < -c_l {
            self.h_l
        } else if self.h_r <= 0.0 {
            if xi < 2.0 * c_l {
                (2.0 * c_l - xi).powi(2) / (9.0 * G)
            } else {
                0.0
            }
        } else if xi <= self.u_m - (G * self.h_m).sqrt() {
            (2.0 * c_l - xi).powi(2) / (9.0 * G)
        } else if xi < self.shock {
            self.h_m
        } else {
            self.h_r
        }
    }
}

/// L1 depth error in m^2 (depth error times length) against the exact solution.
fn dam_break_error(nx: usize, order: Order, h_l: f64, h_r: f64, t_end: f64) -> f64 {
    let dx = LENGTH / nx as f64;
    let grid = Grid::new(nx, 4, dx, dx);
    let bed = Bathymetry::flat(&grid, 0.0);
    let mut s = State::zeros(&grid);
    for (i, j) in grid.interior() {
        let (x, _) = grid.centre(i, j);
        s.h[grid.idx(i, j)] = if x < DAM { h_l } else { h_r };
    }
    let solver = Solver::new(grid, bed).with_order(order);

    let mut t = 0.0;
    while t < t_end {
        let dt = solver.stable_dt(&s).min(t_end - t);
        solver.step(&mut s, dt);
        t += dt;
    }

    let exact = Exact::new(h_l, h_r);
    let j = GHOST + 2;
    let mut err = 0.0;
    for i in GHOST..GHOST + nx {
        let (x, _) = grid.centre(i, j);
        err += (s.h[grid.idx(i, j)] - exact.depth((x - DAM) / t_end)).abs() * dx;
    }
    err
}

#[test]
fn exact_wet_bed_solution_is_self_consistent() {
    let e = Exact::new(2.0, 0.5);
    assert!(e.h_r < e.h_m && e.h_m < e.h_l);
    // Mass flux across the shock: s (h_m - h_r) = h_m u_m.
    assert!((e.shock * (e.h_m - e.h_r) - e.h_m * e.u_m).abs() < 1e-12);
    // Momentum jump across the shock (Rankine-Hugoniot).
    let lhs = e.shock * e.h_m * e.u_m;
    let rhs = e.h_m * e.u_m * e.u_m + 0.5 * G * (e.h_m.powi(2) - e.h_r.powi(2));
    assert!((lhs - rhs).abs() < 1e-9, "{lhs} vs {rhs}");
}

/// Second order must be much more accurate than first order and must converge at
/// close to the rate 1 that is the limit for L1 error across a discontinuity.
fn check_dam_break(label: &str, h_l: f64, h_r: f64) {
    let t = 4.0;
    let first = dam_break_error(400, Order::First, h_l, h_r, t);
    let coarse = dam_break_error(200, Order::Second, h_l, h_r, t);
    let fine = dam_break_error(400, Order::Second, h_l, h_r, t);
    let rate = (coarse / fine).log2();
    println!(
        "{label}: L1 error first(400)={first:.4} second(200)={coarse:.4} second(400)={fine:.4}, rate {rate:.2}"
    );

    assert!(
        fine < 0.25 * first,
        "{label}: second order {fine} should be at least 4x better than first order {first}"
    );
    assert!(
        rate > 0.9,
        "{label}: observed convergence rate {rate:.2} is below 0.9"
    );
}

#[test]
fn wet_bed_dam_break_matches_stoker() {
    check_dam_break("wet bed (Stoker)", 2.0, 0.5);
}

#[test]
fn dry_bed_dam_break_matches_ritter() {
    check_dam_break("dry bed (Ritter)", 2.0, 0.0);
}
