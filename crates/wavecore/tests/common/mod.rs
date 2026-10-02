//! Helpers shared by the wave tests: harmonic analysis of a run.

use wavecore::{GHOST, Grid, Solver, State};

/// Complex amplitude of the component at angular frequency `omega`, per probe cell.
pub struct Harmonic {
    omega: f64,
    window: (f64, f64),
    sums: Vec<(f64, f64)>,
}

impl Harmonic {
    pub fn new(omega: f64, window: (f64, f64), probes: usize) -> Self {
        Self {
            omega,
            window,
            sums: vec![(0.0, 0.0); probes],
        }
    }

    pub fn record(&mut self, t: f64, dt: f64, eta: &[f64]) {
        if t <= self.window.0 || t > self.window.1 {
            return;
        }
        let (c, s) = ((self.omega * t).cos(), (self.omega * t).sin());
        for (acc, e) in self.sums.iter_mut().zip(eta) {
            acc.0 += e * c * dt;
            acc.1 -= e * s * dt;
        }
    }

    /// `(amplitude, phase)` per probe; the window must span whole periods.
    pub fn result(&self) -> Vec<(f64, f64)> {
        let scale = 2.0 / (self.window.1 - self.window.0);
        self.sums
            .iter()
            .map(|&(re, im)| (scale * re.hypot(im), im.atan2(re)))
            .collect()
    }
}

/// Run to `t_end`, analysing the free surface at `probes` (padded cells) over `window`.
pub fn run_and_analyse(
    solver: &Solver,
    s: &mut State,
    t_end: f64,
    omega: f64,
    window: (f64, f64),
    probes: &[(usize, usize)],
) -> Vec<(f64, f64)> {
    let g = &solver.grid;
    let mut harmonic = Harmonic::new(omega, window, probes.len());
    while s.time < t_end {
        let dt = solver.stable_dt(s);
        assert!(
            dt.is_finite(),
            "no water in the domain: nothing can propagate"
        );
        let dt = dt.min(t_end - s.time);
        solver.step(s, dt);
        let eta: Vec<f64> = probes
            .iter()
            .map(|&(i, j)| s.h[g.idx(i, j)] + solver.bed.b[g.idx(i, j)])
            .collect();
        harmonic.record(s.time, dt, &eta);
    }
    harmonic.result()
}

/// Padded index of the cell containing coordinate `x` metres.
pub fn col(grid: &Grid, x: f64) -> usize {
    GHOST + (x / grid.dx) as usize
}
