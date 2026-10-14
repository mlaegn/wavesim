//! Handing a wave from the main model to the slice solver before it breaks.
//!
//! The main model (Serre-Green-Naghdi, depth-averaged) carries a wave across deep water cheaply
//! and accurately, but cannot follow it once its face steepens past what a depth-averaged surface
//! can hold. The slice solver follows the full shape of the break, but costs grow with the cube
//! of its nodes, so it should cover only the last stretch. Here the main model runs the wave, in
//! a strip one line wide, until its crest reaches water `handoff_depth` deep; its surface and
//! depth-averaged velocity then become the slice solver's starting state over the water from a
//! little behind the wave to the shore.
//!
//! The slice solver needs the potential on the surface, whose derivative along it is the water's
//! velocity there, `u + w eta_x`. The main model has the depth-averaged velocity `ubar`; its own
//! profile in depth gives the surface's, `u = ubar - H^2/3 ubar_xx` along and
//! `w = -H ubar_x - ubar h_x` up, `H` the local water depth and `h` the still depth: the same
//! rule [`waveslice::Solitary`] starts from, with the bed's slope added to `w`.
//!
//! The breaking switch of the main model is off during the approach. It turns dispersion off
//! where a wave is tall for its depth, from 0.3 of the depth, which is right for a wave about
//! to break and wrong for one still shoaling, as a solitary wave 0.3 of the depth is from the
//! start; the hand-off comes before the break.

use wavecore::{Bathymetry, Breaking, Dispersion, GHOST, Grid, Order, Solver, State};

use crate::Error;

/// The main model's wave at the moment of hand-off, along the line: cell centres, the surface
/// above still water and the depth-averaged velocity.
#[derive(Clone, Debug)]
pub struct Approach {
    pub time: f64,
    pub x: Vec<f64>,
    pub eta: Vec<f64>,
    pub ubar: Vec<f64>,
    /// Where the crest was.
    pub crest: f64,
    pub steps: usize,
}

/// Run the main model, starting from `eta(x)` and `ubar(x)` over the bed `depth(x)` on cells
/// `dx` wide from 0 to `length`, until the crest is over water `handoff_depth` deep.
pub fn approach(
    length: f64,
    dx: f64,
    depth: &dyn Fn(f64) -> f64,
    eta: &dyn Fn(f64) -> f64,
    ubar: &dyn Fn(f64) -> f64,
    handoff_depth: f64,
) -> Result<Approach, Error> {
    let nx = (length / dx).round() as usize;
    let grid = Grid::new(nx, 4, length / nx as f64, length / nx as f64);
    let never = Breaking {
        eta_ratio: (1e9, 2e9),
        slope: (1e9, 2e9),
    };
    let solver = Solver::new(grid, Bathymetry::from_fn(&grid, |x, _| -depth(x)))
        .with_order(Order::Fifth)
        .with_dispersion(Dispersion {
            breaking: never,
            ..Dispersion::default()
        });
    let mut s = State::zeros(&grid);
    for (i, j) in grid.interior() {
        let (x, _) = grid.centre(i, j);
        let k = grid.idx(i, j);
        s.h[k] = depth(x) + eta(x);
        s.hu[k] = s.h[k] * ubar(x);
    }
    let row = GHOST + 1;
    let line = |s: &State| -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut out = (
            Vec::with_capacity(nx),
            Vec::with_capacity(nx),
            Vec::with_capacity(nx),
        );
        for i in GHOST..GHOST + nx {
            let (x, _) = grid.centre(i, row);
            let k = grid.idx(i, row);
            out.0.push(x);
            out.1.push(s.h[k] - depth(x));
            out.2
                .push(if s.h[k] > 0.0 { s.hu[k] / s.h[k] } else { 0.0 });
        }
        out
    };
    let mut steps = 0;
    loop {
        let (x, eta, ubar) = line(&s);
        let top = (0..nx)
            .max_by(|&a, &b| eta[a].total_cmp(&eta[b]))
            .unwrap_or(0);
        let crest = x[top];
        if depth(crest) <= handoff_depth {
            return Ok(Approach {
                time: s.time,
                x,
                eta,
                ubar,
                crest,
                steps,
            });
        }
        let dt = solver.stable_dt(&s);
        if !dt.is_finite() || steps > 1_000_000 {
            return Err(Error::Unstable {
                time: s.time,
                what: "the approach".into(),
            });
        }
        solver.step(&mut s, dt);
        steps += 1;
    }
}

impl Approach {
    /// The surface above still water at `x`, straight between cell centres.
    pub fn eta_at(&self, x: f64) -> f64 {
        interpolate(&self.x, &self.eta, x)
    }

    /// The potential on the surface at each of `xs` (increasing), relative to its value at the
    /// first: the integral along the surface of `u + w eta_x` from the main model's profile in
    /// depth, over the bed `depth(x)`.
    pub fn potential(&self, xs: &[f64], depth: &dyn Fn(f64) -> f64) -> Vec<f64> {
        let n = self.x.len();
        let dx = self.x[1] - self.x[0];
        let d = |f: &[f64], i: usize| (f[(i + 1).min(n - 1)] - f[i.saturating_sub(1)]) / (2.0 * dx);
        let d2 = |f: &[f64], i: usize| {
            let (a, b) = (i.saturating_sub(1), (i + 1).min(n - 1));
            (f[b] - 2.0 * f[i] + f[a]) / (dx * dx)
        };
        let slope: Vec<f64> = (0..n)
            .map(|i| {
                let x = self.x[i];
                let still = depth(x);
                let big = still + self.eta[i];
                let bed_slope = (depth(x + 0.5 * dx) - depth(x - 0.5 * dx)) / dx;
                let u = self.ubar[i] - big * big / 3.0 * d2(&self.ubar, i);
                let w = -big * d(&self.ubar, i) - self.ubar[i] * bed_slope;
                u + w * d(&self.eta, i)
            })
            .collect();
        // The integral at every cell centre, then at the nodes.
        let mut along = vec![0.0; n];
        for i in 1..n {
            along[i] = along[i - 1] + 0.5 * (slope[i - 1] + slope[i]) * dx;
        }
        let first = interpolate(&self.x, &along, xs[0]);
        xs.iter()
            .map(|&x| interpolate(&self.x, &along, x) - first)
            .collect()
    }
}

/// `values` at `x`, straight between the points `at` (increasing), held past the ends.
fn interpolate(at: &[f64], values: &[f64], x: f64) -> f64 {
    let k = at.partition_point(|&a| a <= x).clamp(1, at.len() - 1);
    let w = ((x - at[k - 1]) / (at[k] - at[k - 1])).clamp(0.0, 1.0);
    values[k - 1] + (values[k] - values[k - 1]) * w
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolation_is_straight_and_holds_past_the_ends() {
        let (at, v) = ([0.0, 1.0, 3.0], [0.0, 2.0, 0.0]);
        assert_eq!(interpolate(&at, &v, 0.5), 1.0);
        assert_eq!(interpolate(&at, &v, 2.0), 1.0);
        assert_eq!(interpolate(&at, &v, -1.0), 0.0);
        assert_eq!(interpolate(&at, &v, 9.0), 0.0);
    }

    #[test]
    fn the_handed_potential_of_a_solitary_wave_matches_its_own_start() {
        // Flat bed, no time stepped: the potential built from the main model's fields must be
        // the solitary wave's own, which uses the same profile in depth.
        let wave = waveslice::Solitary {
            g: 9.81,
            depth: 1.0,
            height: 0.3,
            crest: 10.0,
        };
        let c = wave.speed();
        let flat = |_: f64| 1.0;
        let a = approach(
            30.0,
            0.02,
            &flat,
            &|x| wave.eta(x),
            &|x| c * wave.eta(x) / (1.0 + wave.eta(x)),
            2.0,
        )
        .unwrap();
        assert_eq!(
            a.steps, 0,
            "the crest is already over water shallower than 2 m"
        );
        let xs: Vec<f64> = (0..=200).map(|k| 1.0 + 0.1 * k as f64).collect();
        let handed = a.potential(&xs, &flat);
        let own = wave.potential(&xs);
        let worst = handed
            .iter()
            .zip(&own)
            .map(|(h, o)| (h - (o - own[0])).abs())
            .fold(0.0, f64::max);
        let scale = own[own.len() - 1] - own[0];
        assert!(worst < 2e-3 * scale, "differs by {worst} of {scale}");
    }
}
