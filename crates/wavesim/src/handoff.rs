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

use std::path::PathBuf;

use wavecore::{
    Bathymetry, Breaking, Dispersion, GHOST, Grid, Order, Solver, Sponge, State, WaveMaker,
    linear_wave,
};
use waveio::{BED_FORMAT, Bed, BedHeader, Frame};

use crate::{
    Error, Layout, RAMP_PERIODS, RunOptions, SETTLE_PERIODS, SIGMA_CELLS, SPONGE_STRENGTH, beds,
    plan, solver_for,
};

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

/// A swell `height` high (crest to trough at the wave-maker) and `period` long, run by the main
/// model along one line of real bed, `(x, z)` with `z` up from still water, on cells `cell`
/// wide, exactly as `wavesim run` runs it: the same wave-maker placement, sponge, ramp and
/// breaking switch. Once the waves have reached the shore and the surf zone has had
/// [`SETTLE_PERIODS`] to settle, the next crest to pass the last point where the water is
/// `handoff_depth` deep is handed over, with the main model's whole line at that moment.
///
/// The breaking switch stays on here, unlike in [`approach`]: the waves ahead break, and the
/// surf zone they make is part of the sea the followed wave runs into. At the hand-off the
/// followed crest should still be well under 0.3 of the depth, where the switch starts; the
/// returned [`Approach::eta_ratio`] says how far under.
pub fn swell(
    line: &[(f64, f64)],
    cell: f64,
    height: f64,
    period: f64,
    handoff_depth: f64,
) -> Result<Approach, Error> {
    let (x0, x1) = (line[0].0, line[line.len() - 1].0);
    let nx = ((x1 - x0) / cell).floor() as usize;
    let ny = 4;
    let z_at = |x: f64| {
        interpolate(
            &line.iter().map(|p| p.0).collect::<Vec<_>>(),
            &line.iter().map(|p| p.1).collect::<Vec<_>>(),
            x,
        )
    };
    let row: Vec<f32> = (0..nx)
        .map(|i| z_at(x0 + (i as f64 + 0.5) * cell) as f32)
        .collect();
    let bed = Bed {
        header: BedHeader {
            format: BED_FORMAT.into(),
            version: 1,
            name: "line".into(),
            description: String::new(),
            nx,
            ny,
            dx: cell,
            dy: cell,
            data: String::new(),
            dtype: "f32-le".into(),
            layout: "row-major, y outer (increasing), x inner (increasing)".into(),
            units: "metres above mean sea level; negative is underwater".into(),
            frame: Frame {
                crs: String::new(),
                origin_easting: 0.0,
                origin_northing: 0.0,
                x_bearing_deg: 0.0,
                note: String::new(),
            },
            stats: serde_json::Value::Null,
            source: serde_json::Value::Null,
        },
        elevation: (0..ny).flat_map(|_| row.iter().copied()).collect(),
    };
    let opts = RunOptions {
        bed: PathBuf::new(),
        out: PathBuf::new(),
        height,
        period,
        tide: 0.0,
        duration: None,
        frame_interval: 1.0,
        manning: 0.0,
        dispersive: true,
        nested: false,
        threads: 1,
        background: false,
        budget_seconds: None,
        max_wall_seconds: None,
        handoff: Vec::new(),
    };
    let p = plan(&bed, &opts)?;
    let (cropped, _) = beds(
        &bed,
        &Layout {
            coarse: p,
            fine: None,
        },
    );
    let grid = cropped.grid();
    let offset = x0 + p.crop as f64 * cell;
    let maker = WaveMaker {
        x: p.maker_x - p.crop as f64 * cell,
        amplitude: height / 2.0,
        period,
        angle: 0.0,
        sigma: SIGMA_CELLS * cell,
        ramp_periods: RAMP_PERIODS,
    };
    let solver = solver_for(&cropped, &opts)
        .with_wavemaker(maker)
        .with_sponge(
            Sponge::new(p.sponge_offshore, SPONGE_STRENGTH).edges(true, false, false, false),
        );
    let mut s = State::lake_at_rest(&grid, &solver.bed, 0.0);
    let j = GHOST + ny / 2;
    let depth_of = |i: usize| -solver.bed.b[grid.idx(GHOST + i, j)];

    // The hand-off point: the last cell, coming from offshore, after which the water stays
    // shallower than `handoff_depth` all the way to the shore.
    let wet = (0..grid.nx).take_while(|&i| depth_of(i) > 0.0).count();
    let mut at = wet.saturating_sub(1);
    while at > 0 && depth_of(at - 1) < handoff_depth {
        at -= 1;
    }
    if at == 0 || depth_of(at) >= handoff_depth {
        return Err(Error::Setup(format!(
            "the water along the line never stays shallower than {handoff_depth} m before the shore"
        )));
    }
    let omega = std::f64::consts::TAU / period;
    let maker_cell = (p.maker_x / cell) as usize - p.crop;
    let travel: f64 = (maker_cell..=at)
        .map(|i| cell / linear_wave(omega, depth_of(i).max(0.5), true).group_speed)
        .sum();
    let ready = (RAMP_PERIODS + SETTLE_PERIODS) * period + travel;

    let k = |i: usize| grid.idx(GHOST + i, j);
    let mut history = [0.0_f64; 3];
    let mut steps = 0;
    loop {
        let dt = solver.stable_dt(&s);
        if !dt.is_finite() || s.time > ready + 4.0 * period {
            return Err(Error::Unstable {
                time: s.time,
                what: "the swell approach found no crest to hand over".into(),
            });
        }
        solver.step(&mut s, dt);
        steps += 1;
        history = [history[1], history[2], s.h[k(at)] - depth_of(at)];
        // The crest has just passed: the surface there rose and now falls.
        if s.time > ready
            && history[1] > history[0]
            && history[2] < history[1]
            && history[1] > 0.25 * height
        {
            break;
        }
    }
    let x: Vec<f64> = (0..grid.nx)
        .map(|i| offset + (i as f64 + 0.5) * cell)
        .collect();
    let eta: Vec<f64> = (0..grid.nx)
        .map(|i| (s.h[k(i)] - depth_of(i)).max(-depth_of(i)))
        .collect();
    let ubar: Vec<f64> = (0..grid.nx)
        .map(|i| {
            if s.h[k(i)] > 1e-3 {
                s.hu[k(i)] / s.h[k(i)]
            } else {
                0.0
            }
        })
        .collect();
    Ok(Approach {
        time: s.time,
        crest: x[at],
        x,
        eta,
        ubar,
        steps,
    })
}

/// A `wavesim-handoff` file, as `wavesim run --handoff-y` saves it.
#[derive(serde::Deserialize)]
struct HandoffFile {
    format: String,
    depth: f64,
    time: f64,
    crest: f64,
    x: Vec<f64>,
    eta: Vec<f64>,
    ubar: Vec<f64>,
    /// Where the settled sea broke seaward of the hand-off point, if it did (absent in files
    /// written before runs recorded it).
    #[serde(default)]
    broke_seaward_at: Option<f64>,
}

impl Approach {
    /// The wave a 2D run saved for the slice solver (see `wavesim run --handoff-y`), the depth
    /// it was saved at, and where the run's waves broke seaward of that point, if they did.
    pub fn read(path: &std::path::Path) -> Result<(Approach, f64, Option<f64>), Error> {
        let text = std::fs::read_to_string(path)
            .map_err(|e| Error::Setup(format!("{}: {e}", path.display())))?;
        let f: HandoffFile = serde_json::from_str(&text)
            .map_err(|e| Error::Setup(format!("{}: not a wavesim-handoff: {e}", path.display())))?;
        if f.format != "wavesim-handoff" || f.x.len() != f.eta.len() || f.x.len() != f.ubar.len() {
            return Err(Error::Setup(format!(
                "{}: not a wavesim-handoff",
                path.display()
            )));
        }
        let approach = Approach {
            time: f.time,
            x: f.x,
            eta: f.eta,
            ubar: f.ubar,
            crest: f.crest,
            steps: 0,
        };
        Ok((approach, f.depth, f.broke_seaward_at))
    }

    /// The followed crest's height over the still depth beneath it, the measure the main
    /// model's breaking switch starts at, at 0.3.
    pub fn eta_ratio(&self, depth: &dyn Fn(f64) -> f64) -> f64 {
        self.eta_at(self.crest) / depth(self.crest)
    }

    /// Only the followed wave: the surface and velocity between the troughs either side of the
    /// crest, tapered to still water over `taper` metres beyond each, still water further out.
    /// In a swell the wave ahead is breaking where this one will, and its foam cannot be handed
    /// to the slice solver; the crest's own shape and flow, which make its break, are kept.
    /// Returns the wave and its two troughs.
    pub fn single_wave(&self, wavelength: f64, taper: f64) -> (Approach, f64, f64) {
        let near = |lo: f64, hi: f64| -> usize {
            (0..self.x.len())
                .filter(|&i| self.x[i] >= lo && self.x[i] <= hi)
                .min_by(|&a, &b| self.eta[a].total_cmp(&self.eta[b]))
                .unwrap_or(0)
        };
        let crest = self.crest;
        let behind = self.x[near(crest - 0.75 * wavelength, crest)];
        let ahead = self.x[near(crest, crest + 0.75 * wavelength)];
        let weight = |x: f64| {
            let outside = if x < behind {
                behind - x
            } else if x > ahead {
                x - ahead
            } else {
                0.0
            };
            if outside >= taper {
                0.0
            } else {
                0.5 * (1.0 + (std::f64::consts::PI * outside / taper).cos())
            }
        };
        let eta = self
            .x
            .iter()
            .zip(&self.eta)
            .map(|(&x, &e)| e * weight(x))
            .collect();
        let ubar = self
            .x
            .iter()
            .zip(&self.ubar)
            .map(|(&x, &u)| u * weight(x))
            .collect();
        let wave = Approach {
            time: self.time,
            x: self.x.clone(),
            eta,
            ubar,
            crest: self.crest,
            steps: self.steps,
        };
        (wave, behind, ahead)
    }

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
