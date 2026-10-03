//! Run a monochromatic swell over a bed file and write frames.
//!
//! The swell travels along +x of the bed's frame (towards the beach), so the direction
//! is chosen when the bed is fetched. The wave-maker makes a sinusoid, so it goes in the
//! shallowest water where a sinusoid is still the right wave (see [`MAX_SECOND_HARMONIC`]);
//! a sponge behind it absorbs the seaward-going half, and the bed further offshore than the
//! sponge needs is cropped off. The sides are walls: for a swell travelling along +x a wall is
//! a mirror, so it leaves the incoming swell untouched.

use std::path::PathBuf;
use std::time::Instant;

use wavecore::{
    Dispersion, DispersionStats, G, GHOST, Order, Solver, Sponge, State, WaveMaker, linear_wave,
};
pub use waveio::Bed;
use waveio::RunWriter;

#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error(transparent)]
    Io(#[from] waveio::Error),
    #[error("{0}")]
    Setup(String),
    #[error("the solution blew up at t = {time:.1} s ({what})")]
    Unstable { time: f64, what: String },
}

#[derive(Clone, Debug)]
pub struct RunOptions {
    /// Path of the bed file's JSON header.
    pub bed: PathBuf,
    /// Output directory.
    pub out: PathBuf,
    /// Wave height (trough to crest) in metres at the wave-maker.
    pub height: f64,
    /// Wave period in seconds.
    pub period: f64,
    /// Still-water level in metres above mean sea level.
    pub tide: f64,
    /// Simulated seconds; `None` plans it (see [`Plan::duration`]).
    pub duration: Option<f64>,
    /// Seconds between output frames.
    pub frame_interval: f64,
    /// Manning roughness in s/m^(1/3); zero switches friction off.
    pub manning: f64,
    /// Use the dispersive (Serre-Green-Naghdi) equations with the third-order scheme.
    /// `false` runs plain shallow water with the MC limiter, in which tall waves steepen
    /// into shocks wherever they are.
    pub dispersive: bool,
    /// Stop after this many seconds of wall-clock time and keep what has been simulated, as a
    /// valid run marked `truncated`. `None` runs to the end. This is the guard that stops a
    /// run from using a laptop for longer than intended.
    pub max_wall_seconds: Option<f64>,
}

/// Where the forcing goes and how long to run, worked out from the bed and the options.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plan {
    /// Columns at the offshore end of the bed that the run leaves out: the wave-maker and its
    /// sponge do not need them.
    pub crop: usize,
    /// Position of the wave-maker line in metres along +x of the bed file (not of the cropped
    /// run).
    pub maker_x: f64,
    /// Still-water depth at the maker, along the middle row.
    pub maker_depth: f64,
    /// The shallowest still-water depth along the wave-maker line.
    pub maker_depth_min: f64,
    /// The second harmonic a sinusoid leaves out there, relative to the wave, at the shallowest
    /// point of the line; at most [`MAX_SECOND_HARMONIC`].
    pub second_harmonic: f64,
    /// Wavelength at the maker in metres.
    pub wavelength: f64,
    /// Width in cells of the offshore sponge.
    pub sponge_offshore: usize,
    /// Simulated seconds: the options' duration, or else long enough for the swell to ramp up,
    /// reach the shore along the middle row at its group speed, and break [`BREAKS_SHOWN`] times.
    pub duration: f64,
}

/// The largest second harmonic, relative to the wave, that the sinusoidal wave-maker may leave
/// out where it stands.
///
/// A wave of finite height has sharper crests and flatter troughs than a sinusoid: to second
/// order (Stokes) a second harmonic is bound to it, and it grows quickly as the water gets
/// shallow (the Ursell number). Where it is large, as on a shelf a swell has to cross before
/// it reaches the reef, the sinusoid is the wrong wave: it sheds the harmonic it lacks as a
/// free wave, which beats with the swell and moves where the waves break. So the wave-maker
/// goes no shallower than where that harmonic is a tenth of the wave, and the shoaling that
/// steepens the swell from there on is the model's own.
pub const MAX_SECOND_HARMONIC: f64 = 0.10;

/// How many waves a planned run shows breaking after the first full-height one reaches the shore.
pub const BREAKS_SHOWN: f64 = 4.0;

/// Amplitude of the second harmonic bound to a Stokes wave of amplitude `a` and wavenumber `k`
/// in still-water depth `h`, relative to `a`: `(ka/4) cosh(kh) (2 + cosh 2kh) / sinh^3(kh)`.
/// It tends to `ka/2` in deep water and to `3a / (4 k^2 h^3)` in shallow water.
pub fn second_harmonic(a: f64, k: f64, h: f64) -> f64 {
    let kh = k * h;
    0.25 * k * a * kh.cosh() * (2.0 + (2.0 * kh).cosh()) / kh.sinh().powi(3)
}

pub struct Progress {
    pub time: f64,
    pub duration: f64,
    pub steps: u64,
    pub frames: usize,
    /// Highest water surface above still water among wet cells, in metres.
    pub max_rise: f64,
}

#[derive(Debug)]
pub struct Summary {
    pub steps: u64,
    pub frames: usize,
    pub wall_seconds: f64,
    pub max_rise: f64,
    pub header: PathBuf,
    pub plan: Plan,
    /// How hard the dispersive solves worked, or `None` for shallow water.
    pub dispersion: Option<DispersionStats>,
    /// The run hit its wall-clock limit before reaching `duration`; the frames written so
    /// far are a valid, shorter run.
    pub truncated: bool,
}

const SPONGE_STRENGTH: f64 = 1.5;
const SIGMA_CELLS: f64 = 2.0;
const RAMP_PERIODS: f64 = 2.0;

fn setup(msg: impl Into<String>) -> Error {
    Error::Setup(msg.into())
}

fn validate(opts: &RunOptions) -> Result<(), Error> {
    let positive = [
        ("height", opts.height),
        ("period", opts.period),
        ("duration", opts.duration.unwrap_or(1.0)),
        ("frame interval", opts.frame_interval),
    ];
    for (name, v) in positive {
        if !(v.is_finite() && v > 0.0) {
            return Err(setup(format!("{name} must be a positive number, got {v}")));
        }
    }
    if !(opts.manning.is_finite() && opts.manning >= 0.0 && opts.tide.is_finite()) {
        return Err(setup("manning must be >= 0 and tide must be a number"));
    }
    Ok(())
}

/// Decide where the wave-maker and the sponge go, and how long to run.
///
/// Coming from offshore, the wave-maker goes in the last column before the second harmonic a
/// sinusoid leaves out exceeds [`MAX_SECOND_HARMONIC`] anywhere along the line: as close to the
/// reef as a sinusoid is still the right wave.
pub fn plan(bed: &Bed, opts: &RunOptions) -> Result<Plan, Error> {
    validate(opts)?;
    let grid = bed.grid();
    let shifted = bed.bathymetry(opts.tide);
    let depth = |i: usize, j: usize| -shifted.b[grid.idx(GHOST + i, GHOST + j)];
    let shallowest = |i: usize| {
        (0..grid.ny)
            .map(|j| depth(i, j))
            .fold(f64::INFINITY, f64::min)
    };
    let omega = std::f64::consts::TAU / opts.period;
    let a = opts.height / 2.0;
    // The dispersive model has no wave of this frequency once omega^2 h / g reaches 3.
    let too_short = |h: f64| opts.dispersive && omega * omega * h / G >= 3.0;
    let harmonic = |h: f64| second_harmonic(a, linear_wave(omega, h, opts.dispersive).k, h);
    let sinusoid_fits = |h: f64| h > 0.0 && !too_short(h) && harmonic(h) <= MAX_SECOND_HARMONIC;

    let offshore = shallowest(0);
    if too_short(offshore) {
        return Err(setup(format!(
            "a {} s wave is too short for the dispersive model in the {offshore:.1} m of water \
             offshore; use a longer --period, or a bed clipped shallower offshore",
            opts.period
        )));
    }
    if !sinusoid_fits(offshore) {
        return Err(setup(format!(
            "the bed is only {offshore:.1} m deep at its offshore edge at this tide, where a {} m, \
             {} s swell is far from a sinusoid (its second harmonic is {:.0}% of it, the limit is \
             {:.0}%); fetch a bed that reaches deeper water (offshore_m in spots.toml)",
            opts.height,
            opts.period,
            100.0 * harmonic(offshore.max(1e-3)),
            100.0 * MAX_SECOND_HARMONIC
        )));
    }
    let i = (1..grid.nx)
        .find(|&i| !sinusoid_fits(shallowest(i)))
        .map_or(grid.nx - 1, |i| i - 1);
    let maker_x = (i as f64 + 0.5) * grid.dx;
    let h = depth(i, grid.ny / 2);
    let wavelength = std::f64::consts::TAU / linear_wave(omega, h, opts.dispersive).k;

    let sponge_offshore =
        ((1.5 * wavelength / grid.dx).ceil() as usize).clamp(4, (grid.nx / 3).max(4));
    let needed = sponge_offshore as f64 * grid.dx + 4.0 * SIGMA_CELLS * grid.dx;
    if maker_x < needed {
        return Err(setup(format!(
            "no room for the offshore sponge: the wave-maker would sit at x = {maker_x:.0} m but \
             needs {needed:.0} m behind it; fetch a bed that reaches further offshore"
        )));
    }

    // Time for the swell to reach the shore along the middle row, at its linear group speed
    // (in water at least half a metre deep, so the last cells before the shore do not count
    // for more than they should).
    let travel: f64 = (i..grid.nx)
        .map(|n| depth(n, grid.ny / 2))
        .take_while(|&d| d > 0.0)
        .map(|d| grid.dx / linear_wave(omega, d.max(0.5), opts.dispersive).group_speed)
        .sum();
    let duration = opts
        .duration
        .unwrap_or((RAMP_PERIODS + BREAKS_SHOWN) * opts.period + travel);

    Ok(Plan {
        crop: ((maker_x - needed) / grid.dx).floor() as usize,
        maker_x,
        maker_depth: h,
        maker_depth_min: shallowest(i),
        second_harmonic: harmonic(shallowest(i)),
        wavelength,
        sponge_offshore,
        duration,
    })
}

/// A rough guess of how long a run will take in wall-clock seconds on `threads` threads, so
/// that a long one can be seen coming before it starts. It is calibrated on one Apple M5 (the
/// Pipeline bed, cropped to 616 x 200 cells of 3 m, took 204 ns per cell per step on 4 threads
/// with dispersion), so treat it as the right order of magnitude, not a promise.
pub fn estimate_seconds(bed: &Bed, opts: &RunOptions, plan: &Plan, threads: usize) -> f64 {
    let bed = bed.crop_x(plan.crop);
    let grid = bed.grid();
    let deepest = bed
        .elevation
        .iter()
        .fold(f64::INFINITY, |m, &z| m.min(f64::from(z)));
    let depth = (opts.tide - deepest).max(1.0);
    let c = (G * depth).sqrt();
    let dt = 0.4 / (c / grid.dx + c / grid.dy);
    let steps = plan.duration / dt;
    let cells = (grid.nx * grid.ny) as f64;
    let one_thread = if opts.dispersive { 620e-9 } else { 310e-9 };
    let speedup = (threads.clamp(1, 6) as f64).powf(0.8);
    cells * steps * one_thread / speedup
}

/// Run the simulation and write the run directory. `on_progress` is called after
/// every frame.
pub fn run(opts: &RunOptions, mut on_progress: impl FnMut(&Progress)) -> Result<Summary, Error> {
    let started = Instant::now();
    let full = Bed::read(&opts.bed)?;
    let plan = plan(&full, opts)?;
    let bed = full.crop_x(plan.crop);
    let grid = bed.grid();
    let maker_x = plan.maker_x - plan.crop as f64 * grid.dx;

    let maker = WaveMaker {
        x: maker_x,
        amplitude: opts.height / 2.0,
        period: opts.period,
        angle: 0.0,
        sigma: SIGMA_CELLS * grid.dx,
        ramp_periods: RAMP_PERIODS,
    };
    let mut solver = Solver::new(grid, bed.bathymetry(opts.tide));
    if opts.dispersive {
        solver = solver
            .with_order(Order::Third)
            .with_dispersion(Dispersion::default());
    }
    let solver = solver
        .with_wavemaker(maker)
        .with_sponge(
            Sponge::new(plan.sponge_offshore, SPONGE_STRENGTH).edges(true, false, false, false),
        )
        .with_manning(opts.manning);
    let mut state = State::lake_at_rest(&grid, &solver.bed, 0.0);

    let settings = serde_json::json!({
        "wave_height_m": opts.height,
        "wave_period_s": opts.period,
        "tide_m": opts.tide,
        "duration_s": plan.duration,
        "frame_interval_s": opts.frame_interval,
        "maker_x_m": maker_x,
        "cropped_offshore_m": plan.crop as f64 * grid.dx,
        "maker_sigma_m": SIGMA_CELLS * grid.dx,
        // Where the wave-maker's own bump has died away. Seaward of this the surface is the
        // source and the sponge, not sea, so a viewer should start here.
        "near_field_end_m": maker_x + 6.0 * SIGMA_CELLS * grid.dx,
        "maker_depth_m": plan.maker_depth,
        "maker_depth_min_m": plan.maker_depth_min,
        "maker_second_harmonic": plan.second_harmonic,
        "sponge_offshore_cells": plan.sponge_offshore,
        "manning": opts.manning,
        "direction": "along +x of the bed frame",
        "model": if opts.dispersive {
            "Serre-Green-Naghdi (flat-bed operator), hybrid breaking, third-order upwind-biased scheme"
        } else {
            "nonlinear shallow water, MC limiter, non-dispersive"
        },
        "dispersive": opts.dispersive,
    });
    let mut writer = RunWriter::create(&opts.out, &bed, &["eta", "breaking"], settings)?;

    let mut eta = vec![0.0_f32; grid.nx * grid.ny];
    let mut breaking = vec![0.0_f32; grid.nx * grid.ny];
    let mut frames = 0;
    let mut steps = 0_u64;
    let mut max_rise = 0.0_f64;
    let mut next_frame = 0.0_f64;
    let mut truncated = false;

    let snapshot = |state: &State, eta: &mut [f32]| -> f64 {
        let mut rise = 0.0_f64;
        for (n, (i, j)) in grid.interior().enumerate() {
            let k = grid.idx(i, j);
            eta[n] = (state.h[k] + f64::from(bed.elevation[n])) as f32;
            if state.h[k] > 1e-3 {
                rise = rise.max(f64::from(eta[n]) - opts.tide);
            }
        }
        rise
    };

    loop {
        if state.time >= next_frame - 1e-9 {
            let rise = snapshot(&state, &mut eta);
            if !rise.is_finite() || rise > 100.0 || eta.iter().any(|v| !v.is_finite()) {
                return Err(Error::Unstable {
                    time: state.time,
                    what: format!("surface {rise:.1} m above still water"),
                });
            }
            max_rise = max_rise.max(rise);
            for (b, v) in breaking.iter_mut().zip(solver.breaking_indicator(&state)) {
                *b = v as f32;
            }
            writer.write_frame(state.time, &[&eta, &breaking])?;
            frames += 1;
            next_frame += opts.frame_interval;
            on_progress(&Progress {
                time: state.time,
                duration: plan.duration,
                steps,
                frames,
                max_rise,
            });
        }
        if state.time >= plan.duration - 1e-9 {
            break;
        }
        if opts
            .max_wall_seconds
            .is_some_and(|limit| started.elapsed().as_secs_f64() > limit)
        {
            truncated = true;
            break;
        }
        let stable = solver.stable_dt(&state);
        if !stable.is_finite() {
            return Err(setup("no water in the domain at this tide"));
        }
        let dt = stable
            .min(next_frame - state.time)
            .min(plan.duration - state.time);
        solver.step(&mut state, dt);
        steps += 1;
    }

    if truncated {
        writer.note("truncated", serde_json::json!(true));
        writer.note("stopped_at_s", serde_json::json!(state.time));
    }
    let header = writer.finish()?;
    Ok(Summary {
        steps,
        frames,
        wall_seconds: started.elapsed().as_secs_f64(),
        max_rise,
        header,
        plan,
        dispersion: solver.dispersion_stats(),
        truncated,
    })
}
