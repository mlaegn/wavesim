//! Run a monochromatic swell over a bed file and write frames.
//!
//! The swell travels along +x of the bed's frame (towards the beach), so the direction
//! is chosen when the bed is fetched. A wave-maker line is placed where the still-water
//! depth first falls to `maker_depth` going shoreward; a sponge behind it absorbs the
//! seaward-going half, and sponges on the sides let the swell leave sideways.

use std::path::PathBuf;
use std::time::Instant;

use wavecore::{
    Dispersion, DispersionStats, G, GHOST, Order, Solver, Sponge, State, WaveMaker, linear_wave,
};
use waveio::{Bed, RunWriter};

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
    /// Simulated seconds.
    pub duration: f64,
    /// Seconds between output frames.
    pub frame_interval: f64,
    /// The wave-maker goes where the still-water depth first falls to this, in metres.
    pub maker_depth: f64,
    /// Manning roughness in s/m^(1/3); zero switches friction off.
    pub manning: f64,
    /// Use the dispersive (Serre-Green-Naghdi) equations with the third-order scheme.
    /// `false` runs plain shallow water with the MC limiter, in which tall waves steepen
    /// into shocks wherever they are.
    pub dispersive: bool,
}

/// Where the forcing goes, worked out from the bed and the options.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Plan {
    /// Position of the wave-maker line in metres along +x.
    pub maker_x: f64,
    /// Still-water depth at the maker, along the middle row.
    pub maker_depth: f64,
    /// Long-wave wavelength at the maker in metres.
    pub wavelength: f64,
    /// Width in cells of the offshore sponge.
    pub sponge_offshore: usize,
    /// Width in cells of the sponges on the two sides.
    pub sponge_side: usize,
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
}

const SPONGE_STRENGTH: f64 = 1.5;
const SIGMA_CELLS: f64 = 2.0;

fn setup(msg: impl Into<String>) -> Error {
    Error::Setup(msg.into())
}

fn validate(opts: &RunOptions) -> Result<(), Error> {
    let positive = [
        ("height", opts.height),
        ("period", opts.period),
        ("duration", opts.duration),
        ("frame interval", opts.frame_interval),
        ("maker depth", opts.maker_depth),
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

/// Decide where the wave-maker and the sponges go.
pub fn plan(bed: &Bed, opts: &RunOptions) -> Result<Plan, Error> {
    validate(opts)?;
    let grid = bed.grid();
    let shifted = bed.bathymetry(opts.tide);
    let j = GHOST + grid.ny / 2;
    let depth = |i: usize| -shifted.b[grid.idx(GHOST + i, j)];

    if depth(0) <= opts.maker_depth {
        return Err(setup(format!(
            "the offshore edge is only {:.1} m deep at this tide, not deeper than the requested \
             maker depth of {} m; use a smaller --maker-depth or a bed that extends further offshore",
            depth(0),
            opts.maker_depth
        )));
    }
    let Some(i) = (1..grid.nx).find(|&i| depth(i) <= opts.maker_depth) else {
        return Err(setup(format!(
            "the water never gets shallower than {} m along the middle row; use a larger --maker-depth",
            opts.maker_depth
        )));
    };
    let maker_x = (i as f64 + 0.5) * grid.dx;
    let h = depth(i);
    let omega = std::f64::consts::TAU / opts.period;
    if opts.dispersive && omega * omega * h / G >= 3.0 {
        return Err(setup(format!(
            "a {} s wave is too short for the dispersive model in {h:.1} m of water; \
             use a longer --period or a shallower --maker-depth",
            opts.period
        )));
    }
    if opts.height > 0.6 * h {
        return Err(setup(format!(
            "a {} m wave is too big for {h:.2} m of water at the wave-maker (it would break at \
             once); use a deeper --maker-depth or a smaller --height",
            opts.height
        )));
    }
    let wavelength = std::f64::consts::TAU / linear_wave(omega, h, opts.dispersive).k;

    let cells = |spacing: f64| (1.5 * wavelength / spacing).ceil() as usize;
    let sponge_offshore = cells(grid.dx).clamp(4, (grid.nx / 3).max(4));
    let sponge_side = cells(grid.dy).clamp(4, (grid.ny / 4).max(4));

    let needed = sponge_offshore as f64 * grid.dx + 4.0 * SIGMA_CELLS * grid.dx;
    if maker_x < needed {
        return Err(setup(format!(
            "no room for the offshore sponge: the wave-maker would sit at x = {maker_x:.0} m but \
             needs {needed:.0} m behind it; use a deeper --maker-depth or a longer offshore bed"
        )));
    }
    Ok(Plan {
        maker_x,
        maker_depth: h,
        wavelength,
        sponge_offshore,
        sponge_side,
    })
}

/// Run the simulation and write the run directory. `on_progress` is called after
/// every frame.
pub fn run(opts: &RunOptions, mut on_progress: impl FnMut(&Progress)) -> Result<Summary, Error> {
    let started = Instant::now();
    let bed = Bed::read(&opts.bed)?;
    let plan = plan(&bed, opts)?;
    let grid = bed.grid();

    let maker = WaveMaker {
        x: plan.maker_x,
        amplitude: opts.height / 2.0,
        period: opts.period,
        angle: 0.0,
        sigma: SIGMA_CELLS * grid.dx,
        ramp_periods: 2.0,
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
        .with_sponge(Sponge::new(plan.sponge_side, SPONGE_STRENGTH).edges(false, false, true, true))
        .with_manning(opts.manning);
    let mut state = State::lake_at_rest(&grid, &solver.bed, 0.0);

    let settings = serde_json::json!({
        "wave_height_m": opts.height,
        "wave_period_s": opts.period,
        "tide_m": opts.tide,
        "duration_s": opts.duration,
        "frame_interval_s": opts.frame_interval,
        "maker_x_m": plan.maker_x,
        "maker_depth_m": plan.maker_depth,
        "sponge_offshore_cells": plan.sponge_offshore,
        "sponge_side_cells": plan.sponge_side,
        "manning": opts.manning,
        "direction": "along +x of the bed frame",
        "model": if opts.dispersive {
            "Serre-Green-Naghdi (flat-bed operator), hybrid breaking, third-order upwind-biased scheme"
        } else {
            "nonlinear shallow water, MC limiter, non-dispersive"
        },
        "dispersive": opts.dispersive,
    });
    let mut writer = RunWriter::create(&opts.out, &bed, settings)?;

    let mut eta = vec![0.0_f32; grid.nx * grid.ny];
    let mut frames = 0;
    let mut steps = 0_u64;
    let mut max_rise = 0.0_f64;
    let mut next_frame = 0.0_f64;

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
            writer.write_frame(state.time, &eta)?;
            frames += 1;
            next_frame += opts.frame_interval;
            on_progress(&Progress {
                time: state.time,
                duration: opts.duration,
                steps,
                frames,
                max_rise,
            });
        }
        if state.time >= opts.duration - 1e-9 {
            break;
        }
        let stable = solver.stable_dt(&state);
        if !stable.is_finite() {
            return Err(setup("no water in the domain at this tide"));
        }
        let dt = stable
            .min(next_frame - state.time)
            .min(opts.duration - state.time);
        solver.step(&mut state, dt);
        steps += 1;
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
    })
}
