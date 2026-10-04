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
    Dispersion, DispersionStats, G, GHOST, Order, Relaxation, Solver, Sponge, State, WaveMaker,
    linear_wave,
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
    #[error(
        "this run would take about {} on {threads} thread(s){}, more than its budget of {}; \
         allow it with a larger --max-minutes, or make it smaller (a coarser bed with \
         tools/fetch_spot.py --cell 6, or a shorter --duration)",
        human(*.estimate), if *.background { " on the efficiency cores" } else { "" }, human(*.budget)
    )]
    OverBudget {
        estimate: f64,
        budget: f64,
        threads: usize,
        background: bool,
    },
}

/// "45 seconds", "3.5 minutes".
pub fn human(seconds: f64) -> String {
    if seconds < 90.0 {
        format!("{seconds:.0} seconds")
    } else {
        format!("{:.1} minutes", seconds / 60.0)
    }
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
    /// Use the dispersive (Serre-Green-Naghdi) equations with the fifth-order scheme.
    /// `false` runs plain shallow water with the MC limiter, in which tall waves steepen
    /// into shocks wherever they are.
    pub dispersive: bool,
    /// Run the swell's crossing on cells twice as large and only the reef on the bed's own
    /// cells, driven by the coarse run (see [`Fine`]); `false` runs everything on one grid.
    pub nested: bool,
    /// Threads the solver runs on; only used to estimate how long the run takes.
    pub threads: usize,
    /// The run is on macOS's efficiency cores (background priority), which are cool but about
    /// 3.5 times slower; only used to estimate how long the run takes.
    pub background: bool,
    /// Refuse to start a run estimated to take longer than this many seconds of wall-clock time
    /// ([`estimate_seconds`]). `None` starts anything.
    pub budget_seconds: Option<f64>,
    /// Stop after this many seconds of wall-clock time and keep what has been simulated, as a
    /// valid run marked `truncated`. `None` runs to the end. With [`RunOptions::budget_seconds`]
    /// it is the guard that keeps a run from using a laptop for longer than intended: the
    /// budget refuses what is seen coming, this stops what was not.
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

/// Where the fine grid of a nested run goes, worked out from the bed and the options.
///
/// The coarse run carries the swell from its wave-maker across the shelf on cells twice as
/// large, where the fifth-order scheme keeps its height; the fine run covers the reef and the
/// shore, across the bed's whole width, on the bed's own cells. Along its offshore edge a
/// relaxation zone pulls it towards the coarse run, which brings the swell in and lets what the
/// reef sends back out. Its sides are walls in the same places as the coarse run's, so the two
/// mirror the bed alike; relaxing the sides towards the coarse run instead imposed the coarse
/// surf zone, which it cannot resolve, on the fine one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Fine {
    /// First column of the bed file that the fine grid covers; it runs to the last.
    pub x0: usize,
    /// The still-water depth that sets where the reef begins: the deeper of
    /// [`FINE_DEPTH_PER_HEIGHT`] wave heights, and the depth where the swell is
    /// [`CELLS_PER_WAVELENGTH`] coarse cells long. The relaxation zone lies seaward of it.
    pub edge_depth: f64,
    /// Width in fine cells of the relaxation zone along the offshore edge: a wavelength at
    /// `edge_depth`.
    pub zone_offshore: usize,
    /// When the fine run starts, in simulated seconds: two periods before the swell's leading
    /// edge reaches it at the group speed. Until then its water is still.
    pub start: f64,
}

/// The grids of a run: the one with the wave-maker, and for a nested run the fine one.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Layout {
    /// The grid with the wave-maker: the bed's own cells, or twice as large for a nested run.
    pub coarse: Plan,
    pub fine: Option<Fine>,
}

/// A nested run's fine grid starts no shallower than this many wave heights, so that the coarse
/// run never breaks inside the zone that drives the fine one. The breaking switch starts to act
/// where the crest is 0.43 of the still-water depth; a shoaled crest is about 0.7 wave heights
/// up, so that is about 1.6 wave heights of water, and this leaves half as much again. The run
/// checks it: it warns if the swell broke inside the zone after all.
pub const FINE_DEPTH_PER_HEIGHT: f64 = 2.5;

/// ...and no shallower than where the swell is this many coarse cells long, which the
/// fifth-order scheme carries a kilometre without losing height.
pub const CELLS_PER_WAVELENGTH: f64 = 20.0;

#[derive(Debug)]
pub struct Summary {
    /// Steps of both grids together.
    pub steps: u64,
    /// Steps and wall-clock seconds of the coarse grid, and of the fine one in a nested run.
    pub coarse_steps: u64,
    pub coarse_seconds: f64,
    pub fine_steps: u64,
    pub fine_seconds: f64,
    pub frames: usize,
    pub wall_seconds: f64,
    pub max_rise: f64,
    /// The run's header; for a nested run, the fine run's.
    pub header: PathBuf,
    /// The coarse run's header in a nested run.
    pub coarse_header: Option<PathBuf>,
    pub layout: Layout,
    /// How hard the dispersive solves worked, or `None` for shallow water.
    pub dispersion: Option<DispersionStats>,
    /// The run hit its wall-clock limit before reaching `duration`; the frames written so
    /// far are a valid, shorter run.
    pub truncated: bool,
    /// Things the run noticed that make its answer less trustworthy.
    pub warnings: Vec<String>,
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

/// Lay out the grids of a run: [`plan`] alone, or for a nested run the plan on cells twice
/// as large and the [`Fine`] grid.
pub fn layout(bed: &Bed, opts: &RunOptions) -> Result<Layout, Error> {
    if !opts.nested {
        return Ok(Layout {
            coarse: plan(bed, opts)?,
            fine: None,
        });
    }
    let coarse_bed = bed.coarsen(2);
    let coarse = plan(&coarse_bed, opts)?;
    let (grid, dxc) = (bed.grid(), coarse_bed.header.dx);
    let shifted = bed.bathymetry(opts.tide);
    let depth = |i: usize, j: usize| -shifted.b[grid.idx(GHOST + i, GHOST + j)];
    let shallowest = |i: usize| {
        (0..grid.ny)
            .map(|j| depth(i, j))
            .fold(f64::INFINITY, f64::min)
    };
    let omega = std::f64::consts::TAU / opts.period;
    let wavelength = |h: f64| std::f64::consts::TAU / linear_wave(omega, h, opts.dispersive).k;

    let cells = wavelength(coarse.maker_depth) / dxc;
    if cells < CELLS_PER_WAVELENGTH {
        return Err(setup(format!(
            "on cells of {dxc} m a {} s swell is only {cells:.0} cells long even at the wave-maker, \
             too coarse to carry it ({CELLS_PER_WAVELENGTH} are needed); run on one grid with \
             --single",
            opts.period
        )));
    }
    // The depth where the swell is CELLS_PER_WAVELENGTH coarse cells long; wavelength grows
    // with depth, so step down from the maker's depth until it is shorter.
    let mut resolved = coarse.maker_depth;
    while resolved > 0.5 && wavelength(resolved - 0.1) >= CELLS_PER_WAVELENGTH * dxc {
        resolved -= 0.1;
    }
    let edge_depth = (FINE_DEPTH_PER_HEIGHT * opts.height).max(resolved);
    let lambda = wavelength(edge_depth);
    let zone_offshore = (lambda / grid.dx).ceil() as usize;

    // The fine grid must not reach back into the wave-maker's own bump.
    let near_field = ((coarse.maker_x + 6.0 * SIGMA_CELLS * dxc) / grid.dx).ceil() as usize;
    let edge = (near_field..grid.nx)
        .find(|&i| shallowest(i) < edge_depth)
        .ok_or_else(|| setup("the water never gets shallow enough for a reef to break on"))?;
    let x0 = edge.saturating_sub(zone_offshore).max(near_field);
    if grid.nx - x0 < 2 * zone_offshore {
        return Err(setup(format!(
            "the fine grid would be only {} cells long, with a relaxation zone of {zone_offshore}; \
             use a bed that reaches further onshore, or run on one grid with --single",
            grid.nx - x0
        )));
    }

    // The swell's leading edge leaves the maker at the start and reaches x0 at its group speed.
    let maker = (coarse.maker_x / grid.dx) as usize;
    let travel: f64 = (maker..x0)
        .map(|i| {
            grid.dx
                / linear_wave(omega, depth(i, grid.ny / 2).max(0.5), opts.dispersive).group_speed
        })
        .sum();
    let start = (travel - 2.0 * opts.period).max(0.0).min(coarse.duration);

    Ok(Layout {
        coarse,
        fine: Some(Fine {
            x0,
            edge_depth,
            zone_offshore,
            start,
        }),
    })
}

/// The beds the grids of `layout` run on: the coarse one, cropped, and the fine one.
///
/// Water seaward of the wave-maker deeper than the deepest point along its line is flattened to
/// that depth: it is only the sponge's, and on the Pipeline bed it was 39 m where the maker
/// stands in 27, which set the coarse grid's time step.
fn beds(bed: &Bed, layout: &Layout) -> (Bed, Option<Bed>) {
    let p = layout.coarse;
    let (mut coarse, fine) = match layout.fine {
        None => (bed.crop_x(p.crop), None),
        Some(f) => (bed.coarsen(2).crop_x(p.crop), Some(bed.crop_x(f.x0))),
    };
    let (nx, dx) = (coarse.header.nx, coarse.header.dx);
    let maker = ((p.maker_x / dx) as usize)
        .saturating_sub(p.crop)
        .min(nx - 1);
    let floor = coarse
        .elevation
        .chunks(nx)
        .map(|row| row[maker])
        .fold(f32::INFINITY, f32::min);
    for row in coarse.elevation.chunks_mut(nx) {
        for z in &mut row[..maker] {
            *z = z.max(floor);
        }
    }
    (coarse, fine)
}

/// A rough guess of how long a run will take in wall-clock seconds, so that a long one can be
/// refused before it starts. It counts each grid's cells and its steps at the time step of its
/// deepest water, at a cost per cell and step calibrated on one Apple M5: the nested Pipeline run
/// (319 x 200 cells of 6 m, then 179 x 400 of 3 m) took 5.7 minutes on 4 threads on the
/// efficiency cores, which are 3.5 times slower than the others. Treat it as the right order of
/// magnitude, not a promise; the wall-clock limit catches what it gets wrong.
pub fn estimate_seconds(bed: &Bed, opts: &RunOptions, layout: &Layout) -> f64 {
    let one_thread = if opts.dispersive { 650e-9 } else { 325e-9 };
    let speedup = (opts.threads.clamp(1, 6) as f64).powf(0.8);
    let cores = if opts.background { 3.5 } else { 1.0 };
    let stage = |bed: &Bed, seconds: f64| {
        let grid = bed.grid();
        let deepest = bed
            .elevation
            .iter()
            .fold(f64::INFINITY, |m, &z| m.min(f64::from(z)));
        let c = (G * (opts.tide - deepest).max(1.0)).sqrt();
        let dt = 0.4 / (c / grid.dx + c / grid.dy);
        (grid.nx * grid.ny) as f64 * seconds / dt
    };
    let (coarse, fine) = beds(bed, layout);
    let duration = layout.coarse.duration;
    let mut cell_steps = stage(&coarse, duration);
    if let (Some(fine), Some(f)) = (fine, layout.fine) {
        cell_steps += stage(&fine, duration - f.start);
    }
    cell_steps * one_thread / speedup * cores
}

/// One grid being run: its solver, its water, and the run directory it writes frames to.
struct Stage {
    bed: Bed,
    solver: Solver,
    state: State,
    writer: RunWriter,
    tide: f64,
    frame_interval: f64,
    next_frame: f64,
    frames: usize,
    steps: u64,
    max_rise: f64,
    eta: Vec<f32>,
    breaking: Vec<f32>,
}

impl Stage {
    fn new(
        bed: Bed,
        solver: Solver,
        out: &std::path::Path,
        settings: serde_json::Value,
        start: f64,
        opts: &RunOptions,
    ) -> Result<Self, Error> {
        let grid = bed.grid();
        let mut state = State::lake_at_rest(&grid, &solver.bed, 0.0);
        state.time = start;
        let writer = RunWriter::create(out, &bed, &["eta", "breaking"], settings)?;
        let cells = grid.nx * grid.ny;
        Ok(Self {
            bed,
            solver,
            state,
            writer,
            tide: opts.tide,
            frame_interval: opts.frame_interval,
            next_frame: (start / opts.frame_interval - 1e-9).ceil() * opts.frame_interval,
            frames: 0,
            steps: 0,
            max_rise: 0.0,
            eta: vec![0.0; cells],
            breaking: vec![0.0; cells],
        })
    }

    fn time(&self) -> f64 {
        self.state.time
    }

    /// Write a frame if one is due; the breaking field is returned for a look at it.
    fn frame_if_due(&mut self) -> Result<bool, Error> {
        if self.state.time < self.next_frame - 1e-9 {
            return Ok(false);
        }
        let grid = self.bed.grid();
        let mut rise = 0.0_f64;
        for (n, (i, j)) in grid.interior().enumerate() {
            let k = grid.idx(i, j);
            self.eta[n] = (self.state.h[k] + f64::from(self.bed.elevation[n])) as f32;
            if self.state.h[k] > 1e-3 {
                rise = rise.max(f64::from(self.eta[n]) - self.tide);
            }
        }
        if !rise.is_finite() || rise > 100.0 || self.eta.iter().any(|v| !v.is_finite()) {
            return Err(Error::Unstable {
                time: self.state.time,
                what: format!("surface {rise:.1} m above still water"),
            });
        }
        self.max_rise = self.max_rise.max(rise);
        for (b, v) in self
            .breaking
            .iter_mut()
            .zip(self.solver.breaking_indicator(&self.state))
        {
            *b = v as f32;
        }
        self.writer
            .write_frame(self.state.time, &[&self.eta, &self.breaking])?;
        self.frames += 1;
        self.next_frame += self.frame_interval;
        Ok(true)
    }

    /// The next step: as long as is stable, ending no later than `until` or the next frame.
    fn next_dt(&self, until: f64) -> Result<f64, Error> {
        let stable = self.solver.stable_dt(&self.state);
        if !stable.is_finite() {
            return Err(setup("no water in the domain at this tide"));
        }
        Ok(stable
            .min(self.next_frame - self.state.time)
            .min(until - self.state.time))
    }

    fn advance(&mut self, dt: f64) {
        self.solver.step(&mut self.state, dt);
        self.steps += 1;
    }

    /// One step of [`Stage::next_dt`].
    fn step(&mut self, until: f64) -> Result<(), Error> {
        let dt = self.next_dt(until)?;
        self.advance(dt);
        Ok(())
    }

    fn finish(
        mut self,
        truncated: bool,
        notes: &[(&str, serde_json::Value)],
    ) -> Result<PathBuf, Error> {
        if truncated {
            self.writer.note("truncated", serde_json::json!(true));
            self.writer
                .note("stopped_at_s", serde_json::json!(self.state.time));
        }
        for (key, value) in notes {
            self.writer.note(key, value.clone());
        }
        Ok(self.writer.finish()?)
    }
}

fn solver_for(bed: &Bed, opts: &RunOptions) -> Solver {
    let solver = Solver::new(bed.grid(), bed.bathymetry(opts.tide)).with_manning(opts.manning);
    if opts.dispersive {
        solver
            .with_order(Order::Fifth)
            .with_dispersion(Dispersion::default())
    } else {
        solver
    }
}

fn model(opts: &RunOptions) -> &'static str {
    if opts.dispersive {
        "Serre-Green-Naghdi (flat-bed operator), hybrid breaking, fifth-order upwind-biased scheme"
    } else {
        "nonlinear shallow water, MC limiter, non-dispersive"
    }
}

/// Surface elevation and velocity `(eta, u, v)` that a cell of a relaxation zone is pulled
/// towards.
type Target = (f64, f64, f64);

/// The coarse run's surface and velocity at each cell of the fine grid's relaxation zones,
/// interpolated bilinearly from the wet coarse cells around it (dry ones would put the bed's
/// height in as the surface), or no water where none of them is wet.
struct Sampler {
    /// Per zone cell: four padded coarse indices and their weights.
    corners: Vec<[(usize, f64); 4]>,
}

impl Sampler {
    fn new(zone: &Relaxation, fine: &Bed, f: &Fine, coarse: &Bed, crop: usize) -> Self {
        let (cg, fg) = (coarse.grid(), fine.grid());
        let corners = zone
            .cells()
            .map(|(i, j)| {
                // Centre of the fine cell in the bed file's frame, then in coarse cells.
                let x = (f.x0 + i - GHOST) as f64 * fg.dx + 0.5 * fg.dx;
                let y = (j - GHOST) as f64 * fg.dy + 0.5 * fg.dy;
                let cx = x / cg.dx - crop as f64 - 0.5;
                let cy = y / cg.dy - 0.5;
                let (ix, iy) = (cx.floor(), cy.floor());
                let (wx, wy) = (cx - ix, cy - iy);
                let at = |di: f64, dj: f64| {
                    let pi = ((ix + di) as isize + GHOST as isize)
                        .clamp(GHOST as isize - 1, (cg.nx + GHOST) as isize)
                        as usize;
                    let pj = ((iy + dj) as isize + GHOST as isize)
                        .clamp(GHOST as isize - 1, (cg.ny + GHOST) as isize)
                        as usize;
                    cg.idx(pi, pj)
                };
                [
                    (at(0.0, 0.0), (1.0 - wx) * (1.0 - wy)),
                    (at(1.0, 0.0), wx * (1.0 - wy)),
                    (at(0.0, 1.0), (1.0 - wx) * wy),
                    (at(1.0, 1.0), wx * wy),
                ]
            })
            .collect();
        Self { corners }
    }

    /// `(eta, u, v)` per zone cell from the coarse water `s` over the coarse bed `b`, whose
    /// ghost cells must be filled.
    fn sample(&self, s: &State, b: &[f64]) -> Vec<Target> {
        self.corners
            .iter()
            .map(|corners| {
                let (mut eta, mut u, mut v, mut weight) = (0.0, 0.0, 0.0, 0.0);
                for &(k, w) in corners {
                    if s.h[k] > 1e-3 && w > 0.0 {
                        eta += w * (s.h[k] + b[k]);
                        u += w * s.hu[k] / s.h[k];
                        v += w * s.hv[k] / s.h[k];
                        weight += w;
                    }
                }
                if weight > 0.0 {
                    (eta / weight, u / weight, v / weight)
                } else {
                    (f64::NEG_INFINITY, 0.0, 0.0)
                }
            })
            .collect()
    }
}

/// Write the coarse grid's frame if one is due, and count it if the swell breaks there under the
/// fine grid's offshore zone, in the middle half of the width (`columns`, for a nested run).
fn coarse_frame(
    coarse: &mut Stage,
    columns: Option<(usize, usize)>,
    breaking: &mut usize,
) -> Result<bool, Error> {
    if !coarse.frame_if_due()? {
        return Ok(false);
    }
    if let Some((first, last)) = columns {
        let g = coarse.bed.grid();
        if (g.ny / 4..g.ny - g.ny / 4).any(|j| {
            coarse.breaking[j * g.nx + first..j * g.nx + last]
                .iter()
                .any(|&v| v >= 0.8)
        }) {
            *breaking += 1;
        }
    }
    Ok(true)
}

/// How far the run has got: the coarse grid's time, and the main run's frames.
fn progress(coarse: &Stage, fine: Option<&Stage>, duration: f64) -> Progress {
    let main = fine.unwrap_or(coarse);
    Progress {
        time: coarse.time(),
        duration,
        steps: coarse.steps + fine.map_or(0, |f| f.steps),
        frames: main.frames,
        max_rise: main.max_rise,
    }
}

/// Run the simulation and write the run directory (and, for a nested run, the coarse run in
/// its `coarse` subdirectory). `on_progress` is called after every frame of the main run.
pub fn run(opts: &RunOptions, mut on_progress: impl FnMut(&Progress)) -> Result<Summary, Error> {
    let started = Instant::now();
    let full = Bed::read(&opts.bed)?;
    let layout = layout(&full, opts)?;
    let estimate = estimate_seconds(&full, opts, &layout);
    if let Some(budget) = opts.budget_seconds
        && estimate > budget
    {
        return Err(Error::OverBudget {
            estimate,
            budget,
            threads: opts.threads,
            background: opts.background,
        });
    }
    let plan = layout.coarse;
    let (coarse_bed, fine_bed) = beds(&full, &layout);
    let grid = coarse_bed.grid();
    let maker_x = plan.maker_x - plan.crop as f64 * grid.dx;
    let maker = WaveMaker {
        x: maker_x,
        amplitude: opts.height / 2.0,
        period: opts.period,
        angle: 0.0,
        sigma: SIGMA_CELLS * grid.dx,
        ramp_periods: RAMP_PERIODS,
    };
    let solver = solver_for(&coarse_bed, opts)
        .with_wavemaker(maker)
        .with_sponge(
            Sponge::new(plan.sponge_offshore, SPONGE_STRENGTH).edges(true, false, false, false),
        );
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
        // The sides are walls, which mirror the bed: on the Pipeline bed the middle half of the
        // width was right to within 4% against a domain nearly twice as wide, so a viewer
        // should show that and trim this many rows on each side.
        "side_margin_cells": grid.ny / 4,
        "manning": opts.manning,
        "direction": "along +x of the bed frame",
        "model": model(opts),
        "dispersive": opts.dispersive,
    });
    let coarse_out = match layout.fine {
        Some(_) => opts.out.join("coarse"),
        None => opts.out.clone(),
    };
    let mut coarse = Stage::new(coarse_bed, solver, &coarse_out, settings, 0.0, opts)?;

    // The fine grid, its relaxation zones, and the coarse water they are pulled towards.
    let mut nest = match (layout.fine, fine_bed) {
        (Some(f), Some(bed)) => {
            let fg = bed.grid();
            let zone = Relaxation::new(
                &fg,
                &[Sponge::new(f.zone_offshore, SPONGE_STRENGTH).edges(true, false, false, false)],
            );
            let sampler = Sampler::new(&zone, &bed, &f, &coarse.bed, plan.crop);
            let settings = serde_json::json!({
                "wave_height_m": opts.height,
                "wave_period_s": opts.period,
                "tide_m": opts.tide,
                "duration_s": plan.duration,
                "frame_interval_s": opts.frame_interval,
                "cropped_offshore_m": f.x0 as f64 * fg.dx,
                "side_margin_cells": fg.ny / 4,
                "nested_in": "coarse",
                "coarse_cell_m": coarse.bed.header.dx,
                "fine_start_s": f.start,
                "edge_depth_m": f.edge_depth,
                "zone_offshore_cells": f.zone_offshore,
                // The relaxation zone is the coarse run, not this one: a viewer should start
                // past it.
                "near_field_end_m": f.zone_offshore as f64 * fg.dx,
                "sponge_offshore_cells": f.zone_offshore,
                "manning": opts.manning,
                "direction": "along +x of the bed frame",
                "model": model(opts),
                "dispersive": opts.dispersive,
            });
            let solver = solver_for(&bed, opts);
            let stage = Stage::new(bed, solver, &opts.out, settings, f.start, opts)?;
            Some((stage, zone, sampler, f))
        }
        _ => None,
    };

    let mut truncated = false;
    let mut breaking_in_zone = 0_usize;
    let mut fine_seconds = 0.0;
    // Columns of the coarse grid under the fine grid's offshore relaxation zone.
    let zone_columns = layout.fine.map(|f| {
        let ratio = full.header.dx / coarse.bed.header.dx;
        let first = (f.x0 as f64 * ratio) as usize - plan.crop;
        let last = ((f.x0 + f.zone_offshore) as f64 * ratio).ceil() as usize - plan.crop;
        (first, last.min(coarse.bed.header.nx))
    });
    // The coarse water just before its latest step, and the zone's targets sampled at both ends
    // of that step.
    let mut before: Option<State> = None;
    let mut targets: Option<(f64, Vec<Target>, f64, Vec<Target>)> = None;
    loop {
        if opts
            .max_wall_seconds
            .is_some_and(|limit| started.elapsed().as_secs_f64() > limit)
        {
            truncated = true;
            break;
        }
        let fine_active = nest
            .as_ref()
            .is_some_and(|(_, _, _, f)| coarse.time() >= f.start - 1e-9);
        if !fine_active {
            if coarse_frame(&mut coarse, zone_columns, &mut breaking_in_zone)? {
                on_progress(&progress(
                    &coarse,
                    nest.as_ref().map(|n| &n.0),
                    plan.duration,
                ));
            }
            if coarse.time() >= plan.duration - 1e-9 {
                break;
            }
            before = Some(coarse.state.clone());
            coarse.step(plan.duration)?;
            continue;
        }
        let Some((fine, zone, sampler, _)) = nest.as_mut() else {
            unreachable!("the fine grid is active only in a nested run")
        };
        fine.frame_if_due()?;
        if fine.time() >= plan.duration - 1e-9 {
            coarse_frame(&mut coarse, zone_columns, &mut breaking_in_zone)?;
            break;
        }
        // The fine grid's own step, with the coarse run ahead of it to bracket its end.
        let dt = fine.next_dt(plan.duration)?;
        let end = fine.time() + dt;
        while coarse.time() < end - 1e-9 {
            if coarse_frame(&mut coarse, zone_columns, &mut breaking_in_zone)? {
                on_progress(&progress(&coarse, Some(fine), plan.duration));
            }
            before = Some(coarse.state.clone());
            coarse.step(plan.duration)?;
            targets = None;
        }
        let fine_started = Instant::now();
        if targets.is_none() {
            let cg = coarse.bed.grid();
            let mut then = before.clone().expect("the coarse run has stepped");
            let mut now = coarse.state.clone();
            then.fill_boundaries(&cg, false);
            now.fill_boundaries(&cg, false);
            let b = &coarse.solver.bed.b;
            targets = Some((
                then.time,
                sampler.sample(&then, b),
                now.time,
                sampler.sample(&now, b),
            ));
        }
        let (t0, from, t1, to) = targets.as_ref().expect("sampled above");
        let w = ((end - t0) / (t1 - t0)).clamp(0.0, 1.0);
        let lerp = |m: usize| {
            let (a, c) = (from[m], to[m]);
            if a.0.is_finite() && c.0.is_finite() {
                (
                    a.0 + w * (c.0 - a.0),
                    a.1 + w * (c.1 - a.1),
                    a.2 + w * (c.2 - a.2),
                )
            } else if w < 0.5 {
                a
            } else {
                c
            }
        };
        fine.advance(dt);
        let fg = fine.bed.grid();
        zone.apply(&fg, &fine.solver.bed, &mut fine.state, dt, lerp);
        fine_seconds += fine_started.elapsed().as_secs_f64();
    }

    let mut warnings = Vec::new();
    if breaking_in_zone > 0 {
        warnings.push(format!(
            "the coarse run broke inside the fine grid's offshore relaxation zone in \
             {breaking_in_zone} frames, in the middle half of the width: the fine grid should start \
             in deeper water"
        ));
    }
    let coarse_steps = coarse.steps;
    let coarse_dispersion = coarse.solver.dispersion_stats();
    let (main_out, coarse_header) = match nest {
        Some((fine, ..)) => {
            let coarse_header = coarse.finish(truncated, &[])?;
            let (steps, frames, max_rise) = (fine.steps, fine.frames, fine.max_rise);
            let dispersion = fine.solver.dispersion_stats();
            let notes = [(
                "breaking_in_offshore_zone_frames",
                serde_json::json!(breaking_in_zone),
            )];
            let header = fine.finish(truncated, &notes)?;
            (
                (header, steps, frames, max_rise, dispersion),
                Some(coarse_header),
            )
        }
        None => {
            let (frames, max_rise) = (coarse.frames, coarse.max_rise);
            let header = coarse.finish(truncated, &[])?;
            ((header, 0, frames, max_rise, coarse_dispersion), None)
        }
    };
    let (header, fine_steps, frames, max_rise, dispersion) = main_out;
    let wall_seconds = started.elapsed().as_secs_f64();
    Ok(Summary {
        steps: coarse_steps + fine_steps,
        coarse_steps,
        coarse_seconds: wall_seconds - fine_seconds,
        fine_steps,
        fine_seconds,
        frames,
        wall_seconds,
        max_rise,
        header,
        coarse_header,
        layout,
        dispersion,
        truncated,
        warnings,
    })
}
