//! `wavesim slice`: one wave breaking over a made-up seabed, seen side-on, computed by the slice
//! solver (`waveslice`) and written for the viewer's side view.
//!
//! The wave is a solitary wave, one hump on still water, a stand-in for one big wave of a set:
//! the published breaking computations use it because its height alone defines it. It starts
//! over flat water and runs onto one of a few seabeds, the same wave on each, so the difference
//! between how they break comes from the bed alone. The run ends when the lip lands.

use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::PathBuf;
use std::time::Instant;

use waveio::{Curl, Landing, SliceWriter};
use waveslice::{Solitary, Tank};

use crate::{Error, human};

const G: f64 = 9.81;
/// Steps timed before the run is estimated, and refused if it would take too long.
const CALIBRATION_STEPS: usize = 5;
/// How many more steps a run takes than its first steps' length and [`travel_time`] suggest, as
/// steps shorten while the wave shoals and its lip speeds up: measured on the steep case, 1083
/// steps where they suggested 347.
const SLOWDOWN: f64 = 3.1;

/// The made-up seabeds, all starting with fifteen depths of flat water for the wave to start
/// on, its crest seven and a half depths from the left wall.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Case {
    /// A plane slope rising 1 in 22: the wave shoals for longer and breaks less hard. (1 in 30
    /// would take hours at this resolution: the slope is long, the nodes along it many, and the
    /// solve's cost grows with the cube of their number.)
    Gentle,
    /// A plane slope rising 1 in 15: Grilli et al.'s plunging case, tested in `waveslice`.
    Steep,
    /// A reef: a 1 in 10 ramp up onto a flat shelf 0.12 of the offshore depth deep, six depths
    /// long, so the wave breaks at the reef's edge and its lip lands on the shelf, far from the
    /// wall. (A first version, a 1 in 3 ledge onto a shelf 0.3 of the depth deep, broke nothing:
    /// no solitary wave breaks on a slope steeper than about 12 degrees (Grilli et al. 1997),
    /// and the wave, longer than the ledge, stepped up onto the shelf whole and only curled
    /// against the wall at its far end.)
    Reef,
}

impl Case {
    pub fn name(self) -> &'static str {
        match self {
            Case::Gentle => "gentle",
            Case::Steep => "steep",
            Case::Reef => "reef",
        }
    }

    /// The seabed as `(x, depth)` points, for water `depth` deep offshore.
    pub fn profile(self, depth: f64) -> Vec<(f64, f64)> {
        let toe = 15.0 * depth;
        let mut p = vec![(0.0, depth), (toe, depth)];
        match self {
            Case::Gentle => p.push((toe + 22.0 * 0.95 * depth, 0.05 * depth)),
            // The steep beach flattens into a shelf 0.05 of the depth deep, three depths long,
            // before the wall: without it the lip lands a metre from the wall. The shelf changed
            // where it lands by 2 cm and nothing else, so the wall had not shaped it, but now
            // nothing is left to doubt. The gentle beach's lip lands 23 m short of its wall.
            Case::Steep => {
                let shallows = toe + 15.0 * 0.95 * depth;
                p.push((shallows, 0.05 * depth));
                p.push((shallows + 3.0 * depth, 0.05 * depth));
            }
            Case::Reef => {
                let edge = toe + 10.0 * 0.88 * depth;
                p.push((edge, 0.12 * depth));
                p.push((edge + 6.0 * depth, 0.12 * depth));
            }
        }
        p
    }
}

#[derive(Clone, Debug)]
pub struct SliceOptions {
    pub case: Case,
    /// Offshore water depth, m.
    pub depth: f64,
    /// The wave's height above still water, m.
    pub height: f64,
    /// Surface nodes no closer than this where the water is shallow, m.
    pub finest: f64,
    /// A frame every this many steps.
    pub frame_every: usize,
    pub out: PathBuf,
    /// Refuse a run estimated to take longer, and stop one that does, keeping what it has.
    pub max_wall_seconds: Option<f64>,
    /// Carry the wave with the main model until its crest is over water this deep, then hand
    /// it to the slice solver (see [`crate::handoff`]); `None` runs the slice solver all the way.
    pub handoff_depth: Option<f64>,
    /// A real line of seabed (a `wavesim-line` file) to break a swell on instead of `case`'s
    /// made-up one, with a swell `height` from crest to trough and `period` long; it needs
    /// `handoff_depth`.
    pub line: Option<PathBuf>,
    pub period: f64,
    /// For `line`: the wave a 2D run saved (`wavesim run --handoff-y`), instead of bringing a
    /// swell in along the line alone, which cannot spread energy sideways as the reef does.
    pub handoff_file: Option<PathBuf>,
}

/// What a finished slice run found.
#[derive(Clone, Debug)]
pub struct SliceSummary {
    pub dir: PathBuf,
    pub ended: String,
    pub curl: Option<Curl>,
    pub landing: Option<Landing>,
    pub steps: usize,
    pub frames: usize,
    pub seconds: f64,
    pub estimate: f64,
}

/// The depth of `profile` at `x`, straight between its points.
fn depth_at(profile: &[(f64, f64)], x: f64) -> f64 {
    let k = profile
        .partition_point(|p| p.0 <= x)
        .clamp(1, profile.len() - 1);
    let (a, b) = (profile[k - 1], profile[k]);
    a.1 + (b.1 - a.1) * ((x - a.0) / (b.0 - a.0)).clamp(0.0, 1.0)
}

/// Seconds for a wave `height` high to travel from `from` at the long wave speed
/// `sqrt(g (d + height))` to where the water is half as deep as the wave is high, about where
/// these waves break, or to the end of `profile`.
fn travel_time(profile: &[(f64, f64)], from: f64, height: f64) -> f64 {
    let end = profile[profile.len() - 1].0;
    let n = 1000;
    let dx = (end - from) / n as f64;
    (0..n)
        .map(|k| from + (k as f64 + 0.5) * dx)
        .take_while(|&x| depth_at(profile, x) > 0.5 * height)
        .map(|x| dx / (G * (depth_at(profile, x) + height)).sqrt())
        .sum()
}

/// The slice solver's tank after the main model has carried `wave` over `profile` until its
/// crest is over water `depth` deep: the water from seven and a half offshore depths behind the
/// crest to the shore, the profile cut there, and what the hand-off was, for the header.
fn handed(
    opts: &SliceOptions,
    wave: &Solitary,
    profile: Vec<(f64, f64)>,
    depth: f64,
) -> Result<Handed, Error> {
    let end = profile[profile.len() - 1].0;
    let still = |x: f64| depth_at(&profile, x);
    let c = wave.speed();
    // The main model's cells half the slice's finest node spacing.
    let cell = 0.5 * opts.finest;
    let approach = crate::handoff::approach(
        end,
        cell,
        &still,
        &|x| wave.eta(x),
        &|x| c * wave.eta(x) / (still(x) + wave.eta(x)),
        depth,
    )?;
    let start = (approach.crest - 7.5 * opts.depth).max(0.0);
    let mut cut = vec![(start, still(start))];
    cut.extend(profile.iter().copied().filter(|p| p.0 > start + 1e-9));
    let mut tank = Tank::shaped(
        G,
        &cut,
        opts.finest,
        0.1 * opts.depth,
        &|x| approach.eta_at(x),
        &|xs| approach.potential(xs, &still),
    );
    tank.time = approach.time;
    let record = serde_json::json!({
        "depth": depth,
        "time": approach.time,
        "crest": approach.crest,
        "start": start,
        "cell": cell,
        "approach_steps": approach.steps,
    });
    Ok(Handed {
        tank,
        profile: cut,
        record,
    })
}

/// A `wavesim-line` file, as `tools/fetch_line.py` writes it.
#[derive(serde::Deserialize)]
struct LineFile {
    name: String,
    x: Vec<f64>,
    z: Vec<Option<f64>>,
    z_bed: Vec<f64>,
}

/// The slice solver's tank for one wave of a swell on a real line: the main model brings the
/// swell in along the whole line (the fine data where it has any, the bed file's where not)
/// until a crest passes the last point `depth` deep, and that one wave is handed over (see
/// [`crate::handoff::Approach::single_wave`]) to a slice from a little behind it to a wall where
/// the water is half a metre deep, or two wavelengths ahead if that is nearer.
/// Where the wave on a real line comes from.
enum Source {
    /// A 2D run's saved wave: the sea as the reef shapes it, refraction and all.
    Run(PathBuf),
    /// A swell brought in along the line alone by the main model, handed over at this depth.
    /// The line cannot spread energy sideways: on the Pipeline line at y = 511.5 m the crest
    /// it hands over at 9 m is three times as high as the 2D run's there.
    Line(f64),
}

fn line_tank(opts: &SliceOptions, path: &std::path::Path, source: Source) -> Result<Handed, Error> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| Error::Setup(format!("{}: {e}", path.display())))?;
    let line: LineFile = serde_json::from_str(&text)
        .map_err(|e| Error::Setup(format!("{}: not a wavesim-line: {e}", path.display())))?;
    let points: Vec<(f64, f64)> = line
        .x
        .iter()
        .zip(line.z.iter().zip(&line.z_bed))
        .map(|(&x, (z, &bed))| (x, z.unwrap_or(bed)))
        .collect();
    let still = |x: f64| -> f64 {
        let k = points
            .partition_point(|p| p.0 <= x)
            .clamp(1, points.len() - 1);
        let (a, b) = (points[k - 1], points[k]);
        -(a.1 + (b.1 - a.1) * ((x - a.0) / (b.0 - a.0)).clamp(0.0, 1.0))
    };
    // The main model's cells: 2 m is about 60 to a wavelength of a 14 s swell on the reef.
    let cell = 2.0;
    let (approach, depth, from) = match &source {
        Source::Run(file) => {
            let (a, depth, broke) = crate::handoff::Approach::read(file)?;
            if let Some(x) = broke {
                return Err(Error::Setup(format!(
                    "the 2D run's waves broke at x = {x:.0} m, seaward of where this one was saved \
                     (crest at {:.0} m), so it had already broken once: a slice would start from \
                     a broken wave",
                    a.crest
                )));
            }
            (a, depth, file.display().to_string())
        }
        Source::Line(depth) => (
            crate::handoff::swell(&points, cell, opts.height, opts.period, *depth)?,
            *depth,
            "the line alone".to_string(),
        ),
    };
    let eta_ratio = approach.eta_ratio(&still);
    // A 2D run's wave must be handed over before the run's breaking switch has begun to act on
    // its crest, which it does once the crest is 0.30 of the total depth beneath it: past that
    // the run has turned its dispersion down, and the slice would start from a wave already
    // part broken. On Pipeline at y = 421.5 m a 2.5 m, 16 s swell handed over with its crest
    // at 0.37 closed a tube 2.6 m tall, four times the tube from a hand-off 1 m deeper.
    // (A swell brought in along the line alone has the switch off, so it is not refused.)
    let onset = wavecore::Breaking::default().eta_ratio.0;
    if matches!(source, Source::Run(_)) && eta_ratio / (1.0 + eta_ratio) > onset {
        return Err(Error::Setup(format!(
            "the saved crest is {:.2} of the total depth beneath it, past the {onset:.2} where \
             the 2D run's breaking switch starts to act, so the wave was already part broken \
             when saved; hand off at a deeper --handoff-depth",
            eta_ratio / (1.0 + eta_ratio)
        )));
    }
    let omega = std::f64::consts::TAU / opts.period;
    let wavelength =
        std::f64::consts::TAU / wavecore::linear_wave(omega, still(approach.crest), true).k;
    let taper = 0.15 * wavelength;
    let (wave, behind, ahead) = approach.single_wave(wavelength, taper);
    let start = behind - taper - 10.0;
    // The saved line must reach back past the trough behind the wave and its taper: a 2D run
    // saves only its fine grid, and a hand-off point too close to where that starts leaves the
    // trough outside it, where the lowest point found would only be the edge of the data.
    let first = wave.x[0];
    let spacing = wave.x[1] - wave.x[0];
    if start < first || behind - first < 2.0 * spacing {
        return Err(Error::Setup(format!(
            "the saved line starts at x = {first:.0} m, but the slice must begin at {start:.0} m to \
             hold the trough behind the wave (at {behind:.0} m) and its taper; hand off further \
             shoreward (a shallower --handoff-depth), or on a line whose hand-off point lies \
             deeper inside the run's fine grid"
        )));
    }
    // The wall at the shallow end: where the water is half a metre deep, or two wavelengths
    // ahead of the handed crest if that is nearer. Every metre of shallow water costs nodes at
    // the finest spacing, but one wavelength was too short: on the Pipeline line at y = 511.5 m
    // the wave broke on the edge of the reef flat and its lip landed 2 m from a wall there. A
    // lip that lands near the wall anyway is reported as shaped by it.
    let shallow = points
        .iter()
        .find(|p| p.0 > approach.crest && -p.1 < 0.5)
        .map(|p| p.0)
        .ok_or_else(|| Error::Setup("the line never reaches water shallower than 0.5 m".into()))?;
    let wall = shallow.min(approach.crest + 2.0 * wavelength);
    let raw: Vec<(f64, f64)> = std::iter::once((start, still(start)))
        .chain(
            points
                .iter()
                .filter(|p| p.0 > start && p.0 < wall)
                .map(|p| (p.0, -p.1)),
        )
        .chain(std::iter::once((wall, still(wall))))
        .collect();
    let cut = simplify(&raw, 0.1);
    let mut tank = Tank::shaped(
        G,
        &cut,
        opts.finest,
        0.1 * cut[0].1,
        &|x| wave.eta_at(x),
        &|xs| wave.potential(xs, &still),
    );
    tank.time = wave.time;
    let record = serde_json::json!({
        "line": line.name,
        "from": from,
        "swell_height": opts.height,
        "swell_period": opts.period,
        "depth": depth,
        "time": wave.time,
        "crest": wave.crest,
        "crest_over_depth": eta_ratio,
        "trough_behind": behind,
        "trough_ahead": ahead,
        "wavelength": wavelength,
        "start": start,
        "cell": cell,
        "approach_steps": wave.steps,
        "bed_points": cut.len(),
    });
    Ok(Handed {
        tank,
        profile: cut,
        record,
    })
}

/// The points of `p` that keep it within `tolerance` of the original (Ramer-Douglas-Peucker),
/// first and last kept: a metre-by-metre line of lidar comes down to the points that shape it.
fn simplify(p: &[(f64, f64)], tolerance: f64) -> Vec<(f64, f64)> {
    let mut keep = vec![false; p.len()];
    keep[0] = true;
    keep[p.len() - 1] = true;
    let mut spans = vec![(0, p.len() - 1)];
    while let Some((a, b)) = spans.pop() {
        let (pa, pb) = (p[a], p[b]);
        let off = |q: (f64, f64)| {
            let t = (q.0 - pa.0) / (pb.0 - pa.0);
            (q.1 - (pa.1 + t * (pb.1 - pa.1))).abs()
        };
        if let Some((k, d)) = (a + 1..b)
            .map(|k| (k, off(p[k])))
            .max_by(|x, y| x.1.total_cmp(&y.1))
            && d > tolerance
        {
            keep[k] = true;
            spans.push((a, k));
            spans.push((k, b));
        }
    }
    p.iter()
        .zip(keep)
        .filter(|(_, k)| *k)
        .map(|(q, _)| *q)
        .collect()
}

/// A tank after the hand-off, its cut profile, and the hand-off's record for the header.
struct Handed {
    tank: Tank,
    profile: Vec<(f64, f64)>,
    record: serde_json::Value,
}

/// Run one wave until its lip lands, it reaches the shallow wall, the time limit passes, or the
/// computation fails, writing frames as it goes. A run its first steps say would take longer
/// than the limit is refused before anything is written.
pub fn run_slice(opts: &SliceOptions) -> Result<SliceSummary, Error> {
    let profile = opts.case.profile(opts.depth);
    if opts.line.is_none() && (opts.height <= 0.0 || opts.height > 0.6 * opts.depth) {
        return Err(Error::Setup(format!(
            "a wave {} m high in {} m of water: a solitary wave can be up to about 0.6 of the \
             depth",
            opts.height, opts.depth
        )));
    }
    let crest = 7.5 * opts.depth;
    let wave = Solitary {
        g: G,
        depth: opts.depth,
        height: opts.height,
        crest,
    };
    let (mut tank, profile, handoff) = match (&opts.line, opts.handoff_depth) {
        (Some(path), depth) => {
            let source = match (&opts.handoff_file, depth) {
                (Some(file), _) => Source::Run(file.clone()),
                (None, Some(depth)) => Source::Line(depth),
                (None, None) => {
                    return Err(Error::Setup(
                        "a real line needs --handoff-file (a wave saved by `wavesim run \
                         --handoff-y`) or --handoff-depth (a swell brought in along the line \
                         alone)"
                            .into(),
                    ));
                }
            };
            let h = line_tank(opts, path, source)?;
            (h.tank, h.profile, h.record)
        }
        (None, None) => (
            Tank::over(&wave, &profile, opts.finest),
            profile,
            serde_json::Value::Null,
        ),
        (None, Some(depth)) => {
            let h = handed(opts, &wave, profile, depth)?;
            (h.tank, h.profile, h.record)
        }
    };
    let speeds = |t: &Tank| -> Vec<f64> { t.velocity().iter().map(|v| v.0.hypot(v.1)).collect() };
    let start = (tank.surface.clone(), speeds(&tank));
    let start_time = tank.time;

    // Time the first steps and refuse a run that would go over its limit.
    let clock = Instant::now();
    for _ in 0..CALIBRATION_STEPS {
        tank.march();
    }
    let per_step = clock.elapsed().as_secs_f64() / CALIBRATION_STEPS as f64;
    let mean_dt = (tank.time - start_time) / CALIBRATION_STEPS as f64;
    let from = start
        .0
        .iter()
        .copied()
        .fold((0.0, f64::MIN), |a, p| if p.1 > a.1 { p } else { a });
    let steps = travel_time(&profile, from.0, opts.height) / mean_dt * SLOWDOWN;
    let estimate = steps * per_step;
    if let Some(limit) = opts.max_wall_seconds
        && estimate > limit
    {
        return Err(Error::Setup(format!(
            "this slice would take about {}, more than its limit of {}; allow it with a larger \
             --max-minutes, or space the surface nodes further apart with a larger --finest",
            human(estimate),
            human(limit)
        )));
    }

    let settings = match &opts.line {
        // A real line: the wave is one of a swell, its height and depth the swell's and the
        // line's, recorded with the hand-off.
        Some(_) => serde_json::json!({
            "case": "line",
            "depth": profile[0].1,
            "g": G,
            "finest": opts.finest,
            "profile": profile,
            "wave": "swell",
            "handoff": handoff,
        }),
        None => serde_json::json!({
            "case": opts.case.name(),
            "depth": opts.depth,
            "height": opts.height,
            "crest": crest,
            "g": G,
            "finest": opts.finest,
            "profile": profile,
            "wave": "solitary",
            "handoff": handoff,
        }),
    };
    let mut out = SliceWriter::create(&opts.out, &tank.bed, tank.surface.len(), settings)?;
    out.write_frame(start_time, &start.0, &start.1)?;
    let (mut curl, mut landing) = (None, None);
    let ended = loop {
        if let Some((i, j)) = tank.landing() {
            let s = &tank.surface;
            let (x, z) = (0.5 * (s[i].0 + s[j].0), 0.5 * (s[i].1 + s[j].1));
            let tube = tank.tube().expect("a landed lip closes a tube");
            landing = Some(Landing {
                time: tank.time,
                x,
                z,
                throw: curl.map_or(0.0, |c: Curl| x - c.x),
                tube_area: tube.area,
                tube_width: tube.width,
                tube_height: tube.height,
            });
            out.write_frame(tank.time, &tank.surface, &speeds(&tank))?;
            // A lip that lands against the wall at the shallow end was shaped by the wall.
            let wall = profile[profile.len() - 1].0;
            if x > wall - 2.0 * opts.height {
                break "landed against the wall, which shaped it: not a clean break".to_string();
            }
            break "landed".to_string();
        }
        // The wave has reached the shallow wall when the water there has risen by half its
        // height; it breaks well before that if it breaks at all.
        if tank.surface[tank.surface.len() - 1].1 > 0.5 * opts.height {
            break "wall".to_string();
        }
        if tank.tangled() {
            break "failed: neighbouring surface nodes crossed".to_string();
        }
        if let Some(limit) = opts.max_wall_seconds
            && clock.elapsed().as_secs_f64() > limit
        {
            break "time limit".to_string();
        }
        if catch_unwind(AssertUnwindSafe(|| tank.march())).is_err() {
            break "failed: the solver could not solve".to_string();
        }
        if !tank
            .surface
            .iter()
            .all(|p| p.0.is_finite() && p.1.is_finite())
        {
            break "failed: the surface is no longer finite".to_string();
        }
        let first_curl = curl.is_none() && tank.overturned();
        if first_curl {
            let s = &tank.surface;
            let i = (0..s.len() - 1)
                .find(|&i| s[i + 1].0 < s[i].0)
                .expect("overturned");
            let top = s.iter().map(|p| p.1).fold(f64::MIN, f64::max);
            curl = Some(Curl {
                time: tank.time,
                x: s[i].0,
                crest: top,
                depth: depth_at(&profile, s[i].0),
            });
        }
        if first_curl || tank.steps.is_multiple_of(opts.frame_every) {
            out.write_frame(tank.time, &tank.surface, &speeds(&tank))?;
        }
    };
    let frames = out.frame_count();
    let dir = out.finish(&ended, curl, landing)?;
    Ok(SliceSummary {
        dir: dir.parent().map(PathBuf::from).unwrap_or_default(),
        ended,
        curl,
        landing,
        steps: tank.steps,
        frames,
        seconds: clock.elapsed().as_secs_f64(),
        estimate,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A slice of a 16 s swell from a saved wave on a plane beach rising 1 in 40, its crest
    /// `crest` m above still water where the beach is 5 m deep, the run's waves having broken
    /// at `broke` if they did; refused at once by a limit of no time at all if nothing refuses
    /// it first.
    fn slice_of_saved_wave(name: &str, crest: f64, broke: Option<f64>) -> String {
        let dir = std::env::temp_dir().join(format!("wavesim-slice-{name}-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let x: Vec<f64> = (0..=400).map(f64::from).collect();
        let depth = |x: f64| 10.0 - x / 40.0;
        let z: Vec<f64> = x.iter().map(|&x| -depth(x)).collect();
        let line = dir.join("line.json");
        let doc = serde_json::json!({ "name": "plane", "x": x, "z": z, "z_bed": z });
        std::fs::write(&line, doc.to_string()).unwrap();
        let wavelength = 107.0;
        let eta: Vec<f64> = x
            .iter()
            .map(|&x| crest * (std::f64::consts::TAU * (x - 200.0) / wavelength).cos())
            .collect();
        let ubar: Vec<f64> = x
            .iter()
            .zip(&eta)
            .map(|(&x, &e)| {
                if depth(x) > 0.0 {
                    (G / depth(x)).sqrt() * e
                } else {
                    0.0
                }
            })
            .collect();
        let saved = dir.join("handoff.json");
        let doc = serde_json::json!({
            "format": "wavesim-handoff", "depth": 5.0, "time": 300.0, "crest": 200.0,
            "x": x, "eta": eta, "ubar": ubar, "broke_seaward_at": broke,
        });
        std::fs::write(&saved, doc.to_string()).unwrap();
        let opts = SliceOptions {
            case: Case::Steep,
            depth: 10.0,
            height: 2.5,
            finest: 0.6,
            frame_every: 4,
            out: dir.join("slice"),
            max_wall_seconds: Some(0.0),
            handoff_depth: None,
            line: Some(line),
            period: 16.0,
            handoff_file: Some(saved),
        };
        run_slice(&opts).unwrap_err().to_string()
    }

    #[test]
    fn a_wave_saved_after_the_breaking_switch_began_on_it_is_refused() {
        // 3 m over 5 m is 0.375 of the total depth; 1 m is 0.17, refused only by the time limit.
        let e = slice_of_saved_wave("broken", 3.0, None);
        assert!(e.contains("part broken"), "{e}");
        let e = slice_of_saved_wave("whole", 1.0, None);
        assert!(e.contains("more than its limit"), "{e}");
    }

    #[test]
    fn a_wave_saved_behind_where_the_run_broke_is_refused() {
        // Low again, as a re-formed wave is, so only where the run broke tells.
        let e = slice_of_saved_wave("reformed", 1.0, Some(150.0));
        assert!(e.contains("already broken once"), "{e}");
    }
}
