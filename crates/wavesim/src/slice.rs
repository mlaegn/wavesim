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

/// Run one wave until its lip lands, it reaches the shallow wall, the time limit passes, or the
/// computation fails, writing frames as it goes. A run its first steps say would take longer
/// than the limit is refused before anything is written.
pub fn run_slice(opts: &SliceOptions) -> Result<SliceSummary, Error> {
    let profile = opts.case.profile(opts.depth);
    if opts.height <= 0.0 || opts.height > 0.6 * opts.depth {
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
    let mut tank = Tank::over(&wave, &profile, opts.finest);
    let speeds = |t: &Tank| -> Vec<f64> { t.velocity().iter().map(|v| v.0.hypot(v.1)).collect() };
    let start = (tank.surface.clone(), speeds(&tank));

    // Time the first steps and refuse a run that would go over its limit.
    let clock = Instant::now();
    for _ in 0..CALIBRATION_STEPS {
        tank.march();
    }
    let per_step = clock.elapsed().as_secs_f64() / CALIBRATION_STEPS as f64;
    let mean_dt = tank.time / CALIBRATION_STEPS as f64;
    let steps = travel_time(&profile, crest, opts.height) / mean_dt * SLOWDOWN;
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

    let settings = serde_json::json!({
        "case": opts.case.name(),
        "depth": opts.depth,
        "height": opts.height,
        "crest": crest,
        "g": G,
        "finest": opts.finest,
        "profile": profile,
        "wave": "solitary",
    });
    let mut out = SliceWriter::create(&opts.out, &tank.bed, tank.surface.len(), settings)?;
    out.write_frame(0.0, &start.0, &start.1)?;
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
