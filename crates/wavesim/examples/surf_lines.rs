//! Step 2 of docs/surf-zone.md and its fallback: single cross-shore lines of Pipeline and Coxos,
//! each run without a surf-zone fix and with one, and a table of what the fix changes. The fix is
//! the breaking viscosity, or with `--fix shallow` plain shallow water shoreward of the break
//! line the run without it found.
//!
//!     cargo run --release --example surf_lines                      # every case
//!     cargo run --release --example surf_lines -- coxos_481         # cases whose name has this
//!     cargo run --release --example surf_lines -- --out <dir> ...   # results elsewhere than out/
//!     cargo run --release --example surf_lines -- --cell 1.5 ...    # the Pipeline lines on finer cells
//!     cargo run --release --example surf_lines -- --fix shallow ... # shallow water inside the break
//!
//! Each case is a 2D run's machinery on four rows with walls at the sides, through
//! `wavesim::run` on one grid: the same planning, wave-maker, statistics and break line. Runs go
//! to `<out>/surf_lines/<case>/{plain,viscous,shallow}`, the table to stdout and `summary.json`
//! (`summary_shallow.json` with `--fix shallow`). Run it on one core (`RAYON_NUM_THREADS=1`),
//! under `taskpolicy -b`.

use std::path::{Path, PathBuf};

use serde_json::{Value, json};
use waveio::{Bed, Run};
use wavesim::{BreakLine, RunOptions, run};

/// Rows of a line's bed: a middle row is read, the walls keep them all alike.
const ROWS: usize = 4;
/// The cells the 2D runs use at both spots, and the Coxos bed's own.
const CELL: f64 = 3.0;
/// The row the results are read from.
const ROW: usize = 1;
/// Wall-clock limit per run, so a run that crawls cannot hold the laptop.
const MAX_SECONDS: f64 = 20.0 * 60.0;

/// The bounds fixed in docs/surf-zone.md before the first run.
const MOVE_M: f64 = 10.0;
const DEPTH_SHARE: f64 = 0.10;
const HEIGHT_SHARE: f64 = 0.10;
const SURF_RANGE: (f64, f64) = (0.55, 0.75);
const SURF_MIN_CELLS: usize = 5;

/// The surf-zone fix a case is run with, beside the run without one.
#[derive(Clone, Copy, PartialEq)]
enum Fix {
    Viscosity,
    Shallow,
}

impl Fix {
    fn folder(self) -> &'static str {
        match self {
            Fix::Viscosity => "viscous",
            Fix::Shallow => "shallow",
        }
    }
}

struct Case {
    name: String,
    elevation: Vec<f32>,
    /// Cell size, metres.
    cell: f64,
    height: f64,
    period: f64,
    tide: f64,
}

struct Outcome {
    x: f64,
    depth: f64,
    height: f64,
    surf: Option<(f64, usize)>,
    /// The mean share of the time the wet cells shoreward of the break read as breaking.
    breaking: f64,
    seconds: f64,
    truncated: bool,
}

fn main() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut value = |flag: &str| {
        let k = args.iter().position(|a| a == flag)?;
        let v = args
            .get(k + 1)
            .unwrap_or_else(|| panic!("{flag} needs a value"))
            .clone();
        args.drain(k..k + 2);
        Some(v)
    };
    let out_root = value("--out").map_or(root.join("out"), PathBuf::from);
    let cell: f64 = value("--cell").map_or(CELL, |c| c.parse().expect("--cell in metres"));
    let fix = match value("--fix").as_deref() {
        None | Some("viscosity") => Fix::Viscosity,
        Some("shallow") => Fix::Shallow,
        Some(other) => panic!("--fix is viscosity or shallow, not {other}"),
    };
    // The 3 m runs in surf_lines/, others beside them, so neither overwrites the other.
    let out = if cell == CELL {
        out_root.join("surf_lines")
    } else {
        out_root.join(format!("surf_lines_{cell}m"))
    };
    let data = root.join("data");

    let mut cases = Vec::new();
    for y in ["361.5", "421.5", "511.5"] {
        let line = line_row(&data.join(format!("pipeline_line_{y}.json")), cell);
        for period in [14.0, 16.0] {
            cases.push(Case {
                name: format!("pipeline_{}_2.5m{period}s", y.trim_end_matches(".5")),
                elevation: line.clone(),
                cell,
                height: 2.5,
                period,
                tide: 0.0,
            });
        }
    }
    let coxos = Bed::read(&data.join("coxos/coxos.json")).expect("the Coxos bed in data/");
    let (nx, dy) = (coxos.header.nx, coxos.header.dy);
    assert_eq!(
        (coxos.header.dx, dy),
        (CELL, CELL),
        "the Coxos bed is on 3 m cells"
    );
    // The Coxos bed is on 3 m cells: a finer grid would only interpolate it, so its rows run on
    // the bed's own cells.
    for j in [160_usize, 200, 240].into_iter().filter(|_| cell == CELL) {
        let row = coxos.elevation[j * nx..(j + 1) * nx].to_vec();
        assert!(row.iter().all(|v| v.is_finite()), "Coxos row {j} has gaps");
        for (tide, label) in [(0.0, "mid"), (-1.2, "low")] {
            cases.push(Case {
                name: format!("coxos_{:.0}_2m13s_{label}", (j as f64 + 0.5) * dy),
                elevation: row.clone(),
                cell: CELL,
                height: 2.0,
                period: 13.0,
                tide,
            });
        }
    }
    let estimate_only = args.iter().any(|a| a == "--estimate");
    args.retain(|a| a != "--estimate");
    cases.retain(|c| args.is_empty() || args.iter().any(|a| c.name.contains(a.as_str())));

    // Plan each run and print wavesim's own estimate of its wall-clock time, simulating
    // nothing. The estimate does not count the viscosity's sub-steps.
    if estimate_only {
        let mut total = 0.0;
        for case in &cases {
            let bed = line_bed(case);
            let opts = options(case, Path::new(""), Path::new(""), None);
            let layout = wavesim::layout(&bed, &opts).expect("a run that plans");
            let seconds = wavesim::estimate_seconds(&bed, &opts, &layout);
            total += 2.0 * seconds;
            println!(
                "{}: {} cells, {:.0} s simulated, about {:.1} min per run",
                case.name,
                bed.header.nx,
                layout.coarse.duration,
                seconds / 60.0
            );
        }
        println!(
            "all {} cases, both ways: about {:.0} min",
            cases.len(),
            total / 60.0
        );
        return;
    }

    let mut summary = Vec::new();
    println!(
        "| Case | Break, without | Break, with | Moved | Depth | Height | Inner surf H/h, without → with | Breaking inside, without → with | Verdict | Minutes |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|");
    for case in &cases {
        let dir = out.join(&case.name);
        let bed_path = dir.join("bed.json");
        line_bed(case)
            .write(&bed_path)
            .expect("write the line's bed");
        let plain_dir = dir.join("plain");
        let plain = run_case(case, &bed_path, &plain_dir, None);
        let fixed = plain.as_ref().and_then(|_| {
            let fixed_dir = dir.join(fix.folder());
            run_case(case, &bed_path, &fixed_dir, Some((fix, &plain_dir)))
        });
        let (row, verdict) = compare(&case.name, plain.as_ref(), fixed.as_ref());
        println!("{row}");
        summary.push(json!({ "case": case.name, "verdict": verdict, "plain": outcome_json(plain.as_ref()), fix.folder(): outcome_json(fixed.as_ref()) }));
    }
    std::fs::create_dir_all(&out).expect("the results directory");
    let file = match fix {
        Fix::Viscosity => "summary.json",
        Fix::Shallow => "summary_shallow.json",
    };
    std::fs::write(
        out.join(file),
        serde_json::to_string_pretty(&summary).unwrap() + "\n",
    )
    .expect("write summary.json");
}

/// A `wavesim-line` file put on `cell` cells: the fine data where it has any, the bed's where
/// not, linearly between its points, relative to mean sea level.
fn line_row(path: &Path, cell: f64) -> Vec<f32> {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let doc: Value = serde_json::from_str(&text).expect("a wavesim-line");
    let num = |v: &Value| v.as_f64();
    let x: Vec<f64> = doc["x"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| num(v).unwrap())
        .collect();
    let z_bed = doc["z_bed"].as_array().unwrap();
    let z: Vec<f64> = doc["z"]
        .as_array()
        .unwrap()
        .iter()
        .zip(z_bed)
        .map(|(z, b)| num(z).or_else(|| num(b)).unwrap())
        .collect();
    let nx = (x[x.len() - 1] / cell).floor() as usize;
    (0..nx)
        .map(|i| {
            let at = (i as f64 + 0.5) * cell;
            let k = x.partition_point(|&p| p <= at).clamp(1, x.len() - 1);
            let w = ((at - x[k - 1]) / (x[k] - x[k - 1])).clamp(0.0, 1.0);
            (z[k - 1] + w * (z[k] - z[k - 1])) as f32
        })
        .collect()
}

/// The line as a bed: its row [`ROWS`] times, on the case's cells.
fn line_bed(case: &Case) -> Bed {
    Bed {
        header: Bed::header_for(&case.name, case.elevation.len(), ROWS, case.cell, case.cell),
        elevation: (0..ROWS)
            .flat_map(|_| case.elevation.iter().copied())
            .collect(),
    }
}

/// A case's run: dispersive, one grid, planning its own length, frames every 10 s (the
/// statistics, not the frames, are what is read). A fix comes with the run without it, whose
/// break line [`Fix::Shallow`] takes.
fn options(case: &Case, bed: &Path, out: &Path, fix: Option<(Fix, &Path)>) -> RunOptions {
    let mut opts = RunOptions {
        bed: bed.to_path_buf(),
        out: out.to_path_buf(),
        height: case.height,
        period: case.period,
        tide: case.tide,
        duration: None,
        frame_interval: 10.0,
        manning: 0.0,
        dispersive: true,
        breaking_viscosity: fix.is_some_and(|(f, _)| f == Fix::Viscosity),
        shallow_inside: None,
        nested: false,
        threads: 1,
        background: true,
        budget_seconds: None,
        max_wall_seconds: Some(MAX_SECONDS),
        handoff: Vec::new(),
    };
    if let Some((Fix::Shallow, plain)) = fix {
        let line = BreakLine::from_run(plain, &opts).expect("the plain run's break line");
        opts.shallow_inside = Some(line);
    }
    opts
}

fn run_case(case: &Case, bed: &Path, out: &Path, fix: Option<(Fix, &Path)>) -> Option<Outcome> {
    let opts = options(case, bed, out, fix);
    let summary = match run(&opts, |_| {}) {
        Ok(s) => s,
        Err(e) => {
            eprintln!(
                "{} ({}): {e}",
                case.name,
                fix.map_or("without", |(f, _)| f.folder())
            );
            return None;
        }
    };
    let result = Run::read(out).expect("the run just written");
    let h = &result.header;
    let point = &h.waves["break_line"][ROW];
    if point.is_null() {
        eprintln!("{}: no break point on row {ROW}", case.name);
        return None;
    }
    let (x, depth, height) = (
        point["x_m"].as_f64()?,
        point["depth_m"].as_f64()?,
        point["height_m"].as_f64()?,
    );
    let bed: Vec<f32> = std::fs::read(out.join(&h.bed))
        .expect("the run's bed")
        .as_chunks::<4>()
        .0
        .iter()
        .map(|b| f32::from_le_bytes(*b))
        .collect();
    let wave = result.stat("wave_height")?;
    let start = h.nx * ROW + (x / h.dx - 0.5).round() as usize;
    let mut ratios: Vec<f64> = (start..h.nx * (ROW + 1))
        .map(|k| (case.tide - f64::from(bed[k]), f64::from(wave[k])))
        .take_while(|&(d, _)| d >= 1.0)
        .filter(|&(d, _)| d <= depth)
        .map(|(d, w)| w / d)
        .collect();
    ratios.sort_by(f64::total_cmp);
    let surf = (ratios.len() >= SURF_MIN_CELLS).then(|| (ratios[ratios.len() / 2], ratios.len()));
    let fraction = result.stat("break_fraction")?;
    let wet: Vec<f64> = (start..h.nx * (ROW + 1))
        .filter(|&k| case.tide - f64::from(bed[k]) > 0.05)
        .map(|k| f64::from(fraction[k]))
        .collect();
    let breaking = wet.iter().sum::<f64>() / wet.len().max(1) as f64;
    Some(Outcome {
        x,
        depth,
        height,
        surf,
        breaking,
        seconds: summary.wall_seconds,
        truncated: summary.truncated,
    })
}

/// A table row, and the verdict: "pass", "break moved" or "surf zone out of range".
fn compare(name: &str, plain: Option<&Outcome>, fixed: Option<&Outcome>) -> (String, &'static str) {
    let (Some(p), Some(v)) = (plain, fixed) else {
        return (
            format!("| {name} | failed | | | | | | | failed | |"),
            "failed",
        );
    };
    let moved = v.x - p.x;
    let depth = v.depth / p.depth - 1.0;
    let height = v.height / p.height - 1.0;
    let held = moved.abs() <= MOVE_M && depth.abs() <= DEPTH_SHARE && height.abs() < HEIGHT_SHARE;
    let surf_ok = v
        .surf
        .is_none_or(|(r, _)| (SURF_RANGE.0..=SURF_RANGE.1).contains(&r));
    let verdict = match (held, surf_ok) {
        (false, _) => "break moved",
        (true, false) => "surf zone out of range",
        (true, true) => "pass",
    };
    let surf =
        |s: Option<(f64, usize)>| s.map_or("none".to_string(), |(r, n)| format!("{r:.2} ({n})"));
    let flag = |o: &Outcome| if o.truncated { " (truncated)" } else { "" };
    let row = format!(
        "| {name} | {:.0} m, {:.1} m deep, {:.2} m{} | {:.0} m, {:.1} m deep, {:.2} m{} | {moved:+.0} m | {:+.0}% | {:+.0}% | {} → {} | {:.3} → {:.3} | {verdict} | {:.1} + {:.1} |",
        p.x,
        p.depth,
        p.height,
        flag(p),
        v.x,
        v.depth,
        v.height,
        flag(v),
        100.0 * depth,
        100.0 * height,
        surf(p.surf),
        surf(v.surf),
        p.breaking,
        v.breaking,
        p.seconds / 60.0,
        v.seconds / 60.0,
    );
    (row, verdict)
}

fn outcome_json(o: Option<&Outcome>) -> Value {
    o.map_or(Value::Null, |o| {
        json!({
            "x_m": o.x, "depth_m": o.depth, "height_m": o.height,
            "inner_surf_h_over_d": o.surf.map(|s| s.0), "inner_surf_cells": o.surf.map(|s| s.1),
            "breaking_inside": o.breaking,
            "wall_seconds": o.seconds, "truncated": o.truncated,
        })
    })
}
