//! A whole run on a small synthetic bed: bed file in, frames out, waves in between.

use std::fs;
use std::path::PathBuf;

use wavecore::linear_wave;
use waveio::{Bed, Run};
use wavesim::{
    Error, MAX_SECOND_HARMONIC, RunOptions, estimate_seconds, plan, run, second_harmonic,
};

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("wavesim-test-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

const NX: usize = 400;
const NY: usize = 8;
const DX: f64 = 2.0;

/// Still-water depth of the test beach: 16 m at x = 0, rising 1 m per 45 m, dry from 720 m.
fn beach_depth(x: f64) -> f64 {
    16.0 - x / 45.0
}

/// The test beach, uniform alongshore.
fn beach(dir: &std::path::Path) -> PathBuf {
    let elevation = (0..NX * NY)
        .map(|n| -beach_depth(((n % NX) as f64 + 0.5) * DX) as f32)
        .collect();
    let bed = Bed {
        header: Bed::header_for("beach", NX, NY, DX, DX),
        elevation,
    };
    let path = dir.join("beach.json");
    bed.write(&path).unwrap();
    path
}

/// The bed a run actually used: the beach with the columns it did not need cropped off.
fn run_bed(out: &std::path::Path) -> Vec<f32> {
    fs::read(out.join("bed.f32"))
        .unwrap()
        .as_chunks::<4>()
        .0
        .iter()
        .map(|b| f32::from_le_bytes(*b))
        .collect()
}

fn options(dir: &std::path::Path, bed: PathBuf) -> RunOptions {
    RunOptions {
        bed,
        out: dir.join("run"),
        height: 0.2,
        period: 8.0,
        tide: 0.0,
        duration: Some(60.0),
        frame_interval: 5.0,
        manning: 0.0,
        dispersive: false,
        max_wall_seconds: None,
    }
}

#[test]
fn the_second_harmonic_has_the_stokes_limits() {
    let a = 0.5;
    // Deep water: k a / 2.
    let k = 0.1;
    assert!((second_harmonic(a, k, 200.0) / (k * a / 2.0) - 1.0).abs() < 1e-9);
    // Shallow water: 3 a / (4 k^2 h^3).
    let (k, h) = (0.01, 2.0);
    let shallow = 3.0 * a / (4.0 * k * k * h * h * h);
    assert!((second_harmonic(a, k, h) / shallow - 1.0).abs() < 1e-3);
}

#[test]
fn the_wave_maker_goes_as_shallow_as_a_sinusoid_is_still_the_right_wave() {
    let dir = scratch("plan");
    let opts = options(&dir, beach(&dir));
    let bed = Bed::read(&opts.bed).unwrap();
    let p = plan(&bed, &opts).unwrap();
    let omega = std::f64::consts::TAU / opts.period;
    let harmonic = |x: f64| {
        let h = beach_depth(x);
        second_harmonic(opts.height / 2.0, linear_wave(omega, h, false).k, h)
    };
    // Where it stands the sinusoid lacks at most a tenth of the wave, and one cell shoreward
    // it would lack more.
    assert!((p.maker_depth - beach_depth(p.maker_x)).abs() < 1e-4);
    assert!(p.second_harmonic <= MAX_SECOND_HARMONIC && harmonic(p.maker_x) <= MAX_SECOND_HARMONIC);
    assert!(
        harmonic(p.maker_x + DX) > MAX_SECOND_HARMONIC,
        "maker at {}",
        p.maker_x
    );
    // The run keeps just enough of the bed offshore for the sponge behind the maker.
    let kept = p.maker_x - p.crop as f64 * DX;
    assert!(kept >= p.sponge_offshore as f64 * DX && kept < p.sponge_offshore as f64 * DX + 30.0);
}

#[test]
fn a_higher_tide_lets_the_wave_maker_sit_further_shoreward() {
    let dir = scratch("tide");
    let mut opts = options(&dir, beach(&dir));
    let bed = Bed::read(&opts.bed).unwrap();
    let low = plan(&bed, &opts).unwrap();
    opts.tide = 1.0;
    let high = plan(&bed, &opts).unwrap();
    assert!(
        high.maker_x > low.maker_x,
        "{} vs {}",
        high.maker_x,
        low.maker_x
    );
}

#[test]
fn a_planned_run_lasts_until_the_swell_has_reached_the_shore_and_broken_a_few_times() {
    let dir = scratch("duration");
    let mut opts = options(&dir, beach(&dir));
    opts.duration = None;
    let bed = Bed::read(&opts.bed).unwrap();
    let p = plan(&bed, &opts).unwrap();
    // Long-wave travel time from the maker to the shore at sqrt(g h) on a plane beach of slope
    // s: the integral of dx / sqrt(g s x) is 2 sqrt(L / (g s)) over a distance L. The plan
    // counts the last stretch, shallower than half a metre, at the speed in half a metre.
    let gs = wavecore::G / 45.0;
    let (distance, last) = (720.0 - p.maker_x, 0.5 * 45.0);
    let travel =
        2.0 * ((distance / gs).sqrt() - (last / gs).sqrt()) + last / (0.5 * wavecore::G).sqrt();
    let expected = (2.0 + wavesim::BREAKS_SHOWN) * opts.period + travel;
    assert!(
        (p.duration / expected - 1.0).abs() < 0.05,
        "{} vs {expected}",
        p.duration
    );
}

#[test]
fn impossible_setups_are_explained() {
    let dir = scratch("errors");
    let mut opts = options(&dir, beach(&dir));
    let bed = Bed::read(&opts.bed).unwrap();
    // A 5 m, 14 s swell is far from a sinusoid in 16 m of water.
    opts.height = 5.0;
    opts.period = 14.0;
    let err = plan(&bed, &opts).unwrap_err();
    assert!(matches!(err, Error::Setup(_)));
    assert!(err.to_string().contains("offshore_m"), "{err}");

    opts.height = -1.0;
    assert!(
        plan(&bed, &opts)
            .unwrap_err()
            .to_string()
            .contains("height")
    );
}

#[test]
fn a_run_writes_frames_that_contain_the_requested_wave() {
    let dir = scratch("run");
    let opts = options(&dir, beach(&dir));
    let summary = run(&opts, |_| {}).unwrap();

    let out = Run::read(&opts.out).unwrap();
    assert_eq!(summary.frames, 13); // t = 0, 5, ..., 60
    assert_eq!(out.header.frame_count, 13);
    assert_eq!(out.header.times.first(), Some(&0.0));
    assert!((out.header.times[12] - 60.0).abs() < 1e-6);
    assert!(out.data.iter().all(|v| v.is_finite()));

    // The first frame is still water: level surface where wet, the bed where dry.
    let bed = run_bed(&opts.out);
    assert_eq!(bed.len(), out.header.nx * NY);
    for (n, (&eta, &z)) in out.frame(0).iter().zip(&bed).enumerate() {
        let expected = z.max(0.0);
        assert!(
            (eta - expected).abs() < 1e-5,
            "cell {n}: {eta} vs {expected}"
        );
    }

    // Later, a wave of the requested amplitude (0.1 m) passes between the wave-maker's own
    // bump and 2.5 m of water, before shoaling has grown it much.
    let nx = out.header.nx;
    let start = out.header.waves["near_field_end_m"].as_f64().unwrap();
    let amplitude = 0.1;
    let mut peak: f32 = 0.0;
    for k in 6..13 {
        for (i, &z) in bed[..nx].iter().enumerate() {
            let x = (i as f64 + 0.5) * DX;
            if x > start && z < -2.5 {
                peak = peak.max(out.frame(k)[NY / 2 * nx + i].abs());
            }
        }
    }
    assert!(
        (0.6 * amplitude..=1.3 * amplitude).contains(&f64::from(peak)),
        "peak {peak} m is not within reach of the requested amplitude {amplitude} m"
    );
    assert!(summary.max_rise > 0.0);
}

#[test]
fn an_alongshore_uniform_beach_gives_an_alongshore_uniform_run() {
    // The sides are walls, and for a swell travelling along +x a wall is a mirror: on a beach
    // that does not change alongshore, every row must see exactly the same wave. Absorbing
    // strips at the sides, which this replaced, damp the swell near them and make it diffract.
    let dir = scratch("uniform");
    let opts = options(&dir, beach(&dir));
    run(&opts, |_| {}).unwrap();
    let out = Run::read(&opts.out).unwrap();
    let nx = out.header.nx;
    for k in 0..out.header.frame_count {
        let eta = out.frame(k);
        for j in 1..NY {
            for i in 0..nx {
                let (a, b) = (eta[i], eta[j * nx + i]);
                assert!(
                    (a - b).abs() < 1e-6,
                    "frame {k}, row {j}, column {i}: {b} vs {a}"
                );
            }
        }
    }
}

#[test]
fn the_run_directory_is_self_describing() {
    let dir = scratch("describe");
    let opts = options(&dir, beach(&dir));
    run(&opts, |_| {}).unwrap();
    let out = Run::read(&opts.out).unwrap();
    assert_eq!(out.header.waves["wave_period_s"], 8.0);
    assert_eq!(out.header.waves["wave_height_m"], 0.2);
    assert!(out.header.waves["maker_second_harmonic"].as_f64().unwrap() <= MAX_SECOND_HARMONIC);
    assert!(out.header.waves["cropped_offshore_m"].as_f64().unwrap() > 0.0);
    assert_eq!(out.header.fields, vec!["eta", "breaking"]);
    assert!(opts.out.join("bed.f32").exists());
}

#[test]
fn a_dispersive_run_writes_the_same_kind_of_frames_and_reports_its_solver_work() {
    let dir = scratch("dispersive");
    let mut opts = options(&dir, beach(&dir));
    opts.dispersive = true;
    let summary = run(&opts, |_| {}).unwrap();

    let out = Run::read(&opts.out).unwrap();
    assert_eq!(out.header.frame_count, 13);
    assert_eq!(out.header.waves["dispersive"], true);
    assert!(
        out.header.waves["model"]
            .as_str()
            .unwrap()
            .contains("Serre")
    );
    assert!(out.data.iter().all(|v| v.is_finite()));

    let stats = summary.dispersion.expect("dispersive run reports stats");
    assert!(stats.solves > 0 && stats.unconverged == 0, "{stats:?}");
    assert!(summary.max_rise > 0.0);
}

#[test]
fn a_wave_too_short_for_the_dispersive_model_is_explained() {
    let dir = scratch("short");
    let mut opts = options(&dir, beach(&dir));
    opts.dispersive = true;
    // omega^2 h / g = 64 in 16 m of water for a 1 s wave: far beyond what the dispersive model
    // can carry.
    opts.period = 1.0;
    let bed = Bed::read(&opts.bed).unwrap();
    let err = plan(&bed, &opts).unwrap_err();
    assert!(err.to_string().contains("too short"), "{err}");
}

#[test]
fn the_breaking_field_marks_steep_tall_waves_near_the_shore_and_nothing_else() {
    let dir = scratch("breaking");
    let mut opts = options(&dir, beach(&dir));
    opts.dispersive = true;
    opts.height = 1.5;
    opts.duration = None;
    run(&opts, |_| {}).unwrap();

    let out = Run::read(&opts.out).unwrap();
    assert!(out.header.fields.contains(&"breaking".to_string()));
    let bed = run_bed(&opts.out);

    // Still water at the start: nothing is breaking.
    assert!(out.field("breaking", 0).unwrap().iter().all(|&b| b == 0.0));

    let mut strongest: f32 = 0.0;
    for k in 0..out.header.frame_count {
        for (n, &b) in out.field("breaking", k).unwrap().iter().enumerate() {
            assert!((0.0..=1.0).contains(&b), "frame {k}, cell {n}: {b}");
            strongest = strongest.max(b);
            // The value is how far the dispersive terms have faded, so a wave 40% as tall as the
            // water is deep already reads about 0.5 (the fade starts at eta/h = 0.30) without
            // breaking. A clear reading, above 0.8, must be in shallow water: a 1.5 m wave is
            // nowhere near breaking in 3.5 m or more.
            if b > 0.8 {
                let depth = -f64::from(bed[n]);
                assert!(
                    depth < 3.5,
                    "frame {k}, cell {n}: breaking in {depth} m of water"
                );
            }
        }
    }
    assert!(
        strongest > 0.9,
        "the wave never read as breaking (strongest {strongest})"
    );
}

#[test]
fn a_run_that_hits_its_time_budget_stops_and_keeps_a_valid_shorter_run() {
    let dir = scratch("budget");
    let mut opts = options(&dir, beach(&dir));
    opts.duration = Some(600.0); // far more than the budget allows
    opts.max_wall_seconds = Some(0.3);
    let summary = run(&opts, |_| {}).unwrap();
    assert!(summary.truncated);
    assert!(
        summary.wall_seconds < 3.0,
        "ran for {} s",
        summary.wall_seconds
    );

    // What was written is a complete run in its own right.
    let out = Run::read(&opts.out).unwrap();
    assert!(out.header.frame_count >= 1 && out.header.frame_count < 120);
    assert_eq!(out.header.waves["truncated"], true);
    let stopped = out.header.waves["stopped_at_s"].as_f64().unwrap();
    assert!(stopped < 600.0 && stopped >= *out.header.times.last().unwrap());
    assert!(out.data.iter().all(|v| v.is_finite()));
}

#[test]
fn a_run_inside_its_budget_is_not_marked_truncated() {
    let dir = scratch("inbudget");
    let mut opts = options(&dir, beach(&dir));
    opts.max_wall_seconds = Some(600.0);
    let summary = run(&opts, |_| {}).unwrap();
    assert!(!summary.truncated);
    let out = Run::read(&opts.out).unwrap();
    assert!(out.header.waves.get("truncated").is_none());
}

#[test]
fn the_time_estimate_grows_with_the_work_and_shrinks_with_threads() {
    let dir = scratch("estimate");
    let mut opts = options(&dir, beach(&dir));
    let bed = Bed::read(&opts.bed).unwrap();
    let estimate = |opts: &RunOptions, threads| {
        estimate_seconds(&bed, opts, &plan(&bed, opts).unwrap(), threads)
    };
    let base = estimate(&opts, 4);
    assert!(base > 0.0 && base.is_finite());
    opts.duration = Some(120.0);
    assert!((estimate(&opts, 4) / base - 2.0).abs() < 1e-9);
    opts.duration = Some(60.0);
    assert!(estimate(&opts, 1) > estimate(&opts, 4));
    opts.dispersive = true;
    assert!(estimate(&opts, 4) > base, "dispersion costs more");
}
