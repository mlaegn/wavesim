//! A whole run on a small synthetic bed: bed file in, frames out, waves in between.

use std::fs;
use std::path::PathBuf;

use waveio::{Bed, Run};
use wavesim::{Error, RunOptions, estimate_seconds, plan, run};

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("wavesim-test-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

/// A beach that rises 1 m per 100 m, uniform alongshore: 6 m deep at x = 0, dry at 600 m.
fn beach(dir: &std::path::Path) -> PathBuf {
    let (nx, ny, dx) = (300, 40, 2.0);
    let elevation = (0..nx * ny)
        .map(|n| {
            let x = ((n % nx) as f64 + 0.5) * dx;
            (-6.0 + x / 100.0) as f32
        })
        .collect();
    let bed = Bed {
        header: Bed::header_for("beach", nx, ny, dx, dx),
        elevation,
    };
    let path = dir.join("beach.json");
    bed.write(&path).unwrap();
    path
}

fn options(dir: &std::path::Path, bed: PathBuf) -> RunOptions {
    RunOptions {
        bed,
        out: dir.join("run"),
        height: 0.02,
        period: 8.0,
        tide: 0.0,
        duration: 60.0,
        frame_interval: 5.0,
        maker_depth: 4.0,
        manning: 0.0,
        dispersive: false,
        max_wall_seconds: None,
    }
}

#[test]
fn the_wave_maker_goes_where_the_water_first_gets_shallow_enough() {
    let dir = scratch("plan");
    let opts = options(&dir, beach(&dir));
    let bed = Bed::read(&opts.bed).unwrap();
    let p = plan(&bed, &opts).unwrap();
    // Depth 4 m is reached at x = 200 m; the maker sits in the cell that first satisfies it.
    assert!((p.maker_x - 201.0).abs() <= 2.0, "maker at {}", p.maker_x);
    assert!(
        (p.maker_depth - 4.0).abs() < 0.05,
        "depth {}",
        p.maker_depth
    );
    assert!(p.sponge_offshore >= 4 && p.sponge_side >= 4);
}

#[test]
fn a_tide_moves_the_maker_seaward_because_the_water_is_deeper() {
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
fn impossible_setups_are_explained() {
    let dir = scratch("errors");
    let mut opts = options(&dir, beach(&dir));
    opts.maker_depth = 40.0;
    let bed = Bed::read(&opts.bed).unwrap();
    let err = plan(&bed, &opts).unwrap_err();
    assert!(matches!(err, Error::Setup(_)));
    assert!(err.to_string().contains("--maker-depth"), "{err}");

    opts.maker_depth = 0.01;
    let err = plan(&bed, &opts).unwrap_err();
    assert!(err.to_string().contains("--maker-depth"), "{err}");

    opts.maker_depth = 4.0;
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
    let bed = Bed::read(&opts.bed).unwrap();
    for (n, (&eta, &z)) in out.frame(0).iter().zip(&bed.elevation).enumerate() {
        let expected = z.max(0.0);
        assert!(
            (eta - expected).abs() < 1e-5,
            "cell {n}: {eta} vs {expected}"
        );
    }

    // Later, a wave of the requested amplitude (0.01 m) is running up the beach.
    let (nx, ny) = (300, 40);
    let amplitude = 0.01;
    let mut peak: f32 = 0.0;
    for k in 6..13 {
        for j in ny / 2 - 3..ny / 2 + 3 {
            for i in (250 / 2)..(400 / 2) {
                peak = peak.max(out.frame(k)[j * nx + i].abs());
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
fn the_run_directory_is_self_describing() {
    let dir = scratch("describe");
    let opts = options(&dir, beach(&dir));
    run(&opts, |_| {}).unwrap();
    let out = Run::read(&opts.out).unwrap();
    assert_eq!(out.header.waves["wave_period_s"], 8.0);
    assert_eq!(out.header.waves["wave_height_m"], 0.02);
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
    // omega^2 h / g = 20 at 5 m for a 1 s wave: far beyond what the dispersive model can carry.
    opts.period = 1.0;
    opts.maker_depth = 5.0;
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
    opts.duration = 130.0;
    run(&opts, |_| {}).unwrap();

    let out = Run::read(&opts.out).unwrap();
    assert!(out.header.fields.contains(&"breaking".to_string()));
    let bed = Bed::read(&opts.bed).unwrap();

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
                let depth = -f64::from(bed.elevation[n]);
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
    opts.duration = 600.0; // far more than the budget allows
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
    let base = estimate_seconds(&bed, &opts, 4);
    assert!(base > 0.0 && base.is_finite());
    opts.duration *= 2.0;
    assert!((estimate_seconds(&bed, &opts, 4) / base - 2.0).abs() < 1e-9);
    opts.duration /= 2.0;
    assert!(estimate_seconds(&bed, &opts, 1) > estimate_seconds(&bed, &opts, 4));
    opts.dispersive = true;
    assert!(
        estimate_seconds(&bed, &opts, 4) > base,
        "dispersion costs more"
    );
}
