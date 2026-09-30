//! A whole run on a small synthetic bed: bed file in, frames out, waves in between.

use std::fs;
use std::path::PathBuf;

use waveio::{Bed, Run};
use wavesim::{Error, RunOptions, plan, run};

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
    assert!(out.eta.iter().all(|v| v.is_finite()));

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
    assert_eq!(out.header.fields, vec!["eta"]);
    assert!(opts.out.join("bed.f32").exists());
}
