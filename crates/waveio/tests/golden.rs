//! The viewer's tests read the files this test writes, so the two sides cannot drift.
//!
//! `viewer/tests/fixtures/tiny/` is a committed run. This test writes the same run to a
//! scratch directory and compares every byte. If the writer changes on purpose,
//! regenerate the fixture and commit it:
//!
//!     UPDATE_GOLDEN=1 cargo test -p waveio --test golden

use std::fs;
use std::path::{Path, PathBuf};

use waveio::{Bed, Frame, RunWriter};

fn fixture_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../viewer/tests/fixtures/tiny")
}

/// A 6 x 4 bed rising 0.75 m per cell from -2 m, so the last three columns are dry,
/// with four frames. Every value is a multiple of 1/16, exact in f32 on any platform.
fn write_tiny(dir: &Path) {
    let (nx, ny) = (6, 4);
    let elevation: Vec<f32> = (0..nx * ny)
        .map(|n| -2.0 + 0.75 * (n % nx) as f32)
        .collect();
    let mut header = Bed::header_for("tiny", nx, ny, 3.0, 2.0);
    header.description = "Synthetic fixture for the viewer's tests".into();
    header.frame = Frame {
        crs: "EPSG:32604".into(),
        origin_easting: 1000.0,
        origin_northing: 2000.0,
        x_bearing_deg: 130.0,
        note: "fixture".into(),
    };
    header.source = serde_json::json!({ "name": "Synthetic fixture", "attribution": "none" });
    let bed = Bed { header, elevation };

    let waves = serde_json::json!({
        "wave_height_m": 1.0,
        "wave_period_s": 14.0,
        "tide_m": 0.0,
        "maker_x_m": 6.0,
        "sponge_offshore_cells": 2,
        "sponge_side_cells": 1,
    });
    let mut writer = RunWriter::create(dir, &bed, waves).unwrap();
    for k in 0..4 {
        let eta: Vec<f32> = bed
            .elevation
            .iter()
            .enumerate()
            .map(|(n, &z)| {
                if z < 0.0 {
                    let wave = 0.125 * k as f32 + 0.0625 * (n % nx) as f32 - 0.25;
                    wave.max(z)
                } else {
                    z
                }
            })
            .collect();
        writer.write_frame(2.0 * k as f64, &eta).unwrap();
    }
    writer.finish().unwrap();
}

#[test]
fn the_run_the_viewer_tests_read_is_what_the_writer_produces() {
    let fixture = fixture_dir();
    if std::env::var_os("UPDATE_GOLDEN").is_some() {
        let _ = fs::remove_dir_all(&fixture);
        write_tiny(&fixture);
        return;
    }

    let scratch = std::env::temp_dir().join(format!("waveio-golden-{}", std::process::id()));
    let _ = fs::remove_dir_all(&scratch);
    write_tiny(&scratch);
    for name in ["run.json", "bed.f32", "frames.f32"] {
        let want = fs::read(fixture.join(name))
            .unwrap_or_else(|e| panic!("{name}: {e}; run with UPDATE_GOLDEN=1 to create it"));
        let got = fs::read(scratch.join(name)).unwrap();
        assert_eq!(
            got, want,
            "{name} differs from viewer/tests/fixtures/tiny; if the format change is intended, \
             regenerate with UPDATE_GOLDEN=1 cargo test -p waveio --test golden"
        );
    }
}
