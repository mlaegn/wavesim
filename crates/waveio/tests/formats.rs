//! Bed and run files must round-trip exactly and reject malformed input loudly.

use std::fs;
use std::path::PathBuf;

use waveio::{Bed, Error, Run, RunWriter};

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("waveio-test-{name}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

fn sample_bed(nx: usize, ny: usize) -> Bed {
    let elevation = (0..nx * ny).map(|n| -5.0 + 0.01 * n as f32).collect();
    Bed {
        header: Bed::header_for("sample", nx, ny, 3.0, 2.0),
        elevation,
    }
}

#[test]
fn bed_round_trips_exactly() {
    let dir = scratch("bed");
    let bed = sample_bed(7, 5);
    bed.write(&dir.join("sample.json")).unwrap();

    let back = Bed::read(&dir.join("sample.json")).unwrap();
    assert_eq!(back.header, bed.header);
    assert_eq!(back.elevation, bed.elevation);
    let grid = back.grid();
    assert_eq!((grid.nx, grid.ny, grid.dx, grid.dy), (7, 5, 3.0, 2.0));
}

#[test]
fn bathymetry_moves_still_water_to_zero() {
    let dir = scratch("tide");
    let bed = sample_bed(4, 3);
    bed.write(&dir.join("sample.json")).unwrap();
    let grid = bed.grid();
    let b = bed.bathymetry(1.5);
    let k = grid.idx(wavecore::GHOST, wavecore::GHOST);
    assert!((b.b[k] - (-5.0 - 1.5)).abs() < 1e-6);
}

#[test]
fn wrong_length_is_rejected() {
    let dir = scratch("short");
    let bed = sample_bed(4, 3);
    bed.write(&dir.join("sample.json")).unwrap();
    let data = dir.join("sample.f32");
    let mut bytes = fs::read(&data).unwrap();
    bytes.truncate(bytes.len() - 4);
    fs::write(&data, bytes).unwrap();

    let err = Bed::read(&dir.join("sample.json")).unwrap_err();
    assert!(matches!(err, Error::Format { .. }), "{err}");
    assert!(err.to_string().contains("bytes"), "{err}");
}

#[test]
fn gaps_in_the_bed_are_rejected() {
    let dir = scratch("nan");
    let mut bed = sample_bed(4, 3);
    bed.elevation[5] = f32::NAN;
    bed.write(&dir.join("sample.json")).unwrap();
    let err = Bed::read(&dir.join("sample.json")).unwrap_err();
    assert!(err.to_string().contains("not finite"), "{err}");
}

#[test]
fn a_different_format_is_rejected() {
    let dir = scratch("format");
    let mut bed = sample_bed(4, 3);
    bed.header.format = "something-else".into();
    bed.write(&dir.join("sample.json")).unwrap();
    let err = Bed::read(&dir.join("sample.json")).unwrap_err();
    assert!(err.to_string().contains("wavesim-bed"), "{err}");
}

#[test]
fn run_round_trips_frames_and_times() {
    let dir = scratch("run");
    let bed = sample_bed(6, 4);
    let mut w =
        RunWriter::create(&dir, &bed, &["eta"], serde_json::json!({"period": 12.0})).unwrap();
    for k in 0..3 {
        let eta: Vec<f32> = (0..24).map(|n| k as f32 + n as f32 * 0.5).collect();
        w.write_frame(2.0 * k as f64, &[&eta]).unwrap();
    }
    w.finish().unwrap();

    let run = Run::read(&dir).unwrap();
    assert_eq!(run.header.frame_count, 3);
    assert_eq!(run.header.times, vec![0.0, 2.0, 4.0]);
    assert_eq!(run.header.fields, vec!["eta"]);
    assert_eq!(run.header.waves["period"], 12.0);
    assert_eq!(run.frame(2)[3], 2.0 + 1.5);
    assert_eq!(fs::read(dir.join("bed.f32")).unwrap().len(), 6 * 4 * 4);
}

#[test]
fn a_frame_of_the_wrong_size_is_rejected() {
    let dir = scratch("badframe");
    let bed = sample_bed(6, 4);
    let mut w = RunWriter::create(&dir, &bed, &["eta"], serde_json::Value::Null).unwrap();
    assert!(w.write_frame(0.0, &[&[0.0; 23]]).is_err());
}

#[test]
fn a_run_can_carry_several_fields_per_frame() {
    let dir = scratch("fields");
    let bed = sample_bed(6, 4);
    let mut w =
        RunWriter::create(&dir, &bed, &["eta", "breaking"], serde_json::Value::Null).unwrap();
    for k in 0..3 {
        let eta: Vec<f32> = (0..24).map(|n| 10.0 * k as f32 + n as f32).collect();
        let breaking: Vec<f32> = (0..24).map(|n| 0.5 + 100.0 * k as f32 + n as f32).collect();
        w.write_frame(k as f64, &[&eta, &breaking]).unwrap();
    }
    w.finish().unwrap();

    let run = Run::read(&dir).unwrap();
    assert_eq!(run.header.fields, vec!["eta", "breaking"]);
    assert_eq!(run.frame(2)[5], 25.0);
    assert_eq!(run.field("breaking", 2).unwrap()[5], 205.5);
    assert_eq!(run.field("breaking", 0).unwrap()[0], 0.5);
    assert!(run.field("velocity", 0).is_none());
    assert_eq!(
        fs::read(dir.join("frames.f32")).unwrap().len(),
        3 * 2 * 24 * 4
    );
}

#[test]
fn the_wrong_number_of_fields_or_a_missing_eta_is_rejected() {
    let dir = scratch("badfields");
    let bed = sample_bed(6, 4);
    assert!(RunWriter::create(&dir, &bed, &["breaking"], serde_json::Value::Null).is_err());
    let mut w =
        RunWriter::create(&dir, &bed, &["eta", "breaking"], serde_json::Value::Null).unwrap();
    assert!(
        w.write_frame(0.0, &[&[0.0; 24]]).is_err(),
        "one field where two are expected"
    );
    assert!(
        w.write_frame(0.0, &[&[0.0; 24], &[0.0; 23]]).is_err(),
        "a short field"
    );
}

#[test]
fn a_note_is_recorded_in_the_header() {
    let dir = scratch("note");
    let bed = sample_bed(6, 4);
    let mut w = RunWriter::create(&dir, &bed, &["eta"], serde_json::Value::Null).unwrap();
    w.write_frame(0.0, &[&[0.0; 24]]).unwrap();
    w.note("truncated", serde_json::json!(true));
    w.finish().unwrap();
    assert_eq!(Run::read(&dir).unwrap().header.waves["truncated"], true);
}
