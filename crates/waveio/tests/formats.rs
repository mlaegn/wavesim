//! Bed, run and slice files must round-trip exactly and reject malformed input loudly.

use std::fs;
use std::path::PathBuf;

use waveio::{Bed, Error, Landing, Run, RunWriter, SLICE_FORMAT, Slice, SliceWriter};

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
fn cropping_drops_offshore_columns_and_moves_the_origin_with_them() {
    let mut bed = sample_bed(7, 5);
    bed.header.frame.origin_easting = 1000.0;
    bed.header.frame.origin_northing = 2000.0;
    // +x points east: dropping 2 columns of 3 m moves the corner 6 m east.
    bed.header.frame.x_bearing_deg = 90.0;
    let cropped = bed.crop_x(2);
    assert_eq!((cropped.header.nx, cropped.header.ny), (5, 5));
    for j in 0..5 {
        assert_eq!(
            cropped.elevation[j * 5..j * 5 + 5],
            bed.elevation[j * 7 + 2..j * 7 + 7]
        );
    }
    assert!((cropped.header.frame.origin_easting - 1006.0).abs() < 1e-9);
    assert!((cropped.header.frame.origin_northing - 2000.0).abs() < 1e-9);
}

#[test]
fn cropping_alongshore_moves_the_origin_ninety_degrees_counter_clockwise_from_x() {
    let mut bed = sample_bed(7, 5);
    bed.header.frame.x_bearing_deg = 90.0; // +x east, so +y north
    let cropped = bed.crop(1..4, 2..5);
    assert_eq!((cropped.header.nx, cropped.header.ny), (3, 3));
    assert_eq!(cropped.elevation[0], bed.elevation[2 * 7 + 1]);
    assert!((cropped.header.frame.origin_easting - 3.0).abs() < 1e-9);
    assert!((cropped.header.frame.origin_northing - 4.0).abs() < 1e-9);
}

#[test]
fn coarsening_averages_blocks_and_drops_what_is_left_over() {
    let bed = sample_bed(7, 5);
    let coarse = bed.coarsen(2);
    assert_eq!((coarse.header.nx, coarse.header.ny), (3, 2));
    assert_eq!((coarse.header.dx, coarse.header.dy), (6.0, 4.0));
    let e = |i: usize, j: usize| f64::from(bed.elevation[j * 7 + i]);
    let expected = (e(2, 2) + e(3, 2) + e(2, 3) + e(3, 3)) / 4.0;
    assert!((f64::from(coarse.elevation[3 + 1]) - expected).abs() < 1e-5);
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

#[test]
fn statistics_round_trip_and_are_found_by_name() {
    let dir = scratch("stats");
    let bed = sample_bed(3, 2);
    let mut w = RunWriter::create(&dir, &bed, &["eta"], serde_json::json!({})).unwrap();
    w.write_frame(0.0, &[&[0.0; 6]]).unwrap();
    let fraction = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5];
    let high = [1.0; 6];
    w.write_stats(
        &["break_fraction", "eta_max"],
        &[&fraction, &high],
        10.0,
        20.0,
    )
    .unwrap();
    w.finish().unwrap();
    let back = Run::read(&dir).unwrap();
    assert_eq!(back.stat("break_fraction"), Some(&fraction[..]));
    assert_eq!(back.stat("eta_max"), Some(&high[..]));
    assert_eq!(back.stat("eta_min"), None);
    let st = back.header.stats.unwrap();
    assert_eq!((st.from_s, st.to_s), (10.0, 20.0));
}

#[test]
fn a_run_without_statistics_has_none() {
    let dir = scratch("nostats");
    let bed = sample_bed(3, 2);
    let mut w = RunWriter::create(&dir, &bed, &["eta"], serde_json::json!({})).unwrap();
    w.write_frame(0.0, &[&[0.0; 6]]).unwrap();
    w.finish().unwrap();
    let back = Run::read(&dir).unwrap();
    assert!(back.header.stats.is_none() && back.stats.is_none());
    let text = std::fs::read_to_string(dir.join("run.json")).unwrap();
    assert!(
        !text.contains("stats"),
        "an absent field should not be written"
    );
}

#[test]
fn a_slice_round_trips_and_says_how_it_ended() {
    let dir = scratch("slice");
    let bed = [(0.0, -1.0), (5.0, -1.0), (10.0, -0.2)];
    let mut w = SliceWriter::create(&dir, &bed, 3, serde_json::json!({ "height": 0.3 })).unwrap();
    w.write_frame(
        0.0,
        &[(0.0, 0.0), (5.0, 0.3), (10.0, 0.0)],
        &[0.0, 1.5, 0.0],
    )
    .unwrap();
    w.write_frame(
        0.25,
        &[(0.0, 0.0), (5.5, 0.32), (10.0, 0.01)],
        &[0.0, 1.75, 0.1],
    )
    .unwrap();
    assert!(
        w.write_frame(0.5, &[(0.0, 0.0)], &[0.0]).is_err(),
        "a frame of the wrong size"
    );
    let landing = Landing {
        time: 0.25,
        x: 6.0,
        z: 0.0,
        throw: 0.5,
        tube_area: 0.1,
        tube_width: 0.4,
        tube_height: 0.3,
    };
    w.finish("landed", None, Some(landing)).unwrap();

    let s = Slice::read(&dir).unwrap();
    assert_eq!(s.header.format, SLICE_FORMAT);
    assert_eq!((s.header.frame_count, s.header.nodes), (2, 3));
    assert_eq!(s.header.times, vec![0.0, 0.25]);
    assert_eq!(s.header.ended, "landed");
    assert_eq!(s.header.landing, Some(landing));
    assert_eq!(s.header.curl, None);
    assert_eq!(s.header.bed[2], [10.0, -0.2]);
    assert_eq!(s.field(1, 0), &[0.0, 5.5, 10.0]);
    assert_eq!(s.field(1, 1), &[0.0, 0.32, 0.01]);
    assert_eq!(s.field(1, 2), &[0.0, 1.75, 0.1]);

    // A frames file cut short is caught.
    let frames = dir.join("frames.f32");
    let bytes = std::fs::read(&frames).unwrap();
    std::fs::write(&frames, &bytes[..bytes.len() - 4]).unwrap();
    assert!(Slice::read(&dir).is_err());
    std::fs::remove_dir_all(&dir).unwrap();
}
