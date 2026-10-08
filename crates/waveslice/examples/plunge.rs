//! A solitary wave breaking on a plane slope, for looking at: Grilli et al.'s (1997) setting.
//!
//!     cargo run --release -p waveslice --example plunge -- [slope 1:N] [height/depth] [spacing] [out.csv]
//!
//! Prints the wave as it goes and writes surface snapshots (step, time, x, z) to the CSV.

use std::io::Write;
use waveslice::{Solitary, Tank};

const G: f64 = 9.81;
/// Re-space the surface nodes when neighbouring gaps differ by more than this factor.
const REGRID_ABOVE: f64 = 1.3;
/// Smooth the surface every this many steps.
const SMOOTH_EVERY: usize = 5;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let arg = |i: usize, default: f64| args.get(i).and_then(|a| a.parse().ok()).unwrap_or(default);
    let (run, height, finest) = (arg(1, 15.0), arg(2, 0.3), arg(3, 0.05));
    let out = args.get(4).cloned().unwrap_or_else(|| "plunge.csv".into());
    // A hard stop on computing time, so a run that goes wrong cannot keep the machine busy.
    let limit = arg(5, 1500.0);
    // Print every step from this one on, to see what happens just before a failure.
    let trace_from = arg(6, f64::INFINITY);
    let (depth, shallow, toe) = (1.0, 0.05, 15.0);
    let bed_depth = |x: f64| {
        if x <= toe {
            depth
        } else {
            depth - (x - toe) / run
        }
    };
    let wave = Solitary {
        g: G,
        depth,
        height: height * depth,
        crest: 7.5,
    };
    let mut t = Tank::beach(&wave, toe, run, shallow, finest);
    let s0 = 1.521 / run / height.sqrt();
    println!(
        "slope 1:{run}, H0/h0 {height}, S0 {s0:.3}; {} surface nodes, {} bed nodes",
        t.surface.len(),
        t.bed.len()
    );

    let mut csv = std::fs::File::create(&out).expect("the output file");
    writeln!(csv, "step,time,x,z").unwrap();
    let e0 = t.energy().total();
    let clock = std::time::Instant::now();
    let mut step = 0;
    let mut regrids = 0;
    loop {
        if step % 10 == 0 {
            for p in &t.surface {
                writeln!(csv, "{step},{:.5},{:.6},{:.6}", t.time, p.0, p.1).unwrap();
            }
        }
        if step % 20 == 0 {
            let top = t
                .surface
                .iter()
                .copied()
                .fold((0.0, f64::MIN), |a, p| if p.1 > a.1 { p } else { a });
            let steepest = t
                .surface
                .windows(2)
                .map(|w| ((w[1].1 - w[0].1) / (w[1].0 - w[0].0)).atan().to_degrees())
                .fold(0.0_f64, |a, b| if b.abs() > a.abs() { b } else { a });
            println!(
                "step {step:5} t {:.3} s crest {:.3} m at x {:.2} (depth {:.3}); steepest face {:.0} deg; energy {:+.2e}; {:.0} s",
                t.time,
                top.1,
                top.0,
                bed_depth(top.0),
                steepest,
                t.energy().total() / e0 - 1.0,
                clock.elapsed().as_secs_f64()
            );
        }
        if t.landed() {
            println!("the lip landed at t {:.3} s", t.time);
            break;
        }
        if t.tangled() {
            println!("FAILED: neighbouring nodes crossed at t {:.3} s", t.time);
            break;
        }
        if t.surface.last().unwrap().1 > 0.5 * depth {
            println!("the wave reached the shallow wall before breaking");
            break;
        }
        let before = clock.elapsed().as_secs_f64();
        let dt = t.advance(0.5);
        step += 1;
        if step % SMOOTH_EVERY == 0 {
            t.smooth();
        }
        let mut passes = 0;
        while t.unevenness() > REGRID_ABOVE && passes < 5 {
            t.regrid();
            passes += 1;
            regrids += 1;
        }
        if step as f64 >= trace_from {
            let gap = t
                .surface
                .windows(2)
                .map(|w| (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1))
                .fold(f64::INFINITY, f64::min);
            println!(
                "  step {step} dt {dt:.2e} closest {gap:.2e} unevenness {:.2} passes {passes} took {:.2} s",
                t.unevenness(),
                clock.elapsed().as_secs_f64() - before
            );
        }
        if !t.surface.iter().all(|p| p.0.is_finite() && p.1.is_finite()) {
            println!("FAILED: the surface is no longer finite at step {step}");
            break;
        }
        if clock.elapsed().as_secs_f64() > limit {
            println!("stopped: over the {limit} s computing limit at step {step}");
            break;
        }
        if dt < 1e-5 || step > 20_000 {
            println!("stopped: step {dt:.2e} s after {step} steps");
            break;
        }
    }
    for p in &t.surface {
        writeln!(csv, "{step},{:.5},{:.6},{:.6}", t.time, p.0, p.1).unwrap();
    }
    println!(
        "overturned: {}; {regrids} re-spacing passes; {:.0} s of computing",
        t.overturned(),
        clock.elapsed().as_secs_f64()
    );
}
