//! Stage 3 of the slice solver: solitary waves on a plane slope, against Grilli, Svendsen &
//! Subramanya (1997), who computed them with a potential-flow model checked against laboratory
//! waves. Their slope parameter, `S0 = 1.521 s / sqrt(H0 / h0)` for a slope `s` and a wave `H0`
//! high in water `h0` deep before the slope, sorts the breakers: spilling below 0.025, plunging
//! to 0.30, surging to 0.37, and above 0.37 the wave does not break. The bounds were set before
//! the first run.

use std::time::Instant;
use waveslice::{Solitary, Tank};

const G: f64 = 9.81;
/// The slope rises 1 in this many.
const RUN: f64 = 15.0;
const TOE: f64 = 15.0;
const SHALLOW: f64 = 0.05;
/// Surface nodes no closer than this where the water is shallow, m.
const FINEST: f64 = 0.05;
/// Stop a test that computes for longer than this, s: a run that goes wrong must not keep the
/// machine busy.
const LIMIT: f64 = 1800.0;

/// A solitary wave `height` high in water 1 m deep, on the beach, its crest 7.5 m from the left
/// wall and the toe of the slope at 15 m.
fn beach(height: f64) -> Tank {
    let wave = Solitary {
        g: G,
        depth: 1.0,
        height,
        crest: 7.5,
    };
    Tank::beach(&wave, TOE, RUN, SHALLOW, FINEST)
}

fn depth_at(x: f64) -> f64 {
    if x <= TOE { 1.0 } else { 1.0 - (x - TOE) / RUN }
}

/// The highest point of the surface.
fn crest(t: &Tank) -> (f64, f64) {
    t.surface
        .iter()
        .copied()
        .fold((0.0, f64::MIN), |a, p| if p.1 > a.1 { p } else { a })
}

#[test]
#[ignore = "slow physics check: cargo test --release -- --include-ignored"]
fn a_solitary_wave_plunges_on_a_one_in_fifteen_slope() {
    // Grilli et al.'s own plunging case, which Derakhti et al. (2020) computed again: S0 = 0.19.
    // Its lip lands a metre from the wall that ends this beach; `wavesim slice steep`, the same
    // case eight times larger with the beach ending in a shelf 24 m long instead, lands it 2 cm
    // further on and otherwise the same, so the wall does not shape it.
    let mut t = beach(0.3);
    let e0 = t.energy().total();
    let clock = Instant::now();
    let (mut count, mut overturned_at) = (0, None);
    while !t.landed() {
        assert!(
            !t.tangled(),
            "neighbouring nodes crossed at {:.3} s",
            t.time
        );
        assert!(
            clock.elapsed().as_secs_f64() < LIMIT,
            "over {LIMIT} s of computing"
        );
        assert!(
            t.surface[t.surface.len() - 1].1 < 0.5,
            "the wave reached the shallow wall without landing"
        );
        count += 1;
        t.march();
        if overturned_at.is_none() && t.overturned() {
            let (x, z) = crest(&t);
            overturned_at = Some((t.time, x, z, depth_at(x)));
        }
    }
    let drift = t.energy().total() / e0 - 1.0;
    let (time, x, height, depth) = overturned_at.expect("the lip landed without overturning");
    println!(
        "overturned at {time:.3} s, crest {height:.3} m at x {x:.2} m in {depth:.3} m of water \
         (H/h {:.2}); landed at {:.3} s; energy {:+.3}% at landing; {count} steps, {:.0} s",
        height / depth,
        t.time,
        100.0 * drift,
        clock.elapsed().as_secs_f64()
    );
    assert!(
        drift.abs() < 5e-3,
        "energy changed by {:.2}%",
        100.0 * drift
    );
}

#[test]
#[ignore = "slow physics check: cargo test --release -- --include-ignored"]
fn a_low_solitary_wave_does_not_break_on_the_same_slope() {
    // 0.07 of the depth: S0 = 0.38, above 0.37, and below Grilli et al.'s largest wave that does
    // not break on a slope s, H0/h0 = 16.9 s^2 = 0.075. It must not overturn while its crest is
    // in water deeper than 0.1 m, twice the depth at the wall that ends this beach.
    let mut t = beach(0.07);
    let e0 = t.energy().total();
    let clock = Instant::now();
    let mut count = 0;
    let mut steepest = 0.0_f64;
    while depth_at(crest(&t).0) > 0.1 {
        assert!(
            !t.overturned(),
            "it overturned at {:.3} s, crest at {:?}",
            t.time,
            crest(&t)
        );
        assert!(
            !t.tangled(),
            "neighbouring nodes crossed at {:.3} s",
            t.time
        );
        assert!(
            clock.elapsed().as_secs_f64() < LIMIT,
            "over {LIMIT} s of computing"
        );
        count += 1;
        t.march();
        let face = t
            .surface
            .windows(2)
            .map(|w| (w[0].1 - w[1].1).atan2(w[1].0 - w[0].0).to_degrees())
            .fold(0.0_f64, f64::max);
        steepest = steepest.max(face);
    }
    let drift = t.energy().total() / e0 - 1.0;
    println!(
        "reached 0.1 m of water at {:.3} s without overturning; steepest front face {steepest:.0} \
         deg; energy {:+.3}%; {count} steps, {:.0} s",
        t.time,
        100.0 * drift,
        clock.elapsed().as_secs_f64()
    );
    assert!(
        drift.abs() < 5e-3,
        "energy changed by {:.2}%",
        100.0 * drift
    );
}
