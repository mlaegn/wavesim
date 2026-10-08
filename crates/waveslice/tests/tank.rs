//! Stage 2 of the slice solver: the surface moving in time, against what is known exactly. The
//! bounds were set before the first run.

use std::f64::consts::{PI, TAU};
use waveslice::Tank;

const G: f64 = 9.81;

/// `n` nodes spread evenly from `a` to `b`.
fn spread(n: usize, a: f64, b: f64) -> Vec<f64> {
    (0..n)
        .map(|k| a + (b - a) * k as f64 / (n - 1) as f64)
        .collect()
}

/// A flat-bottomed tank `length` long and `depth` deep, with the surface at `eta(x)` and the
/// potential `phi(x)` on it.
fn tank(
    length: f64,
    depth: f64,
    nodes: (usize, usize, usize),
    eta: impl Fn(f64) -> f64,
    phi: impl Fn(f64) -> f64,
) -> Tank {
    let (surface, bed, wall) = nodes;
    let xs = spread(surface, 0.0, length);
    Tank::new(
        G,
        spread(bed, 0.0, length)
            .into_iter()
            .map(|x| (x, -depth))
            .collect(),
        wall,
        xs.iter().map(|&x| (x, eta(x))).collect(),
        xs.iter().map(|&x| phi(x)).collect(),
    )
}

/// A standing wave one wavelength long in a basin 10 m long and 3 m deep, released from rest.
fn standing(amplitude: f64, surface: usize) -> (Tank, f64) {
    let (length, depth) = (10.0, 3.0);
    let k = TAU / length;
    let period = TAU / (G * k * (k * depth).tanh()).sqrt();
    let t = tank(
        length,
        depth,
        (surface, 41, 13),
        |x| amplitude * (k * x).cos(),
        |_| 0.0,
    );
    (t, period)
}

#[test]
fn a_small_standing_wave_has_the_period_of_linear_theory() {
    // 1 mm on a 10 m wavelength: the nonlinear change of period is of order (ka)^2, 4e-7.
    let (mut t, period) = standing(0.001, 41);
    let dt = period / 80.0;
    // The times the surface at the left wall crosses still level, between steps by a straight
    // line, for three periods.
    let mut crossings = Vec::new();
    let mut before = t.surface[0].1;
    while t.time < 3.2 * period {
        t.step(dt);
        let now = t.surface[0].1;
        if before.signum() != now.signum() {
            crossings.push(t.time - dt * now / (now - before));
        }
        before = now;
    }
    let measured =
        2.0 * (crossings[crossings.len() - 1] - crossings[0]) / (crossings.len() - 1) as f64;
    let error = measured / period - 1.0;
    println!(
        "period {measured:.6} s against {period:.6} s from linear theory ({:+.4}%), {} crossings",
        100.0 * error,
        crossings.len()
    );
    assert!(crossings.len() == 6, "{crossings:?}");
    assert!(error.abs() < 5e-4, "period off by {:.4}%", 100.0 * error);
}

#[test]
#[ignore = "slow physics check: cargo test --release -- --include-ignored"]
fn a_steep_standing_wave_keeps_its_energy_and_volume() {
    // ka = 0.25: crests 0.4 m on a 3 m depth, well into the nonlinear range.
    let amplitude = 0.4;
    let (mut t, period) = standing(amplitude, 61);
    let dt = period / 100.0;
    let (e0, v0) = (t.energy().total(), t.volume());
    // The water above still level in half a wavelength, to measure volume against.
    let crest_volume = amplitude * 10.0 / PI;
    let (mut worst_energy, mut worst_volume) = (0.0_f64, 0.0_f64);
    for _ in 0..4 {
        for _ in 0..100 {
            t.step(dt);
        }
        let (e, v) = (t.energy().total(), t.volume());
        worst_energy = worst_energy.max((e / e0 - 1.0).abs());
        worst_volume = worst_volume.max(((v - v0) / crest_volume).abs());
    }
    println!(
        "over 4 periods: energy within {worst_energy:.2e}, volume within {worst_volume:.2e} of a crest's"
    );
    assert!(worst_energy < 1e-3, "energy drifted by {worst_energy:.2e}");
    assert!(
        worst_volume < 1e-3,
        "volume drifted by {worst_volume:.2e} of a crest's"
    );
}

#[test]
#[ignore = "slow physics check: cargo test --release -- --include-ignored"]
fn a_solitary_wave_travels_at_the_fully_nonlinear_speed() {
    // A solitary wave 0.2 of the depth, started from the Serre-Green-Naghdi solitary wave: its
    // surface, and on it the velocity of that theory's profile in depth, u = ubar - H^2/3 ubar_xx
    // along and w = -H ubar_x up, H the local depth. That is close to but not exactly the
    // potential-flow wave: it sheds a small tail and settles. (Started instead with the
    // depth-averaged velocity as the surface's, about 4% slow under the crest, the crest still
    // swung by 1.8% while measured and the speed came out 0.085% high; the bounds were kept.)
    // Its speed is then
    // compared with Grimshaw's (1971) third-order speed for the amplitude it settles at,
    // c^2 = g h (1 + e - e^2/20 - 3 e^3/70), within 0.01% of the exact speed at this amplitude.
    // The Serre-Green-Naghdi speed, c^2 = g (h + a), is 0.1% higher, so the bound tells them
    // apart.
    let (depth, length, start): (f64, f64, f64) = (1.0, 60.0, 15.0);
    let a = 0.2 * depth;
    let kappa = (3.0 * a / (4.0 * depth * depth * (depth + a))).sqrt();
    let c0 = (G * (depth + a)).sqrt();
    let eta = move |x: f64| a / (kappa * (x - start)).cosh().powi(2);
    // phi along the surface is the integral of u + w eta_x, by the trapezoid rule on a fine grid,
    // with the derivatives by centred differences.
    let ubar = move |x: f64| c0 * eta(x) / (depth + eta(x));
    let e = 1e-3;
    let d1 = |f: &dyn Fn(f64) -> f64, x: f64| (f(x + e) - f(x - e)) / (2.0 * e);
    let d2 = |f: &dyn Fn(f64) -> f64, x: f64| (f(x + e) - 2.0 * f(x) + f(x - e)) / (e * e);
    let slope = |x: f64| {
        let big = depth + eta(x);
        let u = ubar(x) - big * big / 3.0 * d2(&ubar, x);
        let w = -big * d1(&ubar, x);
        u + w * d1(&eta, x)
    };
    let fine = spread(30_001, 0.0, length);
    let mut potential = vec![0.0; fine.len()];
    for i in 1..fine.len() {
        potential[i] = potential[i - 1]
            + 0.5 * (slope(fine[i - 1]) + slope(fine[i])) * (fine[i] - fine[i - 1]);
    }
    let phi = |x: f64| potential[(x / length * 30_000.0).round() as usize];
    let mut t = tank(length, depth, (301, 151, 9), eta, phi);
    let dt = 0.05;
    // The crest, from the parabola through the highest node and its neighbours, in Newton's
    // form z = p + d1 (x - xp) + curve (x - xp)(x - xq).
    let crest = |t: &Tank| {
        let s = &t.surface;
        let i = (1..s.len() - 1)
            .max_by(|&i, &j| s[i].1.total_cmp(&s[j].1))
            .unwrap();
        let (p, q, r) = (s[i - 1], s[i], s[i + 1]);
        let d1 = (q.1 - p.1) / (q.0 - p.0);
        let curve = ((r.1 - q.1) / (r.0 - q.0) - d1) / (r.0 - p.0);
        let x = 0.5 * (p.0 + q.0) - d1 / (2.0 * curve);
        (x, p.1 + d1 * (x - p.0) + curve * (x - p.0) * (x - q.0))
    };
    // Settle over the first 10 depths of travel, then measure over the next 20.
    let mut samples = Vec::new();
    while crest(&t).0 < start + 30.0 {
        t.step(dt);
        let (x, z) = crest(&t);
        if x > start + 10.0 {
            samples.push((t.time, x, z));
        }
    }
    let n = samples.len() as f64;
    let mean = |f: &dyn Fn(&(f64, f64, f64)) -> f64| samples.iter().map(f).sum::<f64>() / n;
    let (tm, xm) = (mean(&|s| s.0), mean(&|s| s.1));
    let speed = samples.iter().map(|s| (s.0 - tm) * (s.1 - xm)).sum::<f64>()
        / samples.iter().map(|s| (s.0 - tm).powi(2)).sum::<f64>();
    let height = mean(&|s| s.2);
    let (low, high) = samples.iter().fold((f64::INFINITY, 0.0_f64), |(l, h), s| {
        (l.min(s.2), h.max(s.2))
    });
    let e = height / depth;
    let grimshaw = (G * depth * (1.0 + e - e * e / 20.0 - 3.0 * e.powi(3) / 70.0)).sqrt();
    let sgn = (G * (depth + height)).sqrt();
    let error = speed / grimshaw - 1.0;
    println!(
        "settled at a/h {e:.4} (from {low:.4} to {high:.4} m); speed {speed:.5} m/s, Grimshaw {grimshaw:.5} ({:+.3}%), Serre-Green-Naghdi {sgn:.5} ({:+.3}%)",
        100.0 * error,
        100.0 * (speed / sgn - 1.0)
    );
    assert!(
        error.abs() < 5e-4,
        "speed off Grimshaw's by {:.3}%",
        100.0 * error
    );
    assert!(
        (high - low) / height < 0.01,
        "the crest changed by {:.2}% while measured",
        100.0 * (high - low) / height
    );
}
