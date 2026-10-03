//! Verification ladder, rungs 3-5: wave generation, shoaling, refraction, friction and
//! runup, each against a known answer.

mod common;

use common::{col, run_and_analyse};
use wavecore::{Bathymetry, G, GHOST, Grid, Solver, Sponge, State, WaveMaker, manning_factor};

#[test]
fn wave_maker_radiates_the_requested_amplitude() {
    let dx = 0.5;
    let grid = Grid::new(800, 4, dx, dx);
    let bed = Bathymetry::flat(&grid, 2.0);
    let period = 20.0;
    let amplitude = 0.05;
    let maker = WaveMaker {
        x: 200.0,
        amplitude,
        period,
        angle: 0.0,
        sigma: 2.0,
        ramp_periods: 2.0,
    };
    let solver = Solver::new(grid, bed)
        .with_wavemaker(maker)
        .with_sponge(Sponge::new(120, 1.0).edges(true, true, false, false));
    let mut s = State::lake_at_rest(&grid, &solver.bed, 0.0);

    let j = GHOST + 2;
    let probes = [(col(&grid, 120.0), j), (col(&grid, 290.0), j)];
    let omega = std::f64::consts::TAU / period;
    let out = run_and_analyse(&solver, &mut s, 180.0, omega, (100.0, 180.0), &probes);

    for (name, (amp, _)) in ["west", "east"].iter().zip(&out) {
        let err = (amp - amplitude).abs() / amplitude;
        println!(
            "{name}: amplitude {amp:.4} (target {amplitude}), error {:.1}%",
            err * 100.0
        );
        assert!(
            err < 0.02,
            "{name} amplitude {amp} is {:.1}% off",
            err * 100.0
        );
    }
}

#[test]
fn shoaling_follows_greens_law() {
    // Long waves climbing a gentle slope: amplitude grows as depth^(-1/4). This is a
    // linear result, so the wave is kept small (a/h < 0.4%): at 3 cm, nonlinear
    // steepening of shallow-water waves drains a few percent of the first harmonic
    // over the probe spacing.
    let dx = 0.5;
    let grid = Grid::new(1200, 4, dx, dx);
    let bed = Bathymetry::from_fn(&grid, |x, _| -3.0 + x / 300.0);
    let period = 15.0;
    let maker = WaveMaker {
        x: 150.0,
        amplitude: 0.005,
        period,
        angle: 0.0,
        sigma: 2.0,
        ramp_periods: 2.0,
    };
    let solver = Solver::new(grid, bed.clone())
        .with_wavemaker(maker)
        .with_sponge(Sponge::new(120, 1.0).edges(true, true, false, false));
    let mut s = State::lake_at_rest(&grid, &solver.bed, 0.0);

    let j = GHOST + 2;
    let (i1, i2) = (col(&grid, 240.0), col(&grid, 400.0));
    let omega = std::f64::consts::TAU / period;
    let out = run_and_analyse(
        &solver,
        &mut s,
        180.0,
        omega,
        (120.0, 180.0),
        &[(i1, j), (i2, j)],
    );

    let (h1, h2) = (-bed.b[grid.idx(i1, j)], -bed.b[grid.idx(i2, j)]);
    let expected = (h1 / h2).powf(0.25);
    let measured = out[1].0 / out[0].0;
    println!(
        "Green's law: depths {h1:.3} -> {h2:.3} m, expected ratio {expected:.4}, measured {measured:.4}"
    );
    assert!(
        (measured / expected - 1.0).abs() < 0.02,
        "amplitude ratio {measured:.4} vs Green's law {expected:.4}"
    );
}

#[test]
fn oblique_waves_obey_snells_law_and_wave_action_conservation() {
    // Long waves on a plane beach, periodic alongshore so the wave is uniform in y.
    // Snell: k_y = k sin(theta) is conserved while the wave turns towards the shore.
    // Wave action: a^2 c cos(theta) is constant, so amplitude follows from the local
    // speed and angle.
    let (dx, dy) = (0.5, 2.5);
    let (nx, ny) = (900, 40);
    let grid = Grid::new(nx, ny, dx, dy);
    let bed = Bathymetry::from_fn(&grid, |x, _| -3.0 + x / 250.0);
    let period = 10.0;
    let omega = std::f64::consts::TAU / period;
    let source_x = 70.0;
    let c_source = (G * -bed.b[grid.idx(col(&grid, source_x), GHOST)]).sqrt();
    // One whole alongshore wavelength fits the period Ly = ny * dy = 100 m.
    let ky = std::f64::consts::TAU / (ny as f64 * dy);
    let angle = (ky * c_source / omega).asin();
    let maker = WaveMaker {
        x: source_x,
        amplitude: 0.005,
        period,
        angle,
        sigma: 1.0,
        ramp_periods: 2.0,
    };
    let solver = Solver::new(grid, bed.clone())
        .with_periodic_y()
        .with_wavemaker(maker)
        .with_sponge(Sponge::new(100, 1.5).edges(true, true, false, false));
    let mut s = State::lake_at_rest(&grid, &bed, 0.0);

    let dj = 10; // 25 m alongshore
    let cells = |x: f64| [(col(&grid, x), GHOST), (col(&grid, x), GHOST + dj)];
    let probes = [cells(140.0), cells(230.0)].concat();
    let out = run_and_analyse(&solver, &mut s, 120.0, omega, (80.0, 120.0), &probes);

    let expected_lag = ky * dj as f64 * dy;
    let depth = |x: f64| -bed.b[grid.idx(col(&grid, x), GHOST)];
    let mut turned = Vec::new();
    for (name, x, pair) in [("x=140", 140.0, &out[0..2]), ("x=230", 230.0, &out[2..4])] {
        let lag = (pair[0].1 - pair[1].1).rem_euclid(std::f64::consts::TAU);
        let c = (G * depth(x)).sqrt();
        let theta = (ky * c / omega).asin();
        println!(
            "{name}: alongshore phase lag {lag:.4} rad (expected {expected_lag:.4}), \
             depth {:.3} m, wave angle {:.1} deg, amplitude {:.4}, y-uniformity {:.4}",
            depth(x),
            theta.to_degrees(),
            pair[0].0,
            pair[1].0 / pair[0].0
        );
        assert!(
            (lag - expected_lag).abs() < 0.05,
            "{name}: k_y not conserved"
        );
        assert!(
            (pair[1].0 / pair[0].0 - 1.0).abs() < 0.01,
            "{name}: not uniform alongshore"
        );
        turned.push((c, theta));
    }

    let (c1, t1) = turned[0];
    let (c2, t2) = turned[1];
    let expected = (c1 * t1.cos() / (c2 * t2.cos())).sqrt();
    let measured = out[2].0 / out[0].0;
    println!(
        "refraction + shoaling: expected amplitude ratio {expected:.4}, measured {measured:.4}"
    );
    assert!(
        (measured / expected - 1.0).abs() < 0.02,
        "amplitude ratio {measured:.4} vs wave-action law {expected:.4}"
    );
}

#[test]
fn friction_attenuates_waves_at_the_predicted_rate() {
    let dx = 0.5;
    let period = 20.0;
    let amplitude = 0.05;
    let depth = 2.0;
    let n = 0.15;
    let ratio = |manning: f64| {
        let grid = Grid::new(800, 4, dx, dx);
        let bed = Bathymetry::flat(&grid, depth);
        let maker = WaveMaker {
            x: 200.0,
            amplitude,
            period,
            angle: 0.0,
            sigma: 2.0,
            ramp_periods: 2.0,
        };
        let solver = Solver::new(grid, bed)
            .with_wavemaker(maker)
            .with_sponge(Sponge::new(120, 1.0).edges(true, true, false, false))
            .with_manning(manning);
        let mut s = State::lake_at_rest(&grid, &solver.bed, 0.0);
        let j = GHOST + 2;
        let probes = [(col(&grid, 230.0), j), (col(&grid, 330.0), j)];
        let omega = std::f64::consts::TAU / period;
        let out = run_and_analyse(&solver, &mut s, 180.0, omega, (100.0, 180.0), &probes);
        out[1].0 / out[0].0
    };

    let attenuation = ratio(n) / ratio(0.0);

    // Quadratic drag acts on the waves like a linear drag r = (8 / 3 pi) k u_max,
    // with k = g n^2 / h^(4/3) and u_max = a c / h. Amplitude decays as exp(-r t / 2)
    // over the travel time t = distance / c.
    let c = (G * depth).sqrt();
    let k = G * n * n / depth.powf(4.0 / 3.0);
    let r = 8.0 / (3.0 * std::f64::consts::PI) * k * amplitude * c / depth;
    let predicted = (-0.5 * r * 100.0 / c).exp();
    println!(
        "friction attenuation over 100 m: measured {attenuation:.4}, predicted {predicted:.4}"
    );
    assert!(
        (attenuation - predicted).abs() < 0.03,
        "attenuation {attenuation:.4} vs predicted {predicted:.4}"
    );
}

#[test]
fn manning_factor_is_the_exact_quadratic_drag_step() {
    // du/dt = -k u |u| has the solution u0 / (1 + k u0 t) in one direction.
    let (h, u0, dt, n) = (1.5_f64, 0.8, 0.2, 0.03);
    let k = G * n * n / h.powf(4.0 / 3.0);
    let exact = u0 / (1.0 + k * u0 * dt);
    assert!((manning_factor(h, u0, dt, n) * u0 - exact).abs() < 1e-15);
    assert_eq!(manning_factor(h, u0, dt, 0.0), 1.0);
    assert!(
        manning_factor(h, 100.0, 100.0, 0.1) > 0.0,
        "friction must never reverse the flow"
    );
}

#[test]
fn solitary_wave_runup_matches_synolakis_law() {
    // Synolakis (1987): a non-breaking solitary wave of height H on a 1:19.85 beach
    // in depth d runs up to R/d = 2.831 sqrt(cot(beta)) (H/d)^(5/4).
    let (d, h_wave, cot_beta) = (1.0, 0.0185, 19.85);
    let toe = 70.0;
    let dx = 0.05;
    let grid = Grid::new((96.0 / dx) as usize, 3, dx, dx);
    let bed = Bathymetry::from_fn(&grid, |x, _| {
        if x < toe {
            -d
        } else {
            -d + (x - toe) / cot_beta
        }
    });

    let gamma = (3.0 * h_wave / (4.0 * d.powi(3))).sqrt();
    let crest = toe - 20.0_f64.sqrt().acosh() / gamma;
    let mut s = State::zeros(&grid);
    for (i, j) in grid.interior() {
        let k = grid.idx(i, j);
        let (x, _) = grid.centre(i, j);
        let eta = h_wave / (gamma * (x - crest)).cosh().powi(2);
        let depth = (eta - bed.b[k]).max(0.0);
        s.h[k] = depth;
        // Right-moving soliton, first-order velocity (Synolakis 1987).
        s.hu[k] = depth * (G / d).sqrt() * eta * (1.0 - eta / (4.0 * d));
    }
    s.fill_walls(&grid);
    let v0 = s.volume(&grid);
    let solver = Solver::new(grid, bed);

    let mut runup: f64 = 0.0;
    while s.time < 30.0 {
        let dt = solver.stable_dt(&s).min(30.0 - s.time);
        solver.step(&mut s, dt);
        for (i, j) in grid.interior().filter(|&(_, j)| j == GHOST + 1) {
            let k = grid.idx(i, j);
            if s.h[k] > 1e-3 {
                runup = runup.max(s.h[k] + solver.bed.b[k]);
            }
        }
    }

    let expected = 2.831 * cot_beta.sqrt() * (h_wave / d).powf(1.25) * d;
    println!(
        "runup: measured {runup:.4} m, Synolakis law {expected:.4} m, error {:.1}%",
        (runup / expected - 1.0) * 100.0
    );
    assert!(
        (s.volume(&grid) - v0).abs() < 1e-9 * v0,
        "volume not conserved"
    );
    assert!(s.min_depth(&grid) >= 0.0);
    assert!(
        (runup / expected - 1.0).abs() < 0.10,
        "runup {runup:.4} m vs law {expected:.4} m"
    );
}

#[test]
fn a_sponge_over_land_never_speeds_the_water_up() {
    // Regression: a sponge relaxed the surface towards its level, which on land above that
    // level drained the swash faster than its momentum, so the water sped up. In the side
    // sponges of the Pipeline run it left films microns deep moving at tens of metres per
    // second, which cut the time step eightfold. A sponge relaxes depth and momentum at the
    // same rate, so a uniform sheet of water on flat land must keep its velocity exactly.
    let grid = Grid::new(20, 20, 2.0, 2.0);
    let bed = Bathymetry::flat(&grid, -0.5);
    let mut s = State::zeros(&grid);
    for (i, j) in grid.interior() {
        let k = grid.idx(i, j);
        s.h[k] = 1.0;
        s.hu[k] = 1.0;
        s.hv[k] = 0.5;
    }
    // A sponge so wide that its rate is the same everywhere here to within 0.004%, so the land
    // drains evenly and no slope builds up to push the water.
    let solver = Solver::new(grid, bed)
        .with_sponge(Sponge::new(1_000_000, 1.5).edges(true, false, false, false));
    for _ in 0..3 {
        solver.step(&mut s, 0.05);
    }
    // The walls reach four cells inward per step; look at the middle, well away from them.
    for i in GHOST + 7..GHOST + 13 {
        for j in GHOST + 7..GHOST + 13 {
            let k = grid.idx(i, j);
            assert!(
                s.h[k] < 1.0,
                "the sponge did not drain the land at ({i}, {j})"
            );
            let (u, v) = (s.hu[k] / s.h[k], s.hv[k] / s.h[k]);
            assert!(
                (u - 1.0).abs() < 1e-6 && (v - 0.5).abs() < 1e-6,
                "the sponge changed the velocity at ({i}, {j}) to ({u}, {v})"
            );
        }
    }
}
