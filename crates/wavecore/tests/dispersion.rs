//! Verification ladder, rung 6: the dispersive (Serre-Green-Naghdi) model.
//!
//! Where the model has an exact answer (its dispersion relation, its solitary wave) the
//! test compares against it. Where it has none (breaking), the acceptance bounds come from
//! the textbook range for the physical quantity and were fixed before the first run.

mod common;

use std::f64::consts::PI;

use common::{col, run_and_analyse};
use wavecore::{
    Bathymetry, Dispersion, G, GHOST, Grid, LinearWave, Order, Solver, Sponge, State, WaveMaker,
    linear_wave,
};

fn sgn_omega(h: f64, k: f64) -> f64 {
    (G * h * k * k / (1.0 + (k * h).powi(2) / 3.0)).sqrt()
}

fn run_to(solver: &Solver, s: &mut State, t_end: f64) {
    while s.time < t_end {
        let dt = solver.stable_dt(s);
        assert!(dt.is_finite(), "no water in the domain");
        solver.step(s, dt.min(t_end - s.time));
    }
}

/// A solver with the given scheme order, with or without the dispersive correction.
fn model(solver: Solver, dispersive: bool, order: Order) -> Solver {
    let solver = solver.with_order(order);
    if dispersive {
        solver.with_dispersion(Dispersion::default())
    } else {
        solver
    }
}

/// Shallow water with the MC limiter (the baseline), or SGN with the third-order scheme.
fn with_model(solver: Solver, dispersive: bool) -> Solver {
    let order = if dispersive {
        Order::Third
    } else {
        Order::Second
    };
    model(solver, dispersive, order)
}

const ORDERS: [Order; 2] = [Order::Second, Order::Third];

// ---- 1. The dispersion relation -------------------------------------------------------

/// Angular frequency of a standing wave in a closed basin, from its zero crossings.
fn standing_wave_omega(dispersive: bool, order: Order) -> f64 {
    let (depth, length, mode, amplitude) = (5.0, 100.0, 8.0, 0.01);
    let dx = 0.5;
    let grid = Grid::new(200, 4, dx, dx);
    let solver = model(
        Solver::new(grid, Bathymetry::flat(&grid, depth)),
        dispersive,
        order,
    );
    let k = mode * PI / length;
    let mut s = State::zeros(&grid);
    for (i, j) in grid.interior() {
        let (x, _) = grid.centre(i, j);
        s.h[grid.idx(i, j)] = depth + amplitude * (k * x).cos();
    }

    let probe = grid.idx(GHOST, GHOST + 2);
    let mut crossings = Vec::new();
    let mut previous = (0.0, amplitude * (k * 0.25).cos());
    while s.time < 14.0 {
        let dt = solver.stable_dt(&s);
        solver.step(&mut s, dt);
        let eta = s.h[probe] - depth;
        if previous.1 * eta < 0.0 {
            let f = previous.1 / (previous.1 - eta);
            crossings.push(previous.0 + f * (s.time - previous.0));
        }
        previous = (s.time, eta);
    }
    assert!(
        crossings.len() >= 5,
        "only {} zero crossings",
        crossings.len()
    );
    (crossings.len() - 1) as f64 * PI / (crossings[crossings.len() - 1] - crossings[0])
}

#[test]
fn a_standing_wave_oscillates_at_the_sgn_frequency() {
    let (depth, k) = (5.0, 8.0 * PI / 100.0);
    let sgn = sgn_omega(depth, k);
    let shallow = (G * depth).sqrt() * k;
    let exact = (G * k * (k * depth).tanh()).sqrt();

    let measured_shallow = standing_wave_omega(false, Order::Second);
    assert!(
        (measured_shallow / shallow - 1.0).abs() < 0.01,
        "shallow water {measured_shallow} vs theory {shallow}"
    );
    for order in ORDERS {
        let measured = standing_wave_omega(true, order);
        println!(
            "omega ({order:?}): SGN model {measured:.4} (theory {sgn:.4}), shallow water \
             {measured_shallow:.4} (theory {shallow:.4}), full linear theory {exact:.4}"
        );
        assert!(
            (measured / sgn - 1.0).abs() < 0.01,
            "{order:?}: dispersive {measured} vs theory {sgn}"
        );
        assert!((measured / exact - 1.0).abs() < 0.03);
    }
    // Dispersion matters here: shallow water is far off the true frequency, SGN is close.
    assert!((shallow / exact - 1.0).abs() > 0.15);
}

#[test]
fn the_linear_wave_relation_inverts_the_dispersion_relation() {
    for (h, period) in [(2.0, 10.0), (5.0, 8.0), (8.0, 14.0), (0.5, 4.0)] {
        let omega = 2.0 * PI / period;
        let LinearWave {
            k,
            phase_speed,
            group_speed,
        } = linear_wave(omega, h, true);
        assert!(
            (sgn_omega(h, k) / omega - 1.0).abs() < 1e-12,
            "h={h} T={period}"
        );
        assert!((phase_speed - omega / k).abs() < 1e-12);
        // The group speed is d(omega)/dk.
        let eps = 1e-6 * k;
        let numerical = (sgn_omega(h, k + eps) - sgn_omega(h, k - eps)) / (2.0 * eps);
        assert!(
            (group_speed / numerical - 1.0).abs() < 1e-7,
            "group speed {group_speed} vs {numerical}"
        );
        assert!(
            group_speed < phase_speed,
            "energy travels slower than crests"
        );

        let shallow = linear_wave(omega, h, false);
        let c0 = (G * h).sqrt();
        assert!((shallow.phase_speed - c0).abs() < 1e-12);
        assert!((shallow.group_speed - c0).abs() < 1e-12);
        assert!((shallow.k * c0 - omega).abs() < 1e-12);
        assert!(phase_speed <= c0, "dispersion must slow waves down");
    }
}

// ---- 2. The solitary wave --------------------------------------------------------------

struct Solitary {
    peak: f64,
    peak_x: f64,
    shape_error: f64,
}

/// A solitary wave of amplitude 0.2 m in 1 m of water, run for 20 s.
fn solitary_wave(dispersive: bool, order: Order) -> (Solitary, f64, f64) {
    let (depth, amplitude, x0, t_end) = (1.0, 0.2, 20.0, 20.0);
    let dx = 0.1;
    let grid = Grid::new(1200, 4, dx, dx);
    let solver = model(
        Solver::new(grid, Bathymetry::flat(&grid, depth)),
        dispersive,
        order,
    );
    // The exact SGN solitary wave: speed sqrt(g (h + a)), width from a and h.
    let c = (G * (depth + amplitude)).sqrt();
    let kappa = (3.0 * amplitude / (4.0 * depth * depth * (depth + amplitude))).sqrt();
    let exact = |x: f64, t: f64| amplitude / (kappa * (x - x0 - c * t)).cosh().powi(2);

    let mut s = State::zeros(&grid);
    for (i, j) in grid.interior() {
        let (x, _) = grid.centre(i, j);
        let eta = exact(x, 0.0);
        let k = grid.idx(i, j);
        s.h[k] = depth + eta;
        s.hu[k] = s.h[k] * c * eta / (depth + eta);
    }
    run_to(&solver, &mut s, t_end);

    let j = GHOST + 2;
    let row: Vec<f64> = (GHOST..GHOST + grid.nx)
        .map(|i| s.h[grid.idx(i, j)] - depth)
        .collect();
    let (n, peak) = row
        .iter()
        .copied()
        .enumerate()
        .fold((0, f64::MIN), |m, (n, v)| if v > m.1 { (n, v) } else { m });
    // Parabolic refinement of the peak position.
    let (a, b, d) = (row[n - 1], row[n], row[n + 1]);
    let offset = 0.5 * (a - d) / (a - 2.0 * b + d);
    let peak_x = (n as f64 + 0.5 + offset) * dx;

    let (mut num, mut den) = (0.0, 0.0);
    for (n, &v) in row.iter().enumerate() {
        let e = exact((n as f64 + 0.5) * dx, t_end);
        num += (v - e).powi(2);
        den += e * e;
    }
    let travelled = c * t_end;
    (
        Solitary {
            peak,
            peak_x,
            shape_error: (num / den).sqrt(),
        },
        x0 + travelled,
        amplitude,
    )
}

#[test]
fn a_solitary_wave_keeps_its_shape_and_speed() {
    let (shallow, expected_x, amplitude) = solitary_wave(false, Order::Second);
    for order in ORDERS {
        let (sgn, _, _) = solitary_wave(true, order);
        println!(
            "solitary wave after 20 s, {order:?} (expected crest at x = {expected_x:.2} m, \
             amplitude {amplitude}):\n  \
             SGN:           peak {:.4} at x = {:.2}, shape error {:.3}\n  \
             shallow water: peak {:.4} at x = {:.2}, shape error {:.3}",
            sgn.peak,
            sgn.peak_x,
            sgn.shape_error,
            shallow.peak,
            shallow.peak_x,
            shallow.shape_error
        );

        // Fixed in advance: amplitude within 5%, crest within 1 m of 68.6 m travelled,
        // shape within 8% in the L2 sense.
        assert!(
            (sgn.peak / amplitude - 1.0).abs() < 0.05,
            "{order:?}: peak {}",
            sgn.peak
        );
        assert!(
            (sgn.peak_x - expected_x).abs() < 1.0,
            "{order:?}: crest at {}",
            sgn.peak_x
        );
        assert!(
            sgn.shape_error < 0.08,
            "{order:?}: shape error {}",
            sgn.shape_error
        );
        // Shallow water has no solitary wave: it steepens, so it is much further off.
        assert!(shallow.shape_error > 2.0 * sgn.shape_error);
    }
}

// ---- 3. The wave-maker, with dispersion ------------------------------------------------

#[test]
fn the_wave_maker_makes_the_right_amplitude_and_wavelength_with_dispersion() {
    let (depth, period, amplitude) = (5.0, 8.0, 0.05);
    let dx = 0.5;
    let grid = Grid::new(1000, 4, dx, dx);
    let maker = WaveMaker {
        x: 200.0,
        amplitude,
        period,
        angle: 0.0,
        sigma: 2.0,
        ramp_periods: 2.0,
    };
    let solver = Solver::new(grid, Bathymetry::flat(&grid, depth))
        .with_dispersion(Dispersion::default())
        .with_wavemaker(maker)
        .with_sponge(Sponge::new(100, 1.0).edges(true, true, false, false));
    let mut s = State::lake_at_rest(&grid, &solver.bed, 0.0);

    let omega = 2.0 * PI / period;
    let (x1, x2) = (250.0, 290.0);
    let j = GHOST + 2;
    let out = run_and_analyse(
        &solver,
        &mut s,
        100.0,
        omega,
        (60.0, 100.0),
        &[(col(&grid, x1), j), (col(&grid, x2), j)],
    );

    let k_sgn = linear_wave(omega, depth, true).k;
    let k_shallow = linear_wave(omega, depth, false).k;
    let lag = (out[0].1 - out[1].1).rem_euclid(2.0 * PI);
    let expected = (k_sgn * (x2 - x1)).rem_euclid(2.0 * PI);
    let shallow = (k_shallow * (x2 - x1)).rem_euclid(2.0 * PI);
    println!(
        "amplitudes {:.4}, {:.4} (target {amplitude}); phase lag over 40 m = {lag:.3} rad \
         (SGN {expected:.3}, shallow water would give {shallow:.3})",
        out[0].0, out[1].0
    );
    for (name, (amp, _)) in ["near", "far"].iter().zip(&out) {
        assert!(
            (amp / amplitude - 1.0).abs() < 0.02,
            "{name} amplitude {amp}"
        );
    }
    assert!(
        (lag - expected).abs() < 0.1,
        "phase lag {lag} vs SGN {expected}"
    );
    assert!(
        (expected - shallow).abs() > 0.2,
        "the test cannot tell the models apart"
    );

    let stats = solver.dispersion_stats().expect("dispersive solver");
    println!(
        "solver work: {stats:?}, {:.1} iterations per solve",
        stats.iterations as f64 / stats.solves as f64
    );
    assert_eq!(stats.unconverged, 0, "some linear solves did not converge");
    assert!(
        (stats.iterations as f64 / stats.solves as f64) < 20.0,
        "linear solves are too slow: {stats:?}"
    );
}

// ---- 4. A big swell in deep water keeps its energy ---------------------------------------

/// What a probe saw over a window: crest-to-trough height, and the wave height
/// `2 sqrt(2) sigma` that has the same energy (it equals the height of a sinusoid).
#[derive(Clone, Copy, Debug)]
struct Seen {
    crest_to_trough: f64,
    equivalent_height: f64,
}

#[derive(Clone, Copy)]
struct Moments {
    n: f64,
    sum: f64,
    sum_sq: f64,
    lo: f64,
    hi: f64,
}

impl Moments {
    fn new() -> Self {
        Self {
            n: 0.0,
            sum: 0.0,
            sum_sq: 0.0,
            lo: f64::MAX,
            hi: f64::MIN,
        }
    }

    fn add(&mut self, v: f64) {
        self.n += 1.0;
        self.sum += v;
        self.sum_sq += v * v;
        self.lo = self.lo.min(v);
        self.hi = self.hi.max(v);
    }

    fn seen(&self) -> Seen {
        let mean = self.sum / self.n;
        let var = self.sum_sq / self.n - mean * mean;
        Seen {
            crest_to_trough: self.hi - self.lo,
            equivalent_height: 2.0 * 2.0_f64.sqrt() * var.sqrt(),
        }
    }
}

/// A 2.5 m, 14 s swell in 8 m of water, seen 100 m and 500 m from the wave-maker.
fn deep_water_swell(dispersive: bool) -> [Seen; 2] {
    let (depth, period, height) = (8.0, 14.0, 2.5);
    let dx = 1.5;
    let grid = Grid::new(1000, 4, dx, dx);
    let maker = WaveMaker {
        x: 400.0,
        amplitude: height / 2.0,
        period,
        angle: 0.0,
        sigma: 2.0 * dx,
        ramp_periods: 2.0,
    };
    let solver = with_model(
        Solver::new(grid, Bathymetry::flat(&grid, depth)),
        dispersive,
    )
    .with_wavemaker(maker)
    .with_sponge(Sponge::new(150, 1.5).edges(true, true, false, false));
    let mut s = State::lake_at_rest(&grid, &solver.bed, 0.0);

    let j = GHOST + 2;
    let probes = [col(&grid, 500.0), col(&grid, 900.0)];
    let mut seen = [Moments::new(); 2];
    while s.time < 170.0 {
        let dt = solver.stable_dt(&s);
        solver.step(&mut s, dt);
        if s.time > 128.0 {
            for (m, &p) in seen.iter_mut().zip(&probes) {
                m.add(s.h[grid.idx(p, j)] - depth);
            }
        }
    }
    [seen[0].seen(), seen[1].seen()]
}

#[test]
fn a_big_swell_keeps_its_energy_in_deep_water() {
    let [near, far] = deep_water_swell(true);
    let [near_swe, far_swe] = deep_water_swell(false);
    let kept = |a: Seen, b: Seen| (b.equivalent_height / a.equivalent_height).powi(2);
    println!(
        "2.5 m swell, 100 m and 500 m from the maker:\n  \
         SGN:           equivalent height {:.2} -> {:.2} m (energy kept {:.0}%), \
         crest-to-trough {:.2} -> {:.2} m\n  \
         shallow water: equivalent height {:.2} -> {:.2} m (energy kept {:.0}%), \
         crest-to-trough {:.2} -> {:.2} m",
        near.equivalent_height,
        far.equivalent_height,
        100.0 * kept(near, far),
        near.crest_to_trough,
        far.crest_to_trough,
        near_swe.equivalent_height,
        far_swe.equivalent_height,
        100.0 * kept(near_swe, far_swe),
        near_swe.crest_to_trough,
        far_swe.crest_to_trough,
    );

    // Fixed in advance: the wave-maker radiates the requested energy (within 10%), SGN
    // keeps at least 90% of it over 400 m, and shallow water loses much more.
    assert!(
        (near.equivalent_height / 2.5 - 1.0).abs() < 0.10,
        "radiated {}",
        near.equivalent_height
    );
    assert!(
        kept(near, far) > 0.90,
        "SGN lost energy: {near:?} -> {far:?}"
    );
    assert!(
        kept(near, far) > kept(near_swe, far_swe) + 0.3,
        "no better than shallow water"
    );
}

// ---- 5. A beach: shoaling, then breaking -------------------------------------------------

/// What happened at each position on the beach over a window of several wave periods.
struct Profile {
    depth: Vec<f64>,
    crest_to_trough: Vec<f64>,
    equivalent_height: Vec<f64>,
}

/// An 8 s, 1.5 m swell from 8 m of water up a 1:30 beach (flat for 350 m, waterline at 590 m).
fn beach_profile(dx: f64) -> Profile {
    let (period, height) = (8.0, 1.5);
    let grid = Grid::new((651.0 / dx) as usize, 4, dx, dx);
    let still_depth = |x: f64| {
        if x < 350.0 {
            8.0
        } else {
            8.0 - (x - 350.0) / 30.0
        }
    };
    let bed = Bathymetry::from_fn(&grid, |x, _| -still_depth(x));
    let maker = WaveMaker {
        x: 260.0,
        amplitude: height / 2.0,
        period,
        angle: 0.0,
        sigma: 2.0 * dx,
        ramp_periods: 2.0,
    };
    let solver = Solver::new(grid, bed)
        .with_order(Order::Third)
        .with_dispersion(Dispersion::default())
        .with_wavemaker(maker)
        .with_sponge(Sponge::new((187.0 / dx) as usize, 1.5).edges(true, false, false, false));
    let mut s = State::lake_at_rest(&grid, &solver.bed, 0.0);

    let j = GHOST + 2;
    let mut seen = vec![Moments::new(); grid.nx];
    while s.time < 170.0 {
        let dt = solver.stable_dt(&s);
        solver.step(&mut s, dt);
        if s.time > 114.0 {
            for (n, m) in seen.iter_mut().enumerate() {
                let k = grid.idx(GHOST + n, j);
                if s.h[k] > 0.05 {
                    m.add(s.h[k] + solver.bed.b[k]);
                }
            }
        }
    }
    assert_eq!(solver.dispersion_stats().unwrap().unconverged, 0);
    let mut p = Profile {
        depth: Vec::new(),
        crest_to_trough: Vec::new(),
        equivalent_height: Vec::new(),
    };
    for (n, m) in seen.iter().enumerate() {
        if m.n > 0.0 {
            let x = (n as f64 + 0.5) * dx;
            p.depth.push(still_depth(x));
            p.crest_to_trough.push(m.seen().crest_to_trough);
            p.equivalent_height.push(m.seen().equivalent_height);
        }
    }
    p
}

#[test]
fn waves_shoal_intact_then_break_in_a_depth_limited_surf_zone() {
    let p = beach_profile(1.5);
    let at = |d: f64| p.depth.iter().position(|&x| x <= d).unwrap();
    let (n6, n35) = (at(6.0), at(3.5));
    let ratio = |n: usize| p.crest_to_trough[n] / p.depth[n];
    let n_max = (0..p.depth.len())
        .filter(|&n| p.depth[n] > 0.5)
        .max_by(|&a, &b| p.crest_to_trough[a].total_cmp(&p.crest_to_trough[b]))
        .unwrap();
    // The saturated surf zone, away from the swash at the waterline.
    let inner = (0..p.depth.len())
        .filter(|&n| (1.0..2.5).contains(&p.depth[n]) && n > n_max)
        .map(ratio)
        .fold(0.0, f64::max);
    let green = (p.depth[n6] / p.depth[n35]).powf(0.25);
    println!(
        "energy height at 6 m: {:.2} m (requested 1.5); crest-to-trough {:.2} m at 6 m, {:.2} m at \
         3.5 m (x{:.2}; Green's law x{green:.2}); tallest {:.2} m at {:.2} m depth (H/h = {:.2}); \
         largest H/h in 1.0-2.5 m of water past the peak: {inner:.2}",
        p.equivalent_height[n6],
        p.crest_to_trough[n6],
        p.crest_to_trough[n35],
        p.crest_to_trough[n35] / p.crest_to_trough[n6],
        p.crest_to_trough[n_max],
        p.depth[n_max],
        ratio(n_max),
    );

    // The swell arrives at 6 m depth with the energy that was asked for (within 10%),
    assert!(
        (p.equivalent_height[n6] / 1.5 - 1.0).abs() < 0.10,
        "energy height at 6 m = {}",
        p.equivalent_height[n6]
    );
    // grows as waves shoal (Green's law gives x1.14; steep waves grow somewhat more),
    let growth = p.crest_to_trough[n35] / p.crest_to_trough[n6];
    assert!((1.05..1.5).contains(&growth), "shoaling ratio {growth}");
    // breaks at a textbook breaking index (H/h between 0.55 and 1.2),
    assert!(
        (0.55..1.2).contains(&ratio(n_max)),
        "H/h at the tallest point = {}",
        ratio(n_max)
    );
    // and the surf zone is depth-limited rather than growing without bound. The first
    // version of this check looked down to 0.3 m, which is the swash, where a bore's height
    // exceeds the still depth by definition; the saturated zone is 1.0 m and deeper.
    assert!(inner < 1.0, "H/h in the surf zone reaches {inner}");
}

// ---- 6. Rest ------------------------------------------------------------------------------

#[test]
fn still_water_stays_still_with_dispersion() {
    for ((mean_depth, amplitude), order) in [(8.0, 3.0), (0.5, 2.5)]
        .into_iter()
        .flat_map(|d| ORDERS.map(|o| (d, o)))
    {
        let grid = Grid::new(64, 48, 3.0, 3.0);
        let bed = Bathymetry::rough(&grid, mean_depth, amplitude);
        let mut s = State::lake_at_rest(&grid, &bed, 0.0);
        let solver = Solver::new(grid, bed)
            .with_order(order)
            .with_dispersion(Dispersion::default());
        let v0 = s.volume(&grid);
        run_to(&solver, &mut s, 30.0);

        assert!(
            s.max_abs_momentum(&grid) < 1e-9,
            "depth {mean_depth}: still water started moving: {}",
            s.max_abs_momentum(&grid)
        );
        assert!((s.volume(&grid) - v0).abs() < 1e-9 * v0);
        assert!(s.min_depth(&grid) >= 0.0);
    }
}
