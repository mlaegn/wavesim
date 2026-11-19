//! The time step and the cell that sets it, including when the water has gone wrong.

use wavecore::{Bathymetry, GHOST, Grid, Solver, State};

fn lake() -> (Solver, State) {
    let grid = Grid::new(20, 10, 2.0, 2.0);
    let bed = Bathymetry::flat(&grid, 2.0);
    let state = State::lake_at_rest(&grid, &bed, 0.0);
    (Solver::new(grid, bed), state)
}

#[test]
fn a_cell_that_is_no_longer_a_number_is_seen_not_skipped() {
    // f64::max drops a NaN, so it used to leave the step as if the cell were not there.
    let (solver, mut s) = lake();
    let g = solver.grid;
    let still = solver.stable_dt(&s);
    s.hu[g.idx(GHOST + 7, GHOST + 3)] = f64::NAN;
    assert!(solver.stable_dt(&s).is_nan(), "{}", solver.stable_dt(&s));
    let c = solver.limiting_cell(&s).unwrap();
    assert_eq!((c.i, c.j), (7, 3));
    assert!(c.speed.is_nan() && c.depth == 2.0);

    // All of it gone: NaN, not the infinite step of a domain with no water.
    let (_, mut s) = lake();
    s.h.fill(f64::NAN);
    assert!(solver.stable_dt(&s).is_nan());
    assert!(still.is_finite() && still > 0.0);
}

#[test]
fn the_limiting_cell_is_the_fastest_water() {
    let (solver, mut s) = lake();
    let g = solver.grid;
    // A film a millimetre deep at 500 m/s.
    let k = g.idx(GHOST + 12, GHOST + 4);
    s.h[k] = 1e-3;
    s.hu[k] = 0.5;
    let c = solver.limiting_cell(&s).unwrap();
    assert_eq!((c.i, c.j), (12, 4));
    assert!((c.speed - 500.0).abs() < 1e-9 && c.depth == 1e-3);
    let film = (500.0 + (9.81e-3_f64).sqrt()) / 2.0 + (9.81e-3_f64).sqrt() / 2.0;
    assert!((solver.stable_dt(&s) - solver.cfl / film).abs() < 1e-9 * solver.cfl / film);

    // Infinitely fast: a step of 0, which a run must not take forever.
    s.hu[k] = f64::INFINITY;
    assert_eq!(solver.stable_dt(&s), 0.0);
}

#[test]
fn still_water_is_limited_by_its_depth_and_dry_land_by_nothing() {
    let (solver, s) = lake();
    let c = solver.limiting_cell(&s).unwrap();
    assert_eq!((c.depth, c.speed), (2.0, 0.0));

    let grid = Grid::new(6, 4, 2.0, 2.0);
    let bed = Bathymetry::flat(&grid, -1.0);
    let dry = State::lake_at_rest(&grid, &bed, 0.0);
    let solver = Solver::new(grid, bed);
    assert!(solver.limiting_cell(&dry).is_none());
    assert_eq!(solver.stable_dt(&dry), f64::INFINITY);
}
