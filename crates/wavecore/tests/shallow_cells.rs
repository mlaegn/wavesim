//! Cells held to plain shallow water (`Solver::with_shallow_cells`), the surf-zone fallback in
//! docs/surf-zone.md: there the dispersive model must be exactly the shallow-water one, and
//! the breaking statistics must not see the difference.

use wavecore::{Bathymetry, Dispersion, GHOST, Grid, Order, Solver, State};

/// A smooth hump 10 cm tall on 2 m of still water: no breaking anywhere, so dispersion is fully
/// on unless the cells are held to shallow water.
fn hump() -> (Grid, Bathymetry, State) {
    let grid = Grid::new(120, 4, 1.0, 1.0);
    let bed = Bathymetry::flat(&grid, 2.0);
    let mut s = State::lake_at_rest(&grid, &bed, 0.0);
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            let (x, _) = grid.centre(i + GHOST, j + GHOST);
            s.h[grid.idx(i + GHOST, j + GHOST)] += 0.1 * (-((x - 40.0) / 5.0).powi(2)).exp();
        }
    }
    (grid, bed, s)
}

/// The surface after `steps` fixed steps, interior only.
fn surface(solver: &Solver, mut s: State, steps: usize) -> Vec<f64> {
    for _ in 0..steps {
        solver.step(&mut s, 0.05);
    }
    let g = solver.grid;
    g.interior().map(|(i, j)| s.h[g.idx(i, j)]).collect()
}

fn largest_gap(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0, f64::max)
}

#[test]
fn held_cells_run_exactly_as_shallow_water() {
    let (grid, bed, s) = hump();
    let dispersive = || {
        Solver::new(grid, bed.clone())
            .with_order(Order::Fifth)
            .with_dispersion(Dispersion::default())
    };
    let shallow = Solver::new(grid, bed.clone()).with_order(Order::Second);
    let all = vec![true; grid.nx * grid.ny];

    let held = surface(&dispersive().with_shallow_cells(all), s.clone(), 200);
    let plain = surface(&shallow, s.clone(), 200);
    let free = surface(&dispersive(), s, 200);

    assert!(
        largest_gap(&held, &plain) < 1e-9,
        "held cells differ from shallow water by {:e} m",
        largest_gap(&held, &plain)
    );
    // And the comparison has teeth: left alone, dispersion changes the hump by far more.
    assert!(
        largest_gap(&free, &plain) > 1e-4,
        "dispersion changed the hump by only {:e} m",
        largest_gap(&free, &plain)
    );
}

#[test]
fn only_the_held_cells_change() {
    let (grid, bed, s) = hump();
    let dispersive = || {
        Solver::new(grid, bed.clone())
            .with_order(Order::Fifth)
            .with_dispersion(Dispersion::default())
    };
    // Held from x = 110 m: the hump, starting at 40 m, gets no further than about 95 m in 10 s
    // (its front travels at 4.4 m/s), so it meets no held cell. The dispersive solve couples
    // every cell at once, so the held ones still change it a little; much less than holding
    // all of them does.
    let far: Vec<bool> = (0..grid.ny)
        .flat_map(|_| (0..grid.nx).map(|i| i >= 110))
        .collect();
    let all = vec![true; grid.nx * grid.ny];
    let held = surface(&dispersive().with_shallow_cells(far), s.clone(), 200);
    let everywhere = surface(&dispersive().with_shallow_cells(all), s.clone(), 200);
    let free = surface(&dispersive(), s, 200);
    let (near, whole) = (largest_gap(&held, &free), largest_gap(&everywhere, &free));
    assert!(
        near < 0.01 * whole,
        "held cells away from the wave changed it by {near:e} m, all held by {whole:e} m"
    );
}

#[test]
fn the_breaking_indicator_reads_the_surface_not_the_held_cells() {
    let (grid, bed, s) = hump();
    let make = || {
        Solver::new(grid, bed.clone())
            .with_order(Order::Fifth)
            .with_dispersion(Dispersion::default())
    };
    let held = make().with_shallow_cells(vec![true; grid.nx * grid.ny]);
    let indicator = held.breaking_indicator(&s);
    assert_eq!(indicator, make().breaking_indicator(&s));
    assert!(
        indicator.iter().all(|&b| b < 0.1),
        "the hump reads as breaking"
    );
}

#[test]
#[should_panic(expected = "one entry per interior cell")]
fn a_mask_of_the_wrong_size_is_refused() {
    let (grid, bed, _) = hump();
    let _ = Solver::new(grid, bed).with_shallow_cells(vec![true; 10]);
}
