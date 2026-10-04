//! Verification ladder, rung 1: still water must stay still, whatever the bed does.

use wavecore::{Bathymetry, GHOST, Grid, Order, Solver, State};

const ORDERS: [Order; 3] = [Order::First, Order::Second, Order::Fifth];

fn run(solver: &Solver, s: &mut State, steps: usize) {
    for _ in 0..steps {
        let dt = solver.stable_dt(s);
        assert!(dt.is_finite() && dt > 0.0, "no valid time step");
        solver.step(s, dt);
    }
}

#[test]
fn lake_at_rest_over_rough_wet_bed() {
    for order in ORDERS {
        let grid = Grid::new(64, 48, 3.0, 3.0);
        let bed = Bathymetry::rough(&grid, 8.0, 3.0);
        let mut s = State::lake_at_rest(&grid, &bed, 0.0);
        let solver = Solver::new(grid, bed).with_order(order);
        let v0 = s.volume(&grid);

        run(&solver, &mut s, 500);

        assert!(
            s.max_abs_momentum(&grid) < 1e-11,
            "{order:?}: still water started moving: {}",
            s.max_abs_momentum(&grid)
        );
        assert!((s.volume(&grid) - v0).abs() < 1e-9 * v0, "{order:?}");
        assert!(s.min_depth(&grid) >= 0.0, "{order:?}");
    }
}

#[test]
fn lake_at_rest_with_dry_bumps() {
    // Bed pokes 2 m above still water in places: wet/dry fronts sit next to
    // steep bed steps, the classic way this scheme fails.
    for order in ORDERS {
        let grid = Grid::new(64, 48, 3.0, 3.0);
        let bed = Bathymetry::rough(&grid, 0.5, 2.5);
        let mut s = State::lake_at_rest(&grid, &bed, 0.0);
        let solver = Solver::new(grid, bed).with_order(order);
        let v0 = s.volume(&grid);
        let dry = grid
            .interior()
            .filter(|&(i, j)| s.h[grid.idx(i, j)] == 0.0)
            .count();
        assert!(dry > 0, "test bed has no dry cells");

        run(&solver, &mut s, 500);

        assert!(
            s.max_abs_momentum(&grid) < 1e-11,
            "{order:?}: still water started moving: {}",
            s.max_abs_momentum(&grid)
        );
        assert!((s.volume(&grid) - v0).abs() < 1e-9 * v0, "{order:?}");
        assert!(s.min_depth(&grid) >= 0.0, "{order:?}");
    }
}

#[test]
fn dam_break_conserves_volume_and_stays_positive() {
    for order in ORDERS {
        let grid = Grid::new(100, 4, 1.0, 1.0);
        let bed = Bathymetry::flat(&grid, 0.0);
        let mut s = State::zeros(&grid);
        for (i, j) in grid.interior() {
            // Dry bed on the right.
            s.h[grid.idx(i, j)] = if i - GHOST < 50 { 2.0 } else { 0.0 };
        }
        let solver = Solver::new(grid, bed).with_order(order);
        let v0 = s.volume(&grid);

        run(&solver, &mut s, 300);

        assert!((s.volume(&grid) - v0).abs() < 1e-9 * v0, "{order:?}");
        assert!(s.min_depth(&grid) >= 0.0, "{order:?}");
        // Water must have advanced into the dry region.
        assert!(s.h[grid.idx(GHOST + 70, GHOST + 2)] > 0.0, "{order:?}");
    }
}
