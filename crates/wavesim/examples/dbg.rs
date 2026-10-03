//! Scratch: find where the dispersive Pipeline run goes wrong. Not part of the repo.
use std::path::Path;
use wavecore::*;
use waveio::Bed;
use wavesim::{RunOptions, plan};

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let bed_path = args
        .get(1)
        .map(String::as_str)
        .unwrap_or("data/pipeline_6m.json");
    let t_from: f64 = args.get(2).and_then(|a| a.parse().ok()).unwrap_or(80.0);
    let t_to: f64 = args.get(3).and_then(|a| a.parse().ok()).unwrap_or(95.0);
    let bed = Bed::read(Path::new(bed_path)).unwrap();
    let opts = RunOptions {
        bed: bed_path.into(),
        out: "out/dbg".into(),
        height: 2.5,
        period: 14.0,
        tide: 0.0,
        duration: t_to,
        frame_interval: 2.0,
        maker_depth: 8.0,
        manning: 0.0,
        dispersive: true,
    };
    let p = plan(&bed, &opts).unwrap();
    let grid = bed.grid();
    let maker = WaveMaker {
        x: p.maker_x,
        amplitude: 1.25,
        period: 14.0,
        angle: 0.0,
        sigma: 2.0 * grid.dx,
        ramp_periods: 2.0,
    };
    let solver = Solver::new(grid, bed.bathymetry(0.0))
        .with_order(Order::Third)
        .with_dispersion(Dispersion::default())
        .with_wavemaker(maker)
        .with_sponge(Sponge::new(p.sponge_offshore, 1.5).edges(true, false, false, false))
        .with_sponge(Sponge::new(p.sponge_side, 1.5).edges(false, false, true, true));
    let mut s = State::lake_at_rest(&grid, &solver.bed, 0.0);
    let mut step = 0u64;
    let mut last_stats = DispersionStats::default();
    while s.time < t_to {
        let dt = solver.stable_dt(&s).min(t_to - s.time);
        solver.step(&mut s, dt);
        step += 1;
        let st = solver.dispersion_stats().unwrap();
        let it = (st.iterations - last_stats.iterations) as f64
            / (st.solves - last_stats.solves).max(1) as f64;
        last_stats = st;
        if s.time >= t_from && (step % 5 == 0 || it > 30.0) {
            // fastest cell
            let (mut vmax, mut at) = (0.0f64, (0, 0));
            let mut hmax = 0.0f64;
            for (i, j) in grid.interior() {
                let k = grid.idx(i, j);
                if s.h[k] > 1e-3 {
                    let sp = (s.hu[k].hypot(s.hv[k])) / s.h[k];
                    if sp > vmax {
                        vmax = sp;
                        at = (i - GHOST, j - GHOST);
                    }
                }
                hmax = hmax.max(s.h[k] + solver.bed.b[k]);
            }
            let k = grid.idx(at.0 + GHOST, at.1 + GHOST);
            println!(
                "t={:.3} dt={:.5} it/solve={:.0} max speed {:.1} m/s at cell {:?} (x={:.0} m, y={:.0} m) h={:.3} bed={:.2}; highest surface {:.2}",
                s.time,
                dt,
                it,
                vmax,
                at,
                at.0 as f64 * grid.dx,
                at.1 as f64 * grid.dy,
                s.h[k],
                solver.bed.b[k],
                hmax
            );
        }
        if !s.h.iter().all(|v| v.is_finite()) {
            println!("NaN at t={:.3}", s.time);
            break;
        }
        if dt < 1e-5 {
            println!("dt collapsed to {dt}");
            break;
        }
    }
}
