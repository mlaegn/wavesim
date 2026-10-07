//! Stage 1 of the slice solver: the potential of a flow in a region of fixed shape, against
//! flows whose potential is known exactly.

use waveslice::{Kind, Side, solve};

/// `n + 1` nodes along a straight or curved line from `t = 0` to `t = 1`: `n / 2` curved
/// elements for an even `n`.
fn line(n: usize, f: impl Fn(f64) -> (f64, f64)) -> Vec<(f64, f64)> {
    (0..=n).map(|k| f(k as f64 / n as f64)).collect()
}

/// A region bounded by a bed `z = bed(x)` below, a free surface `z = top(x)` above, and walls at
/// `x = 0` and `x = length`, with `n + 1` nodes on the bed and the surface and `m + 1` on each
/// wall, counter-clockwise: bed, right wall, surface (right to left), left wall. `phi` and its gradient
/// give the known values: the potential on the surface, and the outward normal derivative on the
/// bed and the walls.
#[allow(clippy::too_many_arguments)]
fn region(
    length: f64,
    bed: &dyn Fn(f64) -> f64,
    bed_slope: &dyn Fn(f64) -> f64,
    top: &dyn Fn(f64) -> f64,
    n: usize,
    m: usize,
    phi: &dyn Fn(f64, f64) -> f64,
    grad: &dyn Fn(f64, f64) -> (f64, f64),
) -> Vec<Side> {
    let flux = |p: (f64, f64), normal: (f64, f64)| {
        let (gx, gz) = grad(p.0, p.1);
        let l = normal.0.hypot(normal.1);
        (gx * normal.0 + gz * normal.1) / l
    };
    let bottom = line(n, |t| (t * length, bed(t * length)));
    let right = line(m, |t| {
        (length, bed(length) + t * (top(length) - bed(length)))
    });
    let surface = line(n, |t| ((1.0 - t) * length, top((1.0 - t) * length)));
    let left = line(m, |t| (0.0, top(0.0) + t * (bed(0.0) - top(0.0))));
    vec![
        Side {
            kind: Kind::Neumann,
            known: bottom
                .iter()
                .map(|&p| flux(p, (bed_slope(p.0), -1.0)))
                .collect(),
            nodes: bottom,
        },
        Side {
            kind: Kind::Neumann,
            known: right.iter().map(|&p| flux(p, (1.0, 0.0))).collect(),
            nodes: right,
        },
        Side {
            kind: Kind::Dirichlet,
            known: surface.iter().map(|&p| phi(p.0, p.1)).collect(),
            nodes: surface,
        },
        Side {
            kind: Kind::Neumann,
            known: left.iter().map(|&p| flux(p, (-1.0, 0.0))).collect(),
            nodes: left,
        },
    ]
}

#[test]
fn a_uniform_flow_comes_out_exact() {
    // phi = x + 2 z is linear, which the elements on straight sides represent exactly.
    let phi = |x: f64, z: f64| x + 2.0 * z;
    let grad = |_: f64, _: f64| (1.0, 2.0);
    let sides = region(5.0, &|_| -2.0, &|_| 0.0, &|_| 0.0, 20, 8, &phi, &grad);
    let sol = solve(&sides);
    // On the surface the normal is +z, so q = 2; on the bed and walls phi is x + 2z.
    for (k, &q) in sol.q[2].iter().enumerate() {
        assert!((q - 2.0).abs() < 1e-9, "surface node {k}: q = {q}");
    }
    for s in [0, 1, 3] {
        for (k, &(x, z)) in sides[s].nodes.iter().enumerate() {
            assert!(
                (sol.phi[s][k] - phi(x, z)).abs() < 1e-9,
                "side {s}, node {k}"
            );
        }
    }
}

/// The largest error in q on the surface, relative to the largest q there.
fn surface_error(sides: &[Side], exact: &dyn Fn(f64, f64) -> f64) -> f64 {
    let sol = solve(sides);
    let worst = sides[2]
        .nodes
        .iter()
        .zip(&sol.q[2])
        .map(|(&(x, z), &q)| (q - exact(x, z)).abs())
        .fold(0.0, f64::max);
    let scale = sides[2]
        .nodes
        .iter()
        .map(|&(x, z)| exact(x, z).abs())
        .fold(0.0, f64::max);
    worst / scale
}

#[test]
fn the_standing_wave_converges_at_second_order() {
    // A standing wave in a closed basin 10 m long and 3 m deep: phi = cosh(k(z + h)) cos(kx),
    // with no flow through the walls or the bed, and q = k sinh(kh) cos(kx) on the surface. The
    // largest error is where the surface meets the walls (1.0e-4 with 81 surface nodes, 4x less
    // for each halving); away from them it falls faster.
    let (length, h) = (10.0, 3.0);
    let k = std::f64::consts::TAU / length;
    let phi = move |x: f64, z: f64| (k * (z + h)).cosh() * (k * x).cos();
    let grad = move |x: f64, z: f64| {
        (
            -k * (k * (z + h)).cosh() * (k * x).sin(),
            k * (k * (z + h)).sinh() * (k * x).cos(),
        )
    };
    let exact = move |x: f64, _: f64| k * (k * h).sinh() * (k * x).cos();
    let errors: Vec<f64> = [20, 40, 80]
        .iter()
        .map(|&n| {
            let sides = region(
                length,
                &|_| -h,
                &|_| 0.0,
                &|_| 0.0,
                n,
                n * 3 / 10,
                &phi,
                &grad,
            );
            surface_error(&sides, &exact)
        })
        .collect();
    println!("standing wave, surface q error with 21, 41, 81 surface nodes: {errors:?}");
    assert!(errors[2] < 2e-4, "{errors:?}");
    for w in errors.windows(2) {
        assert!(
            w[0] / w[1] > 3.0,
            "halving the elements cut the error by only {}",
            w[0] / w[1]
        );
    }
}

#[test]
fn a_curved_surface_over_a_bumpy_bed_converges_too() {
    // A wavy surface over a bumpy bed, as a wave over a reef, with a flow known exactly:
    // phi = exp(0.7 z) cos(0.7 x) + 0.3 x, harmonic. The surface's slope between elements is
    // only as good as the elements' fit to the curve, which keeps this at second order (4.0e-4
    // with 97 surface nodes).
    let length = 6.0;
    let bed = |x: f64| -2.0 + 0.4 * x.sin();
    let bed_slope = |x: f64| 0.4 * x.cos();
    let top = |x: f64| 0.3 * (1.3 * x).sin();
    let top_slope = |x: f64| 0.39 * (1.3 * x).cos();
    let phi = |x: f64, z: f64| (0.7 * z).exp() * (0.7 * x).cos() + 0.3 * x;
    let grad = |x: f64, z: f64| {
        (
            -0.7 * (0.7 * z).exp() * (0.7 * x).sin() + 0.3,
            0.7 * (0.7 * z).exp() * (0.7 * x).cos(),
        )
    };
    let exact = |x: f64, z: f64| {
        let (gx, gz) = grad(x, z);
        let s = top_slope(x);
        (-s * gx + gz) / (1.0 + s * s).sqrt()
    };
    let errors: Vec<f64> = [24, 48, 96]
        .iter()
        .map(|&n| {
            let sides = region(length, &bed, &bed_slope, &top, n, n / 3, &phi, &grad);
            surface_error(&sides, &exact)
        })
        .collect();
    println!("curved surface over a bumpy bed, surface q error with 25, 49, 97 nodes: {errors:?}");
    assert!(errors[2] < 8e-4, "{errors:?}");
    for w in errors.windows(2) {
        assert!(
            w[0] / w[1] > 3.0,
            "halving the elements cut the error by only {}",
            w[0] / w[1]
        );
    }
}
