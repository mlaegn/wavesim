//! The potential `phi` of the flow in a closed region of the vertical plane `(x, z)`, from what is
//! known on its edge.
//!
//! The edge is made of [`Side`]s, open polylines that follow one another counter-clockwise (the
//! water on the left) and close the region. On a [`Kind::Dirichlet`] side, the free surface, the
//! potential is known and the solve finds its normal derivative `q = d phi / dn` (outward): how
//! fast the water crosses the surface, which is what moves it. On a [`Kind::Neumann`] side, the
//! bed and the walls, `q` is known (zero for an impermeable one) and the solve finds `phi`.
//!
//! The boundary integral equation, for a point `P` on the edge where the edge turns through an
//! interior angle `theta`, is
//!
//! ```text
//! (theta / 2 pi) phi(P) = integral over the edge of ( G q - phi dG/dn ) ds,   G = -ln(r) / 2 pi
//! ```
//!
//! with `r` the distance from `P`. Each side is cut into curved elements of three nodes (the
//! side's nodes 0-1-2, 2-3-4, ...), along which the position, `phi` and `q` are quadratic, and
//! the equation is applied at every node. Straight elements, with `phi` and `q` linear, were
//! tried first: they converge at second order along the surface but only at first order where
//! the surface meets a wall, because there the small error elsewhere has to fall to nothing
//! within one element; curved elements make that error one order smaller, and follow a curling
//! lip with fewer nodes.
//!
//! The integrals are Gauss-Legendre quadrature, subdivided near the point. Over an element
//! through the point they are taken on either side of it, with the distance written as a
//! multiple of the step along the element so nothing cancels, and the logarithm smoothed by a
//! change of variable. The coefficient of `phi(P)` follows from requiring that a constant
//! potential with no flow satisfy the equation, which gives the angle at corners without
//! computing it.
//!
//! Where two sides meet, each has its own node at the corner, with its own normal: `phi` is the
//! same for both, `q` is not.

use crate::linalg;

/// What is known on a side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// The potential: a free surface.
    Dirichlet,
    /// The outward normal derivative of the potential: the bed or a wall.
    Neumann,
}

/// One side of the region: its nodes in order (counter-clockwise around the region), an odd
/// number of them, what is known there, and the known value at each node.
#[derive(Clone, Debug)]
pub struct Side {
    pub kind: Kind,
    pub nodes: Vec<(f64, f64)>,
    /// `phi` at each node for a Dirichlet side, `q` for a Neumann side.
    pub known: Vec<f64>,
}

/// Both quantities at every node of every side.
#[derive(Clone, Debug)]
pub struct Solution {
    pub phi: Vec<Vec<f64>>,
    pub q: Vec<Vec<f64>>,
}

/// Gauss-Legendre points and weights on `[0, 1]`, eight of them.
pub(crate) const GAUSS: [(f64, f64); 8] = {
    const X: [f64; 4] = [
        0.183_434_642_495_649_8,
        0.525_532_409_916_329,
        0.796_666_477_413_626_7,
        0.960_289_856_497_536_3,
    ];
    const W: [f64; 4] = [
        0.362_683_783_378_362,
        0.313_706_645_877_887_3,
        0.222_381_034_453_374_5,
        0.101_228_536_290_376_3,
    ];
    let mut g = [(0.0, 0.0); 8];
    let mut k = 0;
    while k < 4 {
        g[2 * k] = (0.5 - 0.5 * X[k], 0.5 * W[k]);
        g[2 * k + 1] = (0.5 + 0.5 * X[k], 0.5 * W[k]);
        k += 1;
    }
    g
};

/// The three quadratic shape functions at `t` in `[0, 1]`, 1 at the element's first, middle and
/// last node respectively.
pub(crate) fn shape(t: f64) -> [f64; 3] {
    [
        (1.0 - t) * (1.0 - 2.0 * t),
        4.0 * t * (1.0 - t),
        t * (2.0 * t - 1.0),
    ]
}

/// A curved element through three nodes.
struct Element([(f64, f64); 3]);

impl Element {
    fn at(&self, t: f64) -> (f64, f64) {
        let n = shape(t);
        let p = &self.0;
        (
            n[0] * p[0].0 + n[1] * p[1].0 + n[2] * p[2].0,
            n[0] * p[0].1 + n[1] * p[1].1 + n[2] * p[2].1,
        )
    }

    /// `d position / dt`.
    fn tangent(&self, t: f64) -> (f64, f64) {
        let d = [4.0 * t - 3.0, 4.0 - 8.0 * t, 4.0 * t - 1.0];
        let p = &self.0;
        (
            d[0] * p[0].0 + d[1] * p[1].0 + d[2] * p[2].0,
            d[0] * p[0].1 + d[1] * p[1].1 + d[2] * p[2].1,
        )
    }

    /// Half the second derivative of the position, the same everywhere on the element.
    fn bend(&self) -> (f64, f64) {
        let p = &self.0;
        (
            2.0 * (p[0].0 - 2.0 * p[1].0 + p[2].0),
            2.0 * (p[0].1 - 2.0 * p[1].1 + p[2].1),
        )
    }
}

/// The integrals over the element `e` of `G N` and `dG/dn N`, for its three shape functions `N`,
/// seen from `p`.
fn element(p: (f64, f64), e: &Element) -> ([f64; 3], [f64; 3]) {
    let tau = std::f64::consts::TAU;
    let size = (e.0[2].0 - e.0[0].0).hypot(e.0[2].1 - e.0[0].1) + 1e-300;
    let mut g = [0.0; 3];
    let mut h = [0.0; 3];
    let own = (0..3).find(|&m| (e.0[m].0 - p.0).hypot(e.0[m].1 - p.1) <= 1e-12 * size);
    if let Some(m) = own {
        // The point is a node of the element, at tp. Moving dt = delta from it, the position
        // changes by delta (d + delta b) exactly, so r = |delta| |d + delta b| with nothing lost
        // to cancellation, and (r . n) / r^2 = (d x b) / (J |d + delta b|^2), finite.
        let tp = m as f64 / 2.0;
        let (d, b) = (e.tangent(tp), e.bend());
        let cross = d.0 * b.1 - d.1 * b.0;
        let jp = d.0.hypot(d.1);
        for length in [1.0 - tp, -tp] {
            if length == 0.0 {
                continue;
            }
            // ln|delta| times the integrand's value at the point, integrated exactly; the rest
            // vanishes there, and delta = length v^4 smooths it further for the quadrature, on
            // two halves of v.
            let l = length.abs();
            g[m] -= jp * (l * l.ln() - l) / tau;
            for half in [(0.0, 0.5), (0.5, 1.0)] {
                for &(u, w) in &GAUSS {
                    let v = half.0 + u * (half.1 - half.0);
                    let delta = length * v.powi(4);
                    let dt = l * 4.0 * v.powi(3) * w * (half.1 - half.0);
                    let step = (d.0 + delta * b.0, d.1 + delta * b.1);
                    let along = e.tangent(tp + delta);
                    let jacobian = along.0.hypot(along.1);
                    let (log_delta, log_step) = (delta.abs().ln(), step.0.hypot(step.1).ln());
                    let flux = -cross / ((step.0 * step.0 + step.1 * step.1) * tau);
                    let n = shape(tp + delta);
                    for k in 0..3 {
                        let at_point = if k == m { jp } else { 0.0 };
                        let green =
                            log_delta * (n[k] * jacobian - at_point) + log_step * n[k] * jacobian;
                        g[k] -= green / tau * dt;
                        h[k] += flux * n[k] * dt;
                    }
                }
            }
        }
        return (g, h);
    }
    // Subdivide while the point is close compared with the piece being integrated.
    // A stack on the stack: depth 12 leaves at most 13 pieces waiting.
    let mut pieces = [(0.0_f64, 1.0_f64, 0_u32); 16];
    let mut waiting = 1;
    while waiting > 0 {
        waiting -= 1;
        let (s0, s1, depth) = pieces[waiting];
        let mid = 0.5 * (s0 + s1);
        let (m, along) = (e.at(mid), e.tangent(mid));
        let distance = (m.0 - p.0).hypot(m.1 - p.1);
        if distance < 1.5 * (s1 - s0) * along.0.hypot(along.1) && depth < 12 {
            pieces[waiting] = (s0, mid, depth + 1);
            pieces[waiting + 1] = (mid, s1, depth + 1);
            waiting += 2;
            continue;
        }
        for &(u, w) in &GAUSS {
            let t = s0 + u * (s1 - s0);
            let (x, along) = (e.at(t), e.tangent(t));
            let (rx, rz) = (x.0 - p.0, x.1 - p.1);
            let r2 = rx * rx + rz * rz;
            let dt = w * (s1 - s0);
            let green = -0.5 * r2.ln() / tau * along.0.hypot(along.1);
            // The outward normal times the Jacobian: the element runs counter-clockwise, so the
            // water is on its left.
            let flux = -(rx * along.1 - rz * along.0) / (r2 * tau);
            let n = shape(t);
            for k in 0..3 {
                g[k] += green * n[k] * dt;
                h[k] += flux * n[k] * dt;
            }
        }
    }
    (g, h)
}

/// Solve for the unknown quantity at every node of `sides`, which must close a region
/// counter-clockwise, each side starting where the one before it ends. Panics where two Dirichlet
/// sides meet (the normal derivative there is not determined) or if the system is singular.
pub fn solve(sides: &[Side]) -> Solution {
    for (s, side) in sides.iter().enumerate() {
        let len = side.nodes.len();
        assert!(
            len >= 3 && len % 2 == 1,
            "side {s} has {len} nodes, not an odd number >= 3"
        );
    }
    // Every node of every side, in order: its side, its index on the side.
    let nodes: Vec<(usize, usize)> = sides
        .iter()
        .enumerate()
        .flat_map(|(s, side)| (0..side.nodes.len()).map(move |k| (s, k)))
        .collect();
    let first: Vec<usize> = sides
        .iter()
        .scan(0, |next, side| {
            let at = *next;
            *next += side.nodes.len();
            Some(at)
        })
        .collect();
    let n = nodes.len();
    let position = |i: usize| sides[nodes[i].0].nodes[nodes[i].1];
    let kind = |i: usize| sides[nodes[i].0].kind;
    let known = |i: usize| sides[nodes[i].0].known[nodes[i].1];

    // phi at every node, known or the index of an unknown: known on the surface and at a wall's
    // corner with it, shared by the two nodes where the bed meets a wall.
    let mut phi_known: Vec<Option<f64>> = (0..n)
        .map(|i| (kind(i) == Kind::Dirichlet).then(|| known(i)))
        .collect();
    let mut shares: Vec<Option<usize>> = vec![None; n];
    for s in 0..sides.len() {
        let next = (s + 1) % sides.len();
        let (last, start) = (first[s] + sides[s].nodes.len() - 1, first[next]);
        let (p, q) = (position(last), position(start));
        assert!(
            (p.0 - q.0).hypot(p.1 - q.1) < 1e-9,
            "side {s} ends at {p:?} but side {next} starts at {q:?}"
        );
        match (kind(last), kind(start)) {
            (Kind::Dirichlet, Kind::Dirichlet) => panic!("two free-surface sides meet at {p:?}"),
            (Kind::Dirichlet, Kind::Neumann) => phi_known[start] = Some(known(last)),
            (Kind::Neumann, Kind::Dirichlet) => phi_known[last] = Some(known(start)),
            (Kind::Neumann, Kind::Neumann) => shares[start] = Some(last),
        }
    }

    // Unknowns: q at a Dirichlet node, phi at a Neumann node unless known or shared. Each has the
    // equation at its own node.
    let mut phi_unknown: Vec<Option<usize>> = vec![None; n];
    let mut q_unknown: Vec<Option<usize>> = vec![None; n];
    let mut equations = Vec::new();
    for i in 0..n {
        match kind(i) {
            Kind::Dirichlet => q_unknown[i] = Some(equations.len()),
            Kind::Neumann if phi_known[i].is_none() && shares[i].is_none() => {
                phi_unknown[i] = Some(equations.len())
            }
            Kind::Neumann => continue,
        }
        equations.push(i);
    }
    for i in 0..n {
        if let Some(twin) = shares[i] {
            phi_unknown[i] = phi_unknown[twin];
        }
    }
    let count = equations.len();

    let elements: Vec<(usize, Element)> = sides
        .iter()
        .enumerate()
        .flat_map(|(s, side)| {
            let at = first[s];
            (0..side.nodes.len() / 2).map(move |e| {
                let k = 2 * e;
                (
                    at + k,
                    Element([side.nodes[k], side.nodes[k + 1], side.nodes[k + 2]]),
                )
            })
        })
        .collect();

    // sum_j H_ij phi_j - sum_j G_ij q_j = 0, with the knowns moved to the right.
    let mut a = vec![0.0; count * count];
    let mut rhs = vec![0.0; count];
    let (mut hrow, mut grow) = (vec![0.0; n], vec![0.0; n]);
    for (row, &i) in equations.iter().enumerate() {
        let p = position(i);
        hrow.fill(0.0);
        grow.fill(0.0);
        for (j, e) in &elements {
            let (g, h) = element(p, e);
            for k in 0..3 {
                grow[j + k] += g[k];
                hrow[j + k] += h[k];
            }
        }
        // A constant potential with no flow: the coefficient of phi(P) is minus the rest.
        hrow[i] -= hrow.iter().sum::<f64>();
        for j in 0..n {
            match (phi_known[j], phi_unknown[j]) {
                (Some(phi), _) => rhs[row] -= hrow[j] * phi,
                (None, Some(u)) => a[row * count + u] += hrow[j],
                (None, None) => unreachable!("node {j} has neither a known nor an unknown phi"),
            }
            match q_unknown[j] {
                Some(u) => a[row * count + u] -= grow[j],
                None => rhs[row] += grow[j] * known(j),
            }
        }
    }
    let x = linalg::solve(a, rhs).expect("the boundary element system is singular");

    let mut phi: Vec<Vec<f64>> = sides.iter().map(|s| vec![0.0; s.nodes.len()]).collect();
    let mut q = phi.clone();
    for i in 0..n {
        let (s, k) = nodes[i];
        phi[s][k] = phi_known[i].unwrap_or_else(|| x[phi_unknown[i].expect("an unknown phi")]);
        q[s][k] = q_unknown[i].map_or_else(|| known(i), |u| x[u]);
    }
    Solution { phi, q }
}
