//! Water moving in time between two vertical walls over a fixed bed.
//!
//! The free surface is a chain of nodes that move with the water, each carrying the potential
//! `phi` (a mixed Eulerian-Lagrangian method, Longuet-Higgins & Cokelet 1976). At every instant
//! the potential on the surface fixes the flow, and [`crate::solve`] gives its normal derivative
//! `q` there. With the derivative of `phi` along the surface that gives the velocity, and the
//! nodes and their potential move by
//!
//! ```text
//! dx/dt = grad phi,   dphi/dt = |grad phi|^2 / 2 - g z
//! ```
//!
//! the second being Bernoulli's equation following the water, with the still surface at `z = 0`
//! and no pressure on it. A step is the classic fourth-order Runge-Kutta method, four solves.
//! Because the nodes follow the water, the surface can fold over itself.
//!
//! The end nodes stay on the walls: the water there moves along the wall.

use crate::bem::{GAUSS, Kind, Side, shape, solve};

/// The water in a tank: a fixed bed between two vertical walls and a moving free surface.
#[derive(Clone, Debug)]
pub struct Tank {
    /// Gravity, m/s^2.
    pub g: f64,
    /// The bed from the left wall to the right one, an odd number of nodes; the walls stand at
    /// its two ends.
    pub bed: Vec<(f64, f64)>,
    /// Nodes on each wall, odd, spread evenly from the bed to the surface.
    pub wall_nodes: usize,
    /// The surface from left to right, an odd number of nodes, the first and last on the walls.
    pub surface: Vec<(f64, f64)>,
    /// The potential at each surface node.
    pub phi: Vec<f64>,
    /// Seconds since the start.
    pub time: f64,
}

/// What a moving surface is worth.
#[derive(Clone, Copy, Debug)]
pub struct Energy {
    /// Kinetic energy per unit width and density, `1/2 integral of |grad phi|^2`.
    pub kinetic: f64,
    /// Potential energy against still water, `g/2 integral of z^2 dx` along the surface.
    pub potential: f64,
}

impl Energy {
    pub fn total(&self) -> f64 {
        self.kinetic + self.potential
    }
}

impl Tank {
    /// A tank, checked: odd node counts, at least five surface nodes, the surface's ends on the
    /// walls.
    pub fn new(
        g: f64,
        bed: Vec<(f64, f64)>,
        wall_nodes: usize,
        surface: Vec<(f64, f64)>,
        phi: Vec<f64>,
    ) -> Self {
        let odd = |n: usize| n >= 3 && n % 2 == 1;
        assert!(
            odd(bed.len()),
            "{} bed nodes, not an odd number >= 3",
            bed.len()
        );
        assert!(
            odd(wall_nodes),
            "{wall_nodes} wall nodes, not an odd number >= 3"
        );
        assert!(
            odd(surface.len()) && surface.len() >= 5,
            "{} surface nodes, not an odd number >= 5",
            surface.len()
        );
        assert_eq!(phi.len(), surface.len(), "one potential per surface node");
        let (left, right) = (bed[0].0, bed[bed.len() - 1].0);
        assert!(
            surface[0].0 == left && surface[surface.len() - 1].0 == right,
            "the surface must end on the walls at x = {left} and {right}"
        );
        Tank {
            g,
            bed,
            wall_nodes,
            surface,
            phi,
            time: 0.0,
        }
    }

    /// The edge of the water counter-clockwise, as the solver takes it: bed, right wall, surface
    /// from right to left, left wall.
    fn sides(&self, surface: &[(f64, f64)], phi: &[f64]) -> Vec<Side> {
        let m = self.wall_nodes;
        let wall = |from: (f64, f64), to: f64| -> Vec<(f64, f64)> {
            (0..m)
                .map(|k| (from.0, from.1 + (to - from.1) * k as f64 / (m - 1) as f64))
                .collect()
        };
        let (left, right) = (self.bed[0], self.bed[self.bed.len() - 1]);
        let (top_left, top_right) = (surface[0], surface[surface.len() - 1]);
        let wall_left = wall((left.0, top_left.1), left.1);
        let wall_right = wall(right, top_right.1);
        let zeros = |n: usize| vec![0.0; n];
        vec![
            Side {
                kind: Kind::Neumann,
                known: zeros(self.bed.len()),
                nodes: self.bed.clone(),
            },
            Side {
                kind: Kind::Neumann,
                known: zeros(m),
                nodes: wall_right,
            },
            Side {
                kind: Kind::Dirichlet,
                nodes: surface.iter().rev().copied().collect(),
                known: phi.iter().rev().copied().collect(),
            },
            Side {
                kind: Kind::Neumann,
                known: zeros(m),
                nodes: wall_left,
            },
        ]
    }

    /// `q`, the outward normal derivative of the potential, at each surface node, left to right.
    fn flux(&self, surface: &[(f64, f64)], phi: &[f64]) -> Vec<f64> {
        let mut q = solve(&self.sides(surface, phi)).q.swap_remove(2);
        q.reverse();
        q
    }

    /// The velocity of the water at each surface node, left to right.
    pub fn velocity(&self) -> Vec<(f64, f64)> {
        velocities(
            &self.surface,
            &self.phi,
            &self.flux(&self.surface, &self.phi),
        )
    }

    /// How fast each node and its potential change.
    fn rates(&self, surface: &[(f64, f64)], phi: &[f64]) -> (Vec<(f64, f64)>, Vec<f64>) {
        let v = velocities(surface, phi, &self.flux(surface, phi));
        let dphi = v
            .iter()
            .zip(surface)
            .map(|(v, x)| 0.5 * (v.0 * v.0 + v.1 * v.1) - self.g * x.1)
            .collect();
        (v, dphi)
    }

    /// Move the water on by `dt` seconds.
    pub fn step(&mut self, dt: f64) {
        let (x0, p0) = (self.surface.clone(), self.phi.clone());
        let advance = |k: &(Vec<(f64, f64)>, Vec<f64>), h: f64| {
            let x: Vec<(f64, f64)> = x0
                .iter()
                .zip(&k.0)
                .map(|(x, v)| (x.0 + h * v.0, x.1 + h * v.1))
                .collect();
            let p: Vec<f64> = p0.iter().zip(&k.1).map(|(p, d)| p + h * d).collect();
            (x, p)
        };
        let k1 = self.rates(&x0, &p0);
        let (x, p) = advance(&k1, 0.5 * dt);
        let k2 = self.rates(&x, &p);
        let (x, p) = advance(&k2, 0.5 * dt);
        let k3 = self.rates(&x, &p);
        let (x, p) = advance(&k3, dt);
        let k4 = self.rates(&x, &p);
        let mean = |a: f64, b: f64, c: f64, d: f64| (a + 2.0 * b + 2.0 * c + d) / 6.0;
        for i in 0..x0.len() {
            self.surface[i].0 += dt * mean(k1.0[i].0, k2.0[i].0, k3.0[i].0, k4.0[i].0);
            self.surface[i].1 += dt * mean(k1.0[i].1, k2.0[i].1, k3.0[i].1, k4.0[i].1);
            self.phi[i] += dt * mean(k1.1[i], k2.1[i], k3.1[i], k4.1[i]);
        }
        self.time += dt;
    }

    /// A time step for which the fastest wave the nodes can carry, and the water itself, cross
    /// at most `courant` of the closest spacing of surface nodes. One solve.
    pub fn time_step(&self, courant: f64) -> f64 {
        let spacing = self
            .surface
            .windows(2)
            .map(|w| (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1))
            .fold(f64::INFINITY, f64::min);
        let fastest = self
            .velocity()
            .iter()
            .map(|v| v.0.hypot(v.1))
            .fold(0.0, f64::max);
        // The shortest wave the nodes carry is two spacings long, and moves at sqrt(g / k).
        let wave = (self.g * spacing / std::f64::consts::PI).sqrt();
        courant * spacing / (wave + fastest)
    }

    /// Kinetic and potential energy. One solve.
    pub fn energy(&self) -> Energy {
        let q = self.flux(&self.surface, &self.phi);
        let kinetic = 0.5
            * along_surface(&self.surface, |k, n, _, d| {
                let phi: f64 = (0..3).map(|j| n[j] * self.phi[k + j]).sum();
                let q: f64 = (0..3).map(|j| n[j] * q[k + j]).sum();
                phi * q * d.0.hypot(d.1)
            });
        let potential = 0.5 * self.g * along_surface(&self.surface, |_, _, x, d| x.1 * x.1 * d.0);
        Energy { kinetic, potential }
    }

    /// The water above still level per unit width: the integral of `z dx` along the surface.
    pub fn volume(&self) -> f64 {
        along_surface(&self.surface, |_, _, x, d| x.1 * d.0)
    }
}

/// The integral along the surface's curved elements of `f(first node, shape functions,
/// position, d position / dt)` dt.
fn along_surface(
    surface: &[(f64, f64)],
    f: impl Fn(usize, [f64; 3], (f64, f64), (f64, f64)) -> f64,
) -> f64 {
    let mut total = 0.0;
    for k in (0..surface.len() - 2).step_by(2) {
        let p = [surface[k], surface[k + 1], surface[k + 2]];
        for &(t, w) in &GAUSS {
            let n = shape(t);
            let d = [4.0 * t - 3.0, 4.0 - 8.0 * t, 4.0 * t - 1.0];
            let x = (0..3).fold((0.0, 0.0), |a, j| {
                (a.0 + n[j] * p[j].0, a.1 + n[j] * p[j].1)
            });
            let dx = (0..3).fold((0.0, 0.0), |a, j| {
                (a.0 + d[j] * p[j].0, a.1 + d[j] * p[j].1)
            });
            total += w * f(k, n, x, dx);
        }
    }
    total
}

/// The velocity at each surface node from the potential along the surface and `q` across it.
/// Along the surface, the derivative is that of the polynomial through the five nearest nodes,
/// in the distance along the chain; across, it is `q`. The end nodes move along their walls.
fn velocities(surface: &[(f64, f64)], phi: &[f64], q: &[f64]) -> Vec<(f64, f64)> {
    let n = surface.len();
    let mut s = vec![0.0; n];
    for i in 1..n {
        s[i] = s[i - 1] + (surface[i].0 - surface[i - 1].0).hypot(surface[i].1 - surface[i - 1].1);
    }
    (0..n)
        .map(|i| {
            let start = i.saturating_sub(2).min(n - 5);
            let w = derivative(&s[start..start + 5], s[i]);
            let along = |f: &dyn Fn(usize) -> f64| (0..5).map(|j| w[j] * f(start + j)).sum::<f64>();
            let (dx, dz) = (along(&|k| surface[k].0), along(&|k| surface[k].1));
            let m = dx.hypot(dz);
            let (tx, tz) = (dx / m, dz / m);
            // Outward, up out of the water, with the chain running left to right.
            let (nx, nz) = (-tz, tx);
            let tangential = along(&|k| phi[k]) / m;
            let u = tangential * tx + q[i] * nx;
            let w = tangential * tz + q[i] * nz;
            if i == 0 || i == n - 1 {
                (0.0, w)
            } else {
                (u, w)
            }
        })
        .collect()
}

/// The weights that give the derivative at `x` of the polynomial through five values at `s`.
fn derivative(s: &[f64], x: f64) -> [f64; 5] {
    let mut w = [0.0; 5];
    for j in 0..5 {
        for m in (0..5).filter(|&m| m != j) {
            let mut term = 1.0 / (s[j] - s[m]);
            for l in (0..5).filter(|&l| l != j && l != m) {
                term *= (x - s[l]) / (s[j] - s[l]);
            }
            w[j] += term;
        }
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_five_point_derivative_is_exact_for_quartics_on_uneven_spacing() {
        let s = [0.0, 0.3, 1.0, 1.2, 2.5];
        let f = |x: f64| 1.0 - 2.0 * x + 0.5 * x * x - x.powi(3) + 0.25 * x.powi(4);
        let df = |x: f64| -2.0 + x - 3.0 * x * x + x.powi(3);
        for x in [0.0, 0.3, 1.2, 2.5] {
            let w = derivative(&s, x);
            let got: f64 = (0..5).map(|j| w[j] * f(s[j])).sum();
            assert!(
                (got - df(x)).abs() < 1e-12,
                "at {x}: {got} against {}",
                df(x)
            );
        }
    }
}
