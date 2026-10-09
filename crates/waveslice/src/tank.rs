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
//! Because the nodes follow the water, the surface can fold over itself. Where it does, two
//! standard aids keep the nodes in order: [`Tank::smooth`] every few steps, against the
//! zig-zag that steep crests grow from node to node, and [`Tank::regrid`] where nodes crowd at a
//! lip's tip. A run ends when [`Tank::landed`]: the lip has reached the water below it.
//!
//! The end nodes stay on the walls: the water there moves along the wall.

use crate::Solitary;
use crate::bem::{GAUSS, Kind, Side, shape, solve};

/// The fraction of the closest node spacing that the fastest short wave, plus the water, may
/// cross in one step of [`Tank::march`].
pub const COURANT: f64 = 0.5;
/// [`Tank::march`] smooths the surface every this many steps.
pub const SMOOTH_EVERY: usize = 5;
/// [`Tank::march`] re-spaces nodes where neighbouring gaps differ by more than this factor.
pub const REGRID_ABOVE: f64 = 1.3;

/// How many segments apart along the surface two crossing segments must be for the crossing to
/// be a lip landing rather than a tangle of neighbouring nodes.
pub const LANDING_GAP: usize = 6;

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
    /// Steps taken by [`Tank::march`].
    pub steps: usize,
}

/// What a moving surface is worth.
#[derive(Clone, Copy, Debug)]
pub struct Energy {
    /// Kinetic energy per unit width and density, `1/2 integral of |grad phi|^2`.
    pub kinetic: f64,
    /// Potential energy against still water, `g/2 integral of z^2 dx` along the surface.
    pub potential: f64,
}

/// The tube a landed lip closes off, seen side-on.
#[derive(Clone, Copy, Debug)]
pub struct Tube {
    /// The air enclosed, m^2 per metre of crest.
    pub area: f64,
    /// Its reach from back to front, m.
    pub width: f64,
    /// Its reach from bottom to top, m.
    pub height: f64,
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
            steps: 0,
        }
    }

    /// Still water over a bed of straight pieces between `profile` points, `(x, depth)` from the
    /// left wall at the first point to the right wall at the last, with the solitary `wave` on
    /// the surface, its depth the first point's. Surface nodes are an eighth of the local depth
    /// apart, kept between `finest` and a tenth of the offshore depth, so they close up as the
    /// wave shoals. Each piece of bed is cut into curved elements about 0.6 of its least depth
    /// long, but no shorter than four finest spacings, and the profile's points are where two
    /// elements meet, so its corners are kept.
    pub fn over(wave: &Solitary, profile: &[(f64, f64)], finest: f64) -> Self {
        assert!(profile.len() >= 2, "a profile needs at least two points");
        assert!(
            profile.windows(2).all(|w| w[1].0 > w[0].0) && profile.iter().all(|p| p.1 > 0.0),
            "the profile's x must increase and its depths be positive"
        );
        let (start, end) = (profile[0].0, profile[profile.len() - 1].0);
        let bed_depth = |x: f64| {
            let k = profile
                .partition_point(|p| p.0 <= x)
                .clamp(1, profile.len() - 1);
            let (a, b) = (profile[k - 1], profile[k]);
            a.1 + (b.1 - a.1) * ((x - a.0) / (b.0 - a.0)).clamp(0.0, 1.0)
        };
        let offshore = profile[0].1;
        let mut xs = vec![start];
        while xs[xs.len() - 1] < end {
            let x = xs[xs.len() - 1];
            xs.push(x + (0.125 * bed_depth(x)).clamp(finest, 0.1 * offshore));
        }
        // Stretch to end exactly on the wall, and split the last gap for an odd number of nodes.
        let scale = (end - start) / (xs[xs.len() - 1] - start);
        xs.iter_mut()
            .for_each(|x| *x = start + (*x - start) * scale);
        if xs.len() % 2 == 0 {
            let k = xs.len() - 1;
            xs.insert(k, 0.5 * (xs[k - 1] + xs[k]));
        }
        let mut bed = vec![(start, -profile[0].1)];
        for w in profile.windows(2) {
            let (a, b) = (w[0], w[1]);
            let length = (0.6 * a.1.min(b.1)).max(4.0 * finest);
            // A hair under, so a piece that is a whole number of elements long is not given one
            // more by rounding.
            let elements = ((b.0 - a.0) / length - 1e-9).ceil().max(1.0) as usize;
            for k in 1..=2 * elements {
                let x = a.0 + (b.0 - a.0) * k as f64 / (2 * elements) as f64;
                bed.push((x, -bed_depth(x)));
            }
        }
        let phi = wave.potential(&xs);
        let surface = xs.iter().map(|&x| (x, wave.eta(x))).collect();
        Tank::new(wave.g, bed, 5, surface, phi)
    }

    /// A plane beach: [`Tank::over`] a flat bed `wave.depth` deep from `x = 0` to `toe`, then a
    /// slope rising 1 in `run` to the right wall, where the water is `shallow` deep.
    pub fn beach(wave: &Solitary, toe: f64, run: f64, shallow: f64, finest: f64) -> Self {
        let end = toe + run * (wave.depth - shallow);
        Tank::over(
            wave,
            &[(0.0, wave.depth), (toe, wave.depth), (end, shallow)],
            finest,
        )
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
        let k1 = self.rates(&self.surface, &self.phi);
        self.step_from(k1, dt);
    }

    /// Move the water on by a step for which the fastest wave the nodes can carry, and the water
    /// itself, cross at most `courant` of the closest spacing of surface nodes; the step taken.
    /// As a lip forms, its nodes crowd together and speed up, and the step shrinks with them.
    pub fn advance(&mut self, courant: f64) -> f64 {
        let k1 = self.rates(&self.surface, &self.phi);
        let dt = self.stable(&k1.0, courant);
        self.step_from(k1, dt);
        dt
    }

    /// One step the way the breaking tests take them, for a wave that may fold over: a step at
    /// [`COURANT`], the surface smoothed every [`SMOOTH_EVERY`] steps, and nodes re-spaced, up to
    /// five times, while neighbouring gaps differ by more than [`REGRID_ABOVE`]. Returns the
    /// step's length.
    pub fn march(&mut self) -> f64 {
        let dt = self.advance(COURANT);
        self.steps += 1;
        if self.steps.is_multiple_of(SMOOTH_EVERY) {
            self.smooth();
        }
        for _ in 0..5 {
            if self.unevenness() <= REGRID_ABOVE {
                break;
            }
            self.regrid();
        }
        dt
    }

    /// The step [`Tank::advance`] would take. One solve.
    pub fn time_step(&self, courant: f64) -> f64 {
        self.stable(&self.velocity(), courant)
    }

    fn stable(&self, velocity: &[(f64, f64)], courant: f64) -> f64 {
        let spacing = self
            .surface
            .windows(2)
            .map(|w| (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1))
            .fold(f64::INFINITY, f64::min);
        let fastest = velocity.iter().map(|v| v.0.hypot(v.1)).fold(0.0, f64::max);
        // The shortest wave the nodes carry is two spacings long, and moves at sqrt(g / k).
        let wave = (self.g * spacing / std::f64::consts::PI).sqrt();
        courant * spacing / (wave + fastest)
    }

    /// The classic fourth-order Runge-Kutta step, given the rates at its start.
    fn step_from(&mut self, k1: (Vec<(f64, f64)>, Vec<f64>), dt: f64) {
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

    /// Whether the surface has folded over: somewhere it runs back towards the left, so a face
    /// has passed vertical.
    pub fn overturned(&self) -> bool {
        self.surface.windows(2).any(|w| w[1].0 < w[0].0)
    }

    /// Whether a lip has landed: some node has come within one node spacing of a part of the
    /// surface at least [`LANDING_GAP`] segments away along the chain, or crossed it. From here on
    /// the water is foam and splash, which this way of computing it cannot follow, and it is also
    /// where the computing itself gives out: two parts of the surface closer than the nodes are
    /// spaced make the integrals between them nearly singular, and the next solve can throw nodes
    /// far across the tank. A run that stopped only at a crossing did exactly that, in its last
    /// four steps; when they began, the lip's tip was 0.8 spacings from the water ahead, while
    /// the two sides of the lip, the closest pair until then, stayed 1.9 spacings apart. A lip
    /// thinner than one spacing is not resolved anyway. Closer neighbours along the chain do not
    /// count: their crossing is a tangle, see [`Tank::tangled`].
    pub fn landed(&self) -> bool {
        self.landing().is_some()
    }

    /// Where a lip has landed, if it has: the node that came close, and the first node of the
    /// segment it came close to (see [`Tank::landed`]).
    pub fn landing(&self) -> Option<(usize, usize)> {
        let s = &self.surface;
        let gaps = gaps(s);
        let spacing = |i: usize| {
            let before = if i > 0 { gaps[i - 1] } else { f64::INFINITY };
            before.min(*gaps.get(i).unwrap_or(&f64::INFINITY))
        };
        let mut closest: Option<(f64, usize, usize)> = None;
        for i in 0..s.len() {
            for j in (0..s.len() - 1).filter(|&j| j + LANDING_GAP <= i || i + LANDING_GAP <= j) {
                let ratio = distance_to_segment(s[i], s[j], s[j + 1]) / spacing(i).min(gaps[j]);
                if ratio < 1.0 && closest.is_none_or(|c| ratio < c.0) {
                    closest = Some((ratio, i, j));
                }
            }
        }
        closest
            .map(|c| (c.1, c.2))
            .or_else(|| self.crossing(LANDING_GAP, usize::MAX))
    }

    /// The air a landed lip has closed off: the loop of surface from where the lip landed round
    /// to its tip. `None` before it lands.
    pub fn tube(&self) -> Option<Tube> {
        let (i, j) = self.landing()?;
        let (lo, hi) = (i.min(j + 1), i.max(j));
        let ring = &self.surface[lo..=hi];
        let area = 0.5
            * ring
                .iter()
                .zip(ring.iter().cycle().skip(1))
                .map(|(a, b)| a.0 * b.1 - b.0 * a.1)
                .sum::<f64>()
                .abs();
        let span = |f: fn(&(f64, f64)) -> f64| {
            let (lo, hi) = ring
                .iter()
                .map(f)
                .fold((f64::MAX, f64::MIN), |(l, h), v| (l.min(v), h.max(v)));
            hi - lo
        };
        Some(Tube {
            area,
            width: span(|p| p.0),
            height: span(|p| p.1),
        })
    }

    /// Whether neighbouring segments of the surface cross: the nodes have folded over each other,
    /// a numerical failure, not water.
    pub fn tangled(&self) -> bool {
        self.crossing(2, LANDING_GAP).is_some()
    }

    /// Two segments `from..until` apart along the chain that cross: the first node of each.
    fn crossing(&self, from: usize, until: usize) -> Option<(usize, usize)> {
        let s = &self.surface;
        let cross = |a: (f64, f64), b: (f64, f64), c: (f64, f64)| {
            (b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)
        };
        (0..s.len() - 1).find_map(|i| {
            (i + from..(i.saturating_add(until)).min(s.len() - 1))
                .find(|&j| {
                    let (a, b, c, d) = (s[i], s[i + 1], s[j], s[j + 1]);
                    cross(a, b, c) * cross(a, b, d) < 0.0 && cross(c, d, a) * cross(c, d, b) < 0.0
                })
                .map(|j| (i, j))
        })
    }

    /// The largest ratio between neighbouring gaps along the surface, 1 for even spacing.
    pub fn unevenness(&self) -> f64 {
        let gaps = gaps(&self.surface);
        gaps.windows(2)
            .map(|w| (w[0] / w[1]).max(w[1] / w[0]))
            .fold(1.0, f64::max)
    }

    /// Move the surface nodes along the surface, without changing its shape or the potential on
    /// it, so that each gap between nodes is the average of itself and its neighbours (weights
    /// 1-2-1): the slow grading of spacing set at the start survives, sudden jumps do not.
    /// Following the water, nodes crowd where it converges, at the tip of a lip, and the uneven
    /// spacing grows into a saw-tooth; re-spacing them is the usual cure. The new positions and
    /// potentials come from the polynomial through the five nearest nodes, in the distance along
    /// the chain, the same one the velocities use.
    pub fn regrid(&mut self) {
        let n = self.surface.len();
        let old = gaps(&self.surface);
        let mut s = vec![0.0; n];
        for i in 1..n {
            s[i] = s[i - 1] + old[i - 1];
        }
        let m = old.len();
        let mut smooth: Vec<f64> = (0..m)
            .map(|i| {
                let (before, after) = (old[i.saturating_sub(1)], old[(i + 1).min(m - 1)]);
                0.25 * (before + 2.0 * old[i] + after)
            })
            .collect();
        let scale = s[n - 1] / smooth.iter().sum::<f64>();
        smooth.iter_mut().for_each(|g| *g *= scale);
        let (mut surface, mut phi) = (self.surface.clone(), self.phi.clone());
        let mut target = 0.0;
        for i in 1..n - 1 {
            target += smooth[i - 1];
            let k = s.partition_point(|&v| v <= target).saturating_sub(1);
            let start = k.saturating_sub(2).min(n - 5);
            let w = interpolation(&s[start..start + 5], target);
            let at = |f: &dyn Fn(usize) -> f64| (0..5).map(|j| w[j] * f(start + j)).sum::<f64>();
            surface[i] = (at(&|k| self.surface[k].0), at(&|k| self.surface[k].1));
            phi[i] = at(&|k| self.phi[k]);
        }
        self.surface = surface;
        self.phi = phi;
    }

    /// Smooth the surface nodes and their potential with Longuet-Higgins & Cokelet's (1976)
    /// five-point filter, `(-f[i-2] + 4 f[i-1] + 10 f[i] + 4 f[i+1] - f[i+2]) / 16`. It removes a
    /// zig-zag from node to node entirely and changes a smooth surface only by its fourth
    /// derivative times the spacing to the fourth: a wave 40 nodes long keeps 99.996% of its
    /// height. Steep crests grow such zig-zags when nodes follow the water,
    /// whatever the resolution; smoothing every few steps keeps them down. The two nodes nearest
    /// each wall are left alone.
    pub fn smooth(&mut self) {
        let n = self.surface.len();
        let filter = |f: &dyn Fn(usize) -> f64, i: usize| {
            (-f(i - 2) + 4.0 * f(i - 1) + 10.0 * f(i) + 4.0 * f(i + 1) - f(i + 2)) / 16.0
        };
        let (x, p) = (self.surface.clone(), self.phi.clone());
        for i in 2..n - 2 {
            self.surface[i] = (filter(&|k| x[k].0, i), filter(&|k| x[k].1, i));
            self.phi[i] = filter(&|k| p[k], i);
        }
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

/// The distance from `p` to the segment from `a` to `b`.
fn distance_to_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (dx, dz) = (b.0 - a.0, b.1 - a.1);
    let along = (((p.0 - a.0) * dx + (p.1 - a.1) * dz) / (dx * dx + dz * dz)).clamp(0.0, 1.0);
    (a.0 + along * dx - p.0).hypot(a.1 + along * dz - p.1)
}

/// The distance between each pair of neighbouring nodes.
fn gaps(surface: &[(f64, f64)]) -> Vec<f64> {
    surface
        .windows(2)
        .map(|w| (w[1].0 - w[0].0).hypot(w[1].1 - w[0].1))
        .collect()
}

/// The weights that give the value at `x` of the polynomial through five values at `s`.
fn interpolation(s: &[f64], x: f64) -> [f64; 5] {
    let mut w = [1.0; 5];
    for j in 0..5 {
        for l in (0..5).filter(|&l| l != j) {
            w[j] *= (x - s[l]) / (s[j] - s[l]);
        }
    }
    w
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
    fn a_beach_keeps_the_toe_and_the_spacing_it_was_tested_with() {
        let wave = Solitary {
            g: 9.81,
            depth: 1.0,
            height: 0.3,
            crest: 7.5,
        };
        let t = Tank::beach(&wave, 15.0, 15.0, 0.05, 0.05);
        // 25 elements on the flat, 72 on the slope, 0.2 long; the toe is where two meet.
        assert_eq!(t.bed.len(), 2 * (25 + 72) + 1);
        assert_eq!(t.bed[50], (15.0, -1.0));
        assert!((t.bed[t.bed.len() - 1].1 + 0.05).abs() < 1e-12);
        let last = t.surface[t.surface.len() - 1].0;
        assert!((last - 29.25).abs() < 1e-9, "the surface ends at {last}");
    }

    #[test]
    fn a_closed_loop_is_a_landed_lip_and_its_tube_is_measured() {
        // Flat water with a circle of radius 0.5 drawn into it, from its bottom round to a tip
        // that comes back down within a hair of the flat water ahead.
        let mut surface: Vec<(f64, f64)> = (0..=20).map(|k| (0.1 * k as f64, 0.0)).collect();
        let points = 80;
        for k in 1..points {
            let a = std::f64::consts::TAU * k as f64 / points as f64;
            surface.push((
                2.0 + 0.5 * a.sin(),
                0.5 - 0.5 * a.cos() + 0.03 * k as f64 / points as f64,
            ));
        }
        let tip = surface.len() - 1;
        surface.extend((0..=19).map(|k| (2.05 + 0.1 * k as f64, -0.01)));
        surface.push((4.1, -0.01));
        if surface.len().is_multiple_of(2) {
            surface.push((4.2, -0.01));
        }
        let end = surface[surface.len() - 1].0;
        let phi = vec![0.0; surface.len()];
        let t = Tank::new(
            9.81,
            vec![(0.0, -1.0), (end / 2.0, -1.0), (end, -1.0)],
            3,
            surface,
            phi,
        );
        // The landing joins the two ends of the loop, the tip and where the loop began (node 20),
        // in either order.
        let (node, segment) = t.landing().expect("the loop has closed");
        let ends = (node.min(segment), node.max(segment));
        assert!(
            ends.0.abs_diff(20) <= 3 && ends.1.abs_diff(tip) <= 3,
            "landed between {ends:?}; the loop runs from 20 to {tip}"
        );
        let tube = t.tube().unwrap();
        let circle = std::f64::consts::PI * 0.25;
        assert!(
            (tube.area / circle - 1.0).abs() < 0.05,
            "area {} against {circle}",
            tube.area
        );
        assert!(
            (tube.width - 1.0).abs() < 0.05 && (tube.height - 1.0).abs() < 0.08,
            "{tube:?}"
        );
    }

    #[test]
    fn five_point_interpolation_is_exact_for_quartics_on_uneven_spacing() {
        let s = [0.0, 0.3, 1.0, 1.2, 2.5];
        let f = |x: f64| 1.0 - 2.0 * x + 0.5 * x * x - x.powi(3) + 0.25 * x.powi(4);
        for x in [0.1, 0.7, 1.1, 2.0] {
            let w = interpolation(&s, x);
            let got: f64 = (0..5).map(|j| w[j] * f(s[j])).sum();
            assert!((got - f(x)).abs() < 1e-12, "at {x}: {got} against {}", f(x));
        }
    }

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
