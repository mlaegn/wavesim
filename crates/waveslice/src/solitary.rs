//! A solitary wave to start from: one hump travelling over still water, the wave the published
//! breaking computations use, because its height alone defines it.
//!
//! This is the Serre-Green-Naghdi solitary wave, `eta = a sech^2(kappa (x - x0))` with
//! `kappa^2 = 3 a / (4 h^2 (h + a))` and speed `c^2 = g (h + a)`, and on its surface the velocity
//! of that theory's profile in depth: along, `ubar - H^2/3 ubar_xx`, up, `-H ubar_x`, with
//! `ubar = c eta / H` the depth-averaged velocity and `H = h + eta`. It is close to but not
//! exactly the potential-flow solitary wave: it sheds a small tail and settles within a few
//! depths. (With the depth-averaged velocity on the surface instead, about 4% slow under the
//! crest at a fifth of the depth, it was still swinging after ten depths.)

/// A solitary wave of `height` over still water `depth` deep, crest at `crest`.
#[derive(Clone, Copy, Debug)]
pub struct Solitary {
    pub g: f64,
    pub depth: f64,
    pub height: f64,
    pub crest: f64,
}

impl Solitary {
    fn kappa(&self) -> f64 {
        (3.0 * self.height / (4.0 * self.depth * self.depth * (self.depth + self.height))).sqrt()
    }

    /// The speed, `sqrt(g (h + a))`.
    pub fn speed(&self) -> f64 {
        (self.g * (self.depth + self.height)).sqrt()
    }

    /// The surface's height above still water at `x`.
    pub fn eta(&self, x: f64) -> f64 {
        self.height / (self.kappa() * (x - self.crest)).cosh().powi(2)
    }

    /// The potential on the surface at each of `xs`, which must increase: the integral along the
    /// surface of `u + w eta_x`, by the trapezoid rule on steps of at most a hundredth of a depth.
    pub fn potential(&self, xs: &[f64]) -> Vec<f64> {
        let (h, c) = (self.depth, self.speed());
        let ubar = |x: f64| c * self.eta(x) / (h + self.eta(x));
        let e = 1e-3 * h;
        let d1 = |f: &dyn Fn(f64) -> f64, x: f64| (f(x + e) - f(x - e)) / (2.0 * e);
        let d2 = |f: &dyn Fn(f64) -> f64, x: f64| (f(x + e) - 2.0 * f(x) + f(x - e)) / (e * e);
        let eta = |x: f64| self.eta(x);
        let slope = |x: f64| {
            let big = h + eta(x);
            let u = ubar(x) - big * big / 3.0 * d2(&ubar, x);
            let w = -big * d1(&ubar, x);
            u + w * d1(&eta, x)
        };
        let mut phi = vec![0.0; xs.len()];
        for i in 1..xs.len() {
            let (a, b) = (xs[i - 1], xs[i]);
            assert!(b > a, "the points must increase");
            let steps = ((b - a) / (0.01 * h)).ceil().max(1.0) as usize;
            let dx = (b - a) / steps as f64;
            phi[i] = phi[i - 1]
                + (0..steps)
                    .map(|k| {
                        let x = a + k as f64 * dx;
                        0.5 * (slope(x) + slope(x + dx)) * dx
                    })
                    .sum::<f64>();
        }
        phi
    }
}
