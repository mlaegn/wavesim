//! Things that push on the water or take energy out of it: an internal wave-maker,
//! absorbing sponge layers, and Manning bottom friction.

use crate::grid::{GHOST, Grid};
use crate::solver::G;

/// Internal wave-maker: a Gaussian-weighted mass source along a line of constant `x`
/// (Wei et al. 1999). Waves radiate both ways from the line, so put a [`Sponge`]
/// behind it to absorb the half that travels away from the domain of interest.
///
/// The amplitude calibration assumes long waves, which is all the non-dispersive
/// equations can carry: a source of strength `Q` per unit length radiates amplitude
/// `Q / (2 c cos(angle))` on each side, with `c = sqrt(g h)` at the source.
#[derive(Clone, Copy, Debug)]
pub struct WaveMaker {
    /// Position of the source line in metres.
    pub x: f64,
    /// Wave amplitude (half the wave height) in metres, radiated on each side.
    pub amplitude: f64,
    /// Wave period in seconds.
    pub period: f64,
    /// Propagation direction in radians, measured from +x towards +y.
    pub angle: f64,
    /// Standard deviation of the Gaussian source profile in metres. Use at least
    /// two cell widths.
    pub sigma: f64,
    /// The amplitude ramps up smoothly over this many periods.
    pub ramp_periods: f64,
}

/// A [`WaveMaker`] resolved against a grid and a still-water depth.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Drive {
    spec: WaveMaker,
    omega: f64,
    ky: f64,
    /// Peak source strength per unit area (m/s) at the centre of the Gaussian.
    strength: f64,
}

impl Drive {
    pub(crate) fn new(spec: WaveMaker, grid: &Grid, still_depth: f64) -> Self {
        assert!(still_depth > 0.0, "wave-maker must sit in water");
        assert!(spec.period > 0.0 && spec.sigma > 0.0 && spec.amplitude >= 0.0);
        let c = (G * still_depth).sqrt();
        let omega = std::f64::consts::TAU / spec.period;
        let ky = omega / c * spec.angle.sin();
        // Total flux per unit length that gives the requested amplitude...
        let flux = 2.0 * c * spec.amplitude * spec.angle.cos();
        // ...spread over the discrete Gaussian, so the amount injected is exact at any
        // resolution.
        let norm: f64 = (0..grid.nx)
            .map(|i| {
                let x = (i as f64 + 0.5) * grid.dx;
                (-(x - spec.x).powi(2) / (2.0 * spec.sigma * spec.sigma)).exp() * grid.dx
            })
            .sum();
        Self {
            spec,
            omega,
            ky,
            strength: flux / norm,
        }
    }

    /// Alongshore wavenumber imposed by the source, in rad/m.
    pub(crate) fn ky(&self) -> f64 {
        self.ky
    }

    /// Mass source in m/s at `(x, y)` and time `t`.
    pub(crate) fn source(&self, x: f64, y: f64, t: f64) -> f64 {
        let s = &self.spec;
        let profile = (-(x - s.x).powi(2) / (2.0 * s.sigma * s.sigma)).exp();
        if profile < 1e-12 {
            return 0.0;
        }
        let ramp_time = s.ramp_periods * s.period;
        let ramp = if t >= ramp_time {
            1.0
        } else {
            0.5 * (1.0 - (std::f64::consts::PI * t / ramp_time).cos())
        };
        self.strength * profile * ramp * (self.omega * t - self.ky * y).sin()
    }
}

/// Absorbing layer along chosen edges. Inside it, momentum decays and the free
/// surface relaxes to `level`; the rate ramps quadratically from `strength` at the
/// domain edge to nearly zero at `width` cells in. Mass is not conserved there, by design.
#[derive(Clone, Copy, Debug)]
pub struct Sponge {
    /// Layer thickness in cells.
    pub width: usize,
    /// Damping rate at the domain edge, in 1/s.
    pub strength: f64,
    /// Still-water surface elevation the layer relaxes towards.
    pub level: f64,
    pub west: bool,
    pub east: bool,
    pub south: bool,
    pub north: bool,
}

impl Sponge {
    /// A sponge of `width` cells and peak rate `strength` (1/s) on no edges yet.
    pub fn new(width: usize, strength: f64) -> Self {
        Self {
            width,
            strength,
            level: 0.0,
            west: false,
            east: false,
            south: false,
            north: false,
        }
    }

    /// Choose the active edges, in the order west, east, south, north.
    pub fn edges(mut self, west: bool, east: bool, south: bool, north: bool) -> Self {
        (self.west, self.east, self.south, self.north) = (west, east, south, north);
        self
    }

    /// Damping rate in 1/s at the padded interior cell `(i, j)`.
    pub(crate) fn rate(&self, grid: &Grid, i: usize, j: usize) -> f64 {
        let (a, b) = (i - GHOST, j - GHOST);
        let w = self.width as f64;
        let ramp = |d: usize| {
            if d < self.width {
                let f = (w - d as f64) / w;
                f * f
            } else {
                0.0
            }
        };
        let mut r: f64 = 0.0;
        if self.west {
            r = r.max(ramp(a));
        }
        if self.east {
            r = r.max(ramp(grid.nx - 1 - a));
        }
        if self.south {
            r = r.max(ramp(b));
        }
        if self.north {
            r = r.max(ramp(grid.ny - 1 - b));
        }
        self.strength * r
    }
}

/// Depths below this get no friction, which would otherwise blow up as `h^(-4/3)`.
pub(crate) const H_FRICTION_MIN: f64 = 1e-3;

/// Factor by which Manning friction scales the velocity over a step of `dt`, for flow
/// of the given `speed` (m/s) in depth `h` with roughness `n` (s/m^(1/3)).
///
/// It integrates `du/dt = -k u |u|` with `k = g n^2 / h^(4/3)` and the speed frozen at
/// its old value, which is the exact solution for flow in one direction and never
/// reverses the velocity.
pub fn manning_factor(h: f64, speed: f64, dt: f64, n: f64) -> f64 {
    let k = G * n * n / h.powf(4.0 / 3.0);
    1.0 / (1.0 + dt * k * speed)
}
