//! Things that push on the water or take energy out of it: an internal wave-maker,
//! absorbing sponge layers, and Manning bottom friction.

use crate::bathymetry::Bathymetry;
use crate::grid::{GHOST, Grid};
use crate::solver::{G, H_DRY};
use crate::state::State;

/// Internal wave-maker: a Gaussian-weighted mass source along a line of constant `x`
/// (Wei et al. 1999). Waves radiate both ways from the line, so put a [`Sponge`]
/// behind it to absorb the half that travels away from the domain of interest.
///
/// The amplitude calibration assumes long waves, which is all the non-dispersive
/// equations can carry: a source of strength `Q` per unit length radiates amplitude
/// `Q / (2 c_g cos(angle))` on each side, with `c_g` the linear group speed at the source:
/// it is the speed energy leaves at, so it sets the amplitude a given mass flux produces.
/// (Using the phase speed instead overshoots by `1 + (kh)^2 / 3` for dispersive waves.)
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

/// A linear wave of a given frequency in water of a given depth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearWave {
    /// Wavenumber in rad/m.
    pub k: f64,
    /// Phase speed `omega / k` in m/s.
    pub phase_speed: f64,
    /// Group speed `d omega / d k` in m/s, the speed energy travels at.
    pub group_speed: f64,
}

/// The linear wave of angular frequency `omega` in water of depth `h`, for shallow water
/// (`dispersive = false`: `omega = k sqrt(g h)`) or for the Serre-Green-Naghdi equations
/// (`omega^2 = g h k^2 / (1 + (k h)^2 / 3)`, whose group speed is the phase speed divided
/// by `1 + (k h)^2 / 3`).
///
/// Panics if the wave is too short for the dispersive model at this depth, which has
/// no real wavenumber once `omega^2 h / g` reaches 3.
pub fn linear_wave(omega: f64, h: f64, dispersive: bool) -> LinearWave {
    if !dispersive {
        let c = (G * h).sqrt();
        return LinearWave {
            k: omega / c,
            phase_speed: c,
            group_speed: c,
        };
    }
    let nu = omega * omega * h / G;
    assert!(
        nu < 3.0,
        "a wave of angular frequency {omega} is too short for the dispersive model in {h} m of water"
    );
    let kh = (nu / (1.0 - nu / 3.0)).sqrt();
    let k = kh / h;
    let phase_speed = omega / k;
    LinearWave {
        k,
        phase_speed,
        group_speed: phase_speed / (1.0 + kh * kh / 3.0),
    }
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
    pub(crate) fn new(spec: WaveMaker, grid: &Grid, still_depth: f64, dispersive: bool) -> Self {
        assert!(still_depth > 0.0, "wave-maker must sit in water");
        assert!(spec.period > 0.0 && spec.sigma > 0.0 && spec.amplitude >= 0.0);
        let omega = std::f64::consts::TAU / spec.period;
        let LinearWave {
            k, group_speed: c, ..
        } = linear_wave(omega, still_depth, dispersive);
        let ky = k * spec.angle.sin();
        let kx = k * spec.angle.cos();
        // Total flux per unit length that gives the requested amplitude. A source of
        // Gaussian width sigma radiates a wave of wavenumber kx with its strength scaled
        // by exp(-kx^2 sigma^2 / 2), the Fourier transform of the profile, so that is
        // undone here.
        let rolloff = (0.5 * kx * kx * spec.sigma * spec.sigma).exp();
        let flux = 2.0 * c * spec.amplitude * spec.angle.cos() * rolloff;
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

/// Absorbing layer along chosen edges. Inside it, the water relaxes towards still water at
/// `level`: the depth towards the still-water depth (none on land above `level`) and the
/// momentum towards zero, at the same rate, so it never speeds water up. The rate ramps
/// quadratically from `strength` at the domain edge to nearly zero at `width` cells in. Mass
/// is not conserved there, by design.
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

/// Relax one cell towards `target = (eta, u, v)` by the factor `f` (1 keeps the cell as it is,
/// 0 replaces it): the depth towards the target's, which is none where the bed is above the
/// target's surface, and the momentum towards the target's, at the same rate. The speed then
/// always lies between the cell's and the target's, so relaxing never speeds water up.
pub(crate) fn relax_cell(
    cell: (&mut f64, &mut f64, &mut f64),
    bed: f64,
    target: (f64, f64, f64),
    f: f64,
) {
    let (h, hu, hv) = cell;
    let (eta, u, v) = target;
    let ht = (eta - bed).max(0.0);
    *h = f * *h + (1.0 - f) * ht;
    *hu = f * *hu + (1.0 - f) * ht * u;
    *hv = f * *hv + (1.0 - f) * ht * v;
}

/// Zones where the water is pulled towards a target the caller supplies at every step, such as a
/// coarser run of the same sea: this is how a fine grid is driven by a coarse one (one-way
/// nesting). The waves the target carries come in, and waves it lacks, such as reflections the
/// fine grid makes itself, are absorbed. The rates are those of [`Sponge`]s, which say where the
/// zones are; a [`Sponge`] is the special case of a target at rest.
pub struct Relaxation {
    /// Padded `(i, j)` of each cell in the zones, and its rate in 1/s.
    cells: Vec<(usize, usize, f64)>,
}

impl Relaxation {
    /// Every cell where one of `zones` has a rate, at the largest of their rates.
    pub fn new(grid: &Grid, zones: &[Sponge]) -> Self {
        let cells = grid
            .interior()
            .filter_map(|(i, j)| {
                let rate = zones.iter().map(|z| z.rate(grid, i, j)).fold(0.0, f64::max);
                (rate > 0.0).then_some((i, j, rate))
            })
            .collect();
        Self { cells }
    }

    /// The cells of the zones, padded `(i, j)`, in the order [`Relaxation::apply`] asks for
    /// their targets.
    pub fn cells(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        self.cells.iter().map(|&(i, j, _)| (i, j))
    }

    /// Relax the water in the zones over `dt` towards `target(m)`, the surface elevation and
    /// velocity `(eta, u, v)` for the `m`-th cell of [`Relaxation::cells`], in the solver's
    /// frame (still water at zero).
    pub fn apply(
        &self,
        grid: &Grid,
        bed: &Bathymetry,
        s: &mut State,
        dt: f64,
        target: impl Fn(usize) -> (f64, f64, f64),
    ) {
        for (m, &(i, j, rate)) in self.cells.iter().enumerate() {
            let k = grid.idx(i, j);
            let f = (-rate * dt).exp();
            relax_cell(
                (&mut s.h[k], &mut s.hu[k], &mut s.hv[k]),
                bed.b[k],
                target(m),
                f,
            );
            if s.h[k] < H_DRY {
                s.h[k] = s.h[k].max(0.0);
                s.hu[k] = 0.0;
                s.hv[k] = 0.0;
            }
        }
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
