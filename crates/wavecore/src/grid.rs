/// Ghost layers on every side. Second-order reconstruction of a face needs two cells on
/// each side of it; the dispersive terms take the divergence of the velocity, then its
/// gradient, then a gradient again, which reaches three cells outward.
pub const GHOST: usize = 3;

/// Uniform Cartesian grid in metres, stored with [`GHOST`] ghost layers on every side.
///
/// Arrays are "padded": `(nx + 2 * GHOST) * (ny + 2 * GHOST)` values, interior cells
/// at `GHOST..GHOST + nx` by `GHOST..GHOST + ny`. Every kernel reads a small stencil
/// around a cell, so each maps directly onto a GPU compute shader.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Grid {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
}

impl Grid {
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Self {
        assert!(
            nx > 0 && ny > 0,
            "grid must have at least one interior cell"
        );
        assert!(dx > 0.0 && dy > 0.0, "cell size must be positive");
        Self { nx, ny, dx, dy }
    }

    /// Width of a padded row.
    pub fn width(&self) -> usize {
        self.nx + 2 * GHOST
    }

    /// Number of values in a padded array.
    pub fn cells(&self) -> usize {
        self.width() * (self.ny + 2 * GHOST)
    }

    /// Index into a padded array; `i` and `j` are padded coordinates.
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * self.width() + i
    }

    /// Interior cells as padded coordinates `(i, j)`.
    pub fn interior(&self) -> impl Iterator<Item = (usize, usize)> + use<> {
        let (nx, ny) = (self.nx, self.ny);
        (GHOST..GHOST + ny).flat_map(move |j| (GHOST..GHOST + nx).map(move |i| (i, j)))
    }

    /// Cell-centre coordinates in metres of the padded cell `(i, j)`, with the origin
    /// at the lower-left corner of the interior.
    pub fn centre(&self, i: usize, j: usize) -> (f64, f64) {
        (
            ((i - GHOST) as f64 + 0.5) * self.dx,
            ((j - GHOST) as f64 + 0.5) * self.dy,
        )
    }
}

/// Fill every ghost layer by mirroring the interior across the domain edge.
pub(crate) fn mirror_ghosts(grid: &Grid, a: &mut [f64]) {
    let (nx, ny, w) = (grid.nx, grid.ny, grid.width());
    for j in GHOST..GHOST + ny {
        for k in 0..GHOST {
            a[grid.idx(GHOST - 1 - k, j)] = a[grid.idx(GHOST + k, j)];
            a[grid.idx(GHOST + nx + k, j)] = a[grid.idx(GHOST + nx - 1 - k, j)];
        }
    }
    for i in 0..w {
        for k in 0..GHOST {
            a[grid.idx(i, GHOST - 1 - k)] = a[grid.idx(i, GHOST + k)];
            a[grid.idx(i, GHOST + ny + k)] = a[grid.idx(i, GHOST + ny - 1 - k)];
        }
    }
}

/// Fill the south and north ghost layers by wrapping around the interior, for a
/// domain that is periodic in `y`.
pub(crate) fn wrap_ghosts_y(grid: &Grid, a: &mut [f64]) {
    let (ny, w) = (grid.ny, grid.width());
    for i in 0..w {
        for k in 0..GHOST {
            a[grid.idx(i, GHOST - 1 - k)] = a[grid.idx(i, GHOST + ny - 1 - k)];
            a[grid.idx(i, GHOST + ny + k)] = a[grid.idx(i, GHOST + k)];
        }
    }
}

/// Ghost cells for a vector field `(vx, vy)` between reflective walls: both components
/// mirror, and the component normal to a wall flips sign. With `periodic_y`, the south
/// and north ghosts wrap around instead.
pub(crate) fn fill_vector_ghosts(grid: &Grid, vx: &mut [f64], vy: &mut [f64], periodic_y: bool) {
    mirror_ghosts(grid, vx);
    mirror_ghosts(grid, vy);
    let (nx, ny, w) = (grid.nx, grid.ny, grid.width());
    for j in 0..ny + 2 * GHOST {
        for k in 0..GHOST {
            let (l, r) = (grid.idx(GHOST - 1 - k, j), grid.idx(GHOST + nx + k, j));
            vx[l] = -vx[l];
            vx[r] = -vx[r];
        }
    }
    for i in 0..w {
        for k in 0..GHOST {
            let (lo, hi) = (grid.idx(i, GHOST - 1 - k), grid.idx(i, GHOST + ny + k));
            vy[lo] = -vy[lo];
            vy[hi] = -vy[hi];
        }
    }
    if periodic_y {
        wrap_ghosts_y(grid, vx);
        wrap_ghosts_y(grid, vy);
    }
}
