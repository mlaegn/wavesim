/// Uniform Cartesian grid in metres, stored with one ghost layer on every side.
///
/// Arrays are "padded": `(nx + 2) * (ny + 2)` values, interior cells at
/// `1..=nx` by `1..=ny`. Every kernel reads a cell and its four neighbours, so
/// each maps directly onto a GPU compute shader later.
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

    /// Number of values in a padded array.
    pub fn cells(&self) -> usize {
        (self.nx + 2) * (self.ny + 2)
    }

    /// Index into a padded array. `i` runs over `0..nx + 2`, `j` over `0..ny + 2`.
    pub fn idx(&self, i: usize, j: usize) -> usize {
        j * (self.nx + 2) + i
    }

    /// Interior cell indices `(i, j)`, each starting at 1.
    pub fn interior(&self) -> impl Iterator<Item = (usize, usize)> + use<> {
        let (nx, ny) = (self.nx, self.ny);
        (1..=ny).flat_map(move |j| (1..=nx).map(move |i| (i, j)))
    }

    /// Cell-centre coordinates in metres, origin at the lower-left corner of the domain.
    pub fn centre(&self, i: usize, j: usize) -> (f64, f64) {
        ((i as f64 - 0.5) * self.dx, (j as f64 - 0.5) * self.dy)
    }
}
