//! Row loops that run across threads when the `parallel` feature is on and the grid is big
//! enough to benefit, and one after another otherwise.
//!
//! Every loop here writes only inside its own row and reads shared data, and the reductions
//! add the per-row results in row order, so every path produces bit-identical numbers.
//! Whether a loop is threaded can therefore never change a result, only how long it takes.

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// Below this many cells a loop runs in place: handing rows to threads costs more than the
/// loop saves on a small grid.
#[cfg(feature = "parallel")]
const MIN_CELLS: usize = 40_000;

#[cfg(feature = "parallel")]
fn threaded(cells: usize) -> bool {
    cells >= MIN_CELLS
}

/// Call `f(j, row)` for every row `j` (of `width` values) of `a`.
pub(crate) fn rows1<F>(width: usize, a: &mut [f64], f: F)
where
    F: Fn(usize, &mut [f64]) + Sync + Send,
{
    #[cfg(feature = "parallel")]
    if threaded(a.len()) {
        a.par_chunks_mut(width)
            .enumerate()
            .for_each(|(j, ra)| f(j, ra));
        return;
    }
    a.chunks_mut(width).enumerate().for_each(|(j, ra)| f(j, ra));
}

/// Like [`rows1`] over two arrays of the same shape.
pub(crate) fn rows2<F>(width: usize, a: &mut [f64], b: &mut [f64], f: F)
where
    F: Fn(usize, &mut [f64], &mut [f64]) + Sync + Send,
{
    #[cfg(feature = "parallel")]
    if threaded(a.len()) {
        a.par_chunks_mut(width)
            .zip(b.par_chunks_mut(width))
            .enumerate()
            .for_each(|(j, (ra, rb))| f(j, ra, rb));
        return;
    }
    a.chunks_mut(width)
        .zip(b.chunks_mut(width))
        .enumerate()
        .for_each(|(j, (ra, rb))| f(j, ra, rb));
}

/// Like [`rows1`] over three arrays of the same shape.
pub(crate) fn rows3<F>(width: usize, a: &mut [f64], b: &mut [f64], c: &mut [f64], f: F)
where
    F: Fn(usize, &mut [f64], &mut [f64], &mut [f64]) + Sync + Send,
{
    #[cfg(feature = "parallel")]
    if threaded(a.len()) {
        a.par_chunks_mut(width)
            .zip(b.par_chunks_mut(width))
            .zip(c.par_chunks_mut(width))
            .enumerate()
            .for_each(|(j, ((ra, rb), rc))| f(j, ra, rb, rc));
        return;
    }
    a.chunks_mut(width)
        .zip(b.chunks_mut(width))
        .zip(c.chunks_mut(width))
        .enumerate()
        .for_each(|(j, ((ra, rb), rc))| f(j, ra, rb, rc));
}

/// `f(j)` for each row `j` in `j0..j1` (rows of `width` values), in row order.
fn per_row<F>(j0: usize, j1: usize, width: usize, f: F) -> Vec<f64>
where
    F: Fn(usize) -> f64 + Sync + Send,
{
    #[cfg(feature = "parallel")]
    if threaded((j1 - j0) * width) {
        return (j0..j1).into_par_iter().map(&f).collect();
    }
    #[cfg(not(feature = "parallel"))]
    let _ = width;
    (j0..j1).map(&f).collect()
}

/// The sum of `f(j)` for `j` in `j0..j1`, added in row order.
pub(crate) fn sum_rows<F>(j0: usize, j1: usize, width: usize, f: F) -> f64
where
    F: Fn(usize) -> f64 + Sync + Send,
{
    per_row(j0, j1, width, f)
        .into_iter()
        .fold(0.0, |a, b| a + b)
}

/// The largest `f(j)` for `j` in `j0..j1` (`0.0` if all are smaller).
pub(crate) fn max_rows<F>(j0: usize, j1: usize, width: usize, f: F) -> f64
where
    F: Fn(usize) -> f64 + Sync + Send,
{
    per_row(j0, j1, width, f).into_iter().fold(0.0, f64::max)
}
