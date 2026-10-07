//! A dense linear solve: the boundary element system is a full matrix of a few hundred rows.

/// Solve `a x = b` for an `n x n` matrix `a` stored row by row, by LU factorisation with partial
/// pivoting. `None` if the matrix is singular to working precision.
pub(crate) fn solve(mut a: Vec<f64>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let n = b.len();
    assert_eq!(
        a.len(),
        n * n,
        "matrix of {} values for {n} unknowns",
        a.len()
    );
    let scale = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    for k in 0..n {
        let p = (k..n).max_by(|&i, &j| a[i * n + k].abs().total_cmp(&a[j * n + k].abs()))?;
        if a[p * n + k].abs() <= 1e-14 * scale {
            return None;
        }
        if p != k {
            for c in 0..n {
                a.swap(k * n + c, p * n + c);
            }
            b.swap(k, p);
        }
        let pivot = a[k * n + k];
        for i in k + 1..n {
            let f = a[i * n + k] / pivot;
            if f != 0.0 {
                a[i * n + k] = 0.0;
                for c in k + 1..n {
                    a[i * n + c] -= f * a[k * n + c];
                }
                b[i] -= f * b[k];
            }
        }
    }
    for k in (0..n).rev() {
        let s: f64 = (k + 1..n).map(|c| a[k * n + c] * b[c]).sum();
        b[k] = (b[k] - s) / a[k * n + k];
    }
    Some(b)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solves_a_system_with_a_known_answer_and_needs_pivoting() {
        // The first pivot is zero, so this fails without row exchanges.
        let a = vec![0.0, 2.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0];
        let x = [1.0, -2.0, 3.0];
        let b = (0..3)
            .map(|i| (0..3).map(|j| a[i * 3 + j] * x[j]).sum())
            .collect();
        let got = solve(a, b).unwrap();
        for (g, e) in got.iter().zip(x) {
            assert!((g - e).abs() < 1e-12, "{got:?}");
        }
        assert!(solve(vec![1.0, 2.0, 2.0, 4.0], vec![1.0, 2.0]).is_none());
    }
}
