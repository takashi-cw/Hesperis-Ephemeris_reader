// chebyshev.rs — Chebyshev 多項式評価（Clenshaw algorithm）
//
// JPL DE の各セグメントには天体位置が Chebyshev 多項式係数として
// 格納されている。このモジュールはその係数列から位置（および速度）を復元する。
//
// アルゴリズム出典:
//   - Clenshaw (1955) "A note on the summation of Chebyshev series"
//   - jplephem (Brandon Rhodes, MIT License) の設計を参考に Rust で再実装
//
// ライセンス: MIT

// MARK: - 単軸評価

/// Chebyshev 多項式を Clenshaw algorithm で評価する（位置のみ）
///
/// f(x) = Σ c_k * T_k(x)
///
/// - `coeffs`: 係数配列 [c0, c1, ..., cn]
/// - `x`: 評価点（[-1, 1] に正規化済みであること）
pub fn chebyshev_eval(coeffs: &[f64], x: f64) -> f64 {
    let n = coeffs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return coeffs[0];
    }

    let mut b2 = 0.0f64;
    let mut b1 = 0.0f64;
    for i in (1..n).rev() {
        let b = coeffs[i] + 2.0 * x * b1 - b2;
        b2 = b1;
        b1 = b;
    }
    coeffs[0] + x * b1 - b2
}

/// Chebyshev 多項式の位置と x に関する導関数を同時に計算する
///
/// - Returns: `(position, dpdx)`
pub fn chebyshev_eval_with_deriv(coeffs: &[f64], x: f64) -> (f64, f64) {
    let n = coeffs.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    if n == 1 {
        return (coeffs[0], 0.0);
    }

    let mut b2 = 0.0f64;
    let mut b1 = 0.0f64;
    let mut d2 = 0.0f64;
    let mut d1 = 0.0f64;

    for i in (1..n).rev() {
        let b = coeffs[i] + 2.0 * x * b1 - b2;
        let d = 2.0 * b1 + 2.0 * x * d1 - d2;
        b2 = b1;
        b1 = b;
        d2 = d1;
        d1 = d;
    }

    let position = coeffs[0] + x * b1 - b2;
    let dpdx = b1 + x * d1 - d2;
    (position, dpdx)
}

/// 位置と速度（km/day）を計算する
///
/// - `interval_days`: セグメントが対応する期間（日数）
/// - Returns: `(position_km, velocity_km_per_day)`
pub fn chebyshev_eval_with_velocity(coeffs: &[f64], x: f64, interval_days: f64) -> (f64, f64) {
    let (position, dpdx) = chebyshev_eval_with_deriv(coeffs, x);
    let velocity = dpdx * (2.0 / interval_days);
    (position, velocity)
}

// MARK: - 3 成分まとめて評価

/// 3 成分（X, Y, Z）まとめて Chebyshev 評価する（位置のみ）
///
/// - `coeffs_xyz`: [coeffs_x, coeffs_y, coeffs_z]
/// - Returns: `[px, py, pz]` (km)
pub fn chebyshev_eval3(coeffs_xyz: [&[f64]; 3], x: f64) -> [f64; 3] {
    [
        chebyshev_eval(coeffs_xyz[0], x),
        chebyshev_eval(coeffs_xyz[1], x),
        chebyshev_eval(coeffs_xyz[2], x),
    ]
}

/// 3 成分まとめて位置と速度を計算する
///
/// - Returns: `(position: [px, py, pz], velocity: [vx, vy, vz])` (km, km/day)
pub fn chebyshev_eval3_with_velocity(
    coeffs_xyz: [&[f64]; 3],
    x: f64,
    interval_days: f64,
) -> ([f64; 3], [f64; 3]) {
    let (px, vx) = chebyshev_eval_with_velocity(coeffs_xyz[0], x, interval_days);
    let (py, vy) = chebyshev_eval_with_velocity(coeffs_xyz[1], x, interval_days);
    let (pz, vz) = chebyshev_eval_with_velocity(coeffs_xyz[2], x, interval_days);
    ([px, py, pz], [vx, vy, vz])
}

// MARK: - テスト

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chebyshev_eval_constant() {
        // T0(x) = 1 なので係数 [c0] の場合は c0 を返す
        assert_eq!(chebyshev_eval(&[5.0], 0.3), 5.0);
    }

    #[test]
    fn test_chebyshev_eval_linear() {
        // T0=1, T1=x なので [a, b] → a + b*x
        let coeffs = [2.0, 3.0];
        let x = 0.5;
        let expected = 2.0 + 3.0 * x;
        let result = chebyshev_eval(&coeffs, x);
        assert!((result - expected).abs() < 1e-14, "got {result}");
    }

    #[test]
    fn test_chebyshev_eval_quadratic() {
        // T2(x) = 2x^2 - 1 なので [0, 0, 1] → 2x^2 - 1
        let coeffs = [0.0, 0.0, 1.0];
        let x = 0.5;
        let expected = 2.0 * x * x - 1.0;
        let result = chebyshev_eval(&coeffs, x);
        assert!((result - expected).abs() < 1e-14, "got {result}");
    }

    #[test]
    fn test_chebyshev_eval3() {
        let cx = vec![1.0, 0.5];
        let cy = vec![2.0, 0.0];
        let cz = vec![0.0, 1.0];
        let x = 0.5;
        let pos = chebyshev_eval3([&cx, &cy, &cz], x);
        assert!((pos[0] - chebyshev_eval(&cx, x)).abs() < 1e-14);
        assert!((pos[1] - chebyshev_eval(&cy, x)).abs() < 1e-14);
        assert!((pos[2] - chebyshev_eval(&cz, x)).abs() < 1e-14);
    }
}
