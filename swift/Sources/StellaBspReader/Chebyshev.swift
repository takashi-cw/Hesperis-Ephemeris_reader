// Chebyshev.swift — Chebyshev 多項式評価（Clenshaw algorithm）
//
// Layer 1: core（依存なし）
//
// JPL DE440s の各セグメントには天体位置が Chebyshev 多項式係数として
// 格納されている。このモジュールはその係数列から位置（および速度）を復元する。
//
// アルゴリズム出典:
//   - Clenshaw (1955) "A note on the summation of Chebyshev series"
//   - Meeus "Astronomical Algorithms" 2nd ed., Ch.3
//   - jplephem (Brandon Rhodes, MIT License) の設計を参考に Swift で再実装
//
// ライセンス: MIT

// MARK: - 単軸評価

/// Chebyshev 多項式を Clenshaw algorithm で評価する（位置のみ）
///
/// f(x) = sum_{k=0}^{n} c_k * T_k(x)
///
/// - Parameters:
///   - coeffs: Chebyshev 係数配列 [c0, c1, ..., cn]
///   - x: 評価点。[-1, 1] に正規化済みであること
/// - Returns: 多項式の値
public func chebyshevEval(_ coeffs: [Double], _ x: Double) -> Double {
    let n = coeffs.count
    if n == 0 { return 0.0 }
    if n == 1 { return coeffs[0] }

    var b2 = 0.0
    var b1 = 0.0
    for i in stride(from: n - 1, through: 1, by: -1) {
        let b = coeffs[i] + 2.0 * x * b1 - b2
        b2 = b1
        b1 = b
    }
    return coeffs[0] + x * b1 - b2
}

/// Chebyshev 多項式の位置と x に関する導関数を同時に計算する
///
/// - Parameters:
///   - coeffs: Chebyshev 係数配列 [c0, c1, ..., cn]
///   - x: 評価点 [-1, 1]
/// - Returns: (position: 多項式の値, dpdx: x に関する微分値)
public func chebyshevEvalWithDeriv(_ coeffs: [Double], _ x: Double) -> (position: Double, dpdx: Double) {
    let n = coeffs.count
    if n == 0 { return (0.0, 0.0) }
    if n == 1 { return (coeffs[0], 0.0) }

    var b2 = 0.0; var b1 = 0.0
    var d2 = 0.0; var d1 = 0.0

    for i in stride(from: n - 1, through: 1, by: -1) {
        let b = coeffs[i] + 2.0 * x * b1 - b2
        let d = 2.0 * b1 + 2.0 * x * d1 - d2
        b2 = b1; b1 = b
        d2 = d1; d1 = d
    }

    let position = coeffs[0] + x * b1 - b2
    let dpdx     = b1 + x * d1 - d2
    return (position, dpdx)
}

/// 位置と速度（km/day）を計算する
///
/// - Parameters:
///   - coeffs: Chebyshev 係数配列
///   - x: 評価点 [-1, 1]
///   - intervalDays: セグメントが対応する期間（日数）
/// - Returns: (position: km, velocity: km/day)
public func chebyshevEvalWithVelocity(_ coeffs: [Double], _ x: Double, intervalDays: Double) -> (position: Double, velocity: Double) {
    let (position, dpdx) = chebyshevEvalWithDeriv(coeffs, x)
    let velocity = dpdx * (2.0 / intervalDays)
    return (position, velocity)
}

// MARK: - 3 成分まとめて評価

/// 3 成分（X, Y, Z）まとめて Chebyshev 評価する（位置のみ）
///
/// - Parameters:
///   - coeffsXYZ: [[coeffsX], [coeffsY], [coeffsZ]]
///   - x: 評価点 [-1, 1]
/// - Returns: (px, py, pz) 位置ベクトル（km）
public func chebyshevEval3(_ coeffsXYZ: [[Double]], _ x: Double) -> (Double, Double, Double) {
    return (
        chebyshevEval(coeffsXYZ[0], x),
        chebyshevEval(coeffsXYZ[1], x),
        chebyshevEval(coeffsXYZ[2], x)
    )
}

/// 3 成分まとめて位置と速度を計算する
///
/// - Parameters:
///   - coeffsXYZ: [[coeffsX], [coeffsY], [coeffsZ]]
///   - x: 評価点 [-1, 1]
///   - intervalDays: セグメント期間（日数）
/// - Returns: position: (px, py, pz)、velocity: (vx, vy, vz)（km/day）
public func chebyshevEval3WithVelocity(
    _ coeffsXYZ: [[Double]],
    _ x: Double,
    intervalDays: Double
) -> (position: (Double, Double, Double), velocity: (Double, Double, Double)) {
    let (px, vx) = chebyshevEvalWithVelocity(coeffsXYZ[0], x, intervalDays: intervalDays)
    let (py, vy) = chebyshevEvalWithVelocity(coeffsXYZ[1], x, intervalDays: intervalDays)
    let (pz, vz) = chebyshevEvalWithVelocity(coeffsXYZ[2], x, intervalDays: intervalDays)
    return (position: (px, py, pz), velocity: (vx, vy, vz))
}

// MARK: - 時刻正規化

/// 評価点を [-1, 1] に正規化する（ユーティリティ関数）
///
/// - Note: BspReader 内部では使用していない。
///   BspReader は SPK Type 2 レコードヘッダから直接 mid/radius を読み取り、
///   `x = (tSeconds - mid) / radius` でインライン正規化する（仕様準拠）。
///   本関数はセグメント境界（JD）から正規化する別経路であり、
///   テスト・デバッグ用途または将来の利用のために保持している。
///
/// - Parameters:
///   - jd: ユリウス日
///   - tStart: セグメント開始 JD
///   - tEnd: セグメント終了 JD
/// - Returns: 正規化された評価点 [-1, 1]
public func normalizeTime(jd: Double, tStart: Double, tEnd: Double) -> Double {
    return (2.0 * jd - (tStart + tEnd)) / (tEnd - tStart)
}
