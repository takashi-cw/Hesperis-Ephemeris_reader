# chebyshev.jl — Chebyshev 多項式評価（Clenshaw algorithm）
#
# JPL DE の各セグメントには天体位置が Chebyshev 多項式係数として
# 格納されている。このモジュールはその係数列から位置（および速度）を復元する。
#
# アルゴリズム出典:
#   - Clenshaw (1955) "A note on the summation of Chebyshev series"
#   - jplephem (Brandon Rhodes, MIT License) の設計を参考に Julia で再実装
#
# ライセンス: MIT

# --- 単軸評価 ---

"""
    chebyshev_eval(coeffs, x) -> Float64

Chebyshev 多項式を Clenshaw algorithm で評価する（位置のみ）。

f(x) = Σ c_k * T_k(x)

- `coeffs`: 係数配列 [c0, c1, ..., cn]
- `x`: 評価点（[-1, 1] に正規化済みであること）
"""
function chebyshev_eval(coeffs::Vector{Float64}, x::Float64)::Float64
    n = length(coeffs)
    n == 0 && return 0.0
    n == 1 && return coeffs[1]

    b2 = 0.0
    b1 = 0.0
    for i in n:-1:2
        b = coeffs[i] + 2.0 * x * b1 - b2
        b2 = b1
        b1 = b
    end
    return coeffs[1] + x * b1 - b2
end

"""
    chebyshev_eval_with_deriv(coeffs, x) -> (position, dpdx)

Chebyshev 多項式の位置と x に関する導関数を同時に計算する。
"""
function chebyshev_eval_with_deriv(
    coeffs::Vector{Float64},
    x::Float64,
)::Tuple{Float64,Float64}
    n = length(coeffs)
    n == 0 && return (0.0, 0.0)
    n == 1 && return (coeffs[1], 0.0)

    b2 = 0.0; b1 = 0.0
    d2 = 0.0; d1 = 0.0

    for i in n:-1:2
        b = coeffs[i] + 2.0 * x * b1 - b2
        d = 2.0 * b1 + 2.0 * x * d1 - d2
        b2 = b1; b1 = b
        d2 = d1; d1 = d
    end

    position = coeffs[1] + x * b1 - b2
    dpdx     = b1 + x * d1 - d2
    return (position, dpdx)
end

"""
    chebyshev_eval_with_velocity(coeffs, x, interval_days) -> (position_km, velocity_km_per_day)

位置（km）と速度（km/day）を計算する。

- `interval_days`: セグメントが対応する期間（日数）
"""
function chebyshev_eval_with_velocity(
    coeffs::Vector{Float64},
    x::Float64,
    interval_days::Float64,
)::Tuple{Float64,Float64}
    (position, dpdx) = chebyshev_eval_with_deriv(coeffs, x)
    velocity = dpdx * (2.0 / interval_days)
    return (position, velocity)
end

# --- 3 成分まとめて評価 ---

"""
    chebyshev_eval3(cx, cy, cz, x) -> NTuple{3,Float64}

3 成分（X, Y, Z）まとめて Chebyshev 評価する（位置のみ）。
"""
function chebyshev_eval3(
    cx::Vector{Float64},
    cy::Vector{Float64},
    cz::Vector{Float64},
    x::Float64,
)::NTuple{3,Float64}
    return (
        chebyshev_eval(cx, x),
        chebyshev_eval(cy, x),
        chebyshev_eval(cz, x),
    )
end

"""
    chebyshev_eval3_with_velocity(cx, cy, cz, x, interval_days)
        -> (position::NTuple{3,Float64}, velocity::NTuple{3,Float64})

3 成分まとめて位置（km）と速度（km/day）を計算する。
"""
function chebyshev_eval3_with_velocity(
    cx::Vector{Float64},
    cy::Vector{Float64},
    cz::Vector{Float64},
    x::Float64,
    interval_days::Float64,
)::Tuple{NTuple{3,Float64},NTuple{3,Float64}}
    (px, vx) = chebyshev_eval_with_velocity(cx, x, interval_days)
    (py, vy) = chebyshev_eval_with_velocity(cy, x, interval_days)
    (pz, vz) = chebyshev_eval_with_velocity(cz, x, interval_days)
    return ((px, py, pz), (vx, vy, vz))
end
