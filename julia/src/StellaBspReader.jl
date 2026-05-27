"""
StellaBspReader — NASA JPL BSP/SPK バイナリ天体暦ファイルの Pure Julia リーダー

対応 SPK セグメントタイプ:
  Type 2  — Chebyshev 多項式（位置のみ）        de440s.bsp / de440.bsp / de441.bsp
  Type 3  — Chebyshev 多項式（位置＋速度）       de440.bsp / de441.bsp の月秤動角セグメント
  Type 13 — Hermite 補間（小天体）               スコープ外（BspError を送出）

外部依存: なし（Julia 標準ライブラリのみ）

ライセンス: MIT
"""
module StellaBspReader

include("constants.jl")
include("chebyshev.jl")
include("bsp.jl")

export BspFile, BspSegment, BspError
export load_bsp, bsp_from_bytes
export get_position, get_position_and_velocity, compute_position
export chebyshev_eval, chebyshev_eval_with_deriv, chebyshev_eval_with_velocity
export chebyshev_eval3, chebyshev_eval3_with_velocity
export Naif
export J2000_JD, SECS_PER_DAY, AU_KM

end
