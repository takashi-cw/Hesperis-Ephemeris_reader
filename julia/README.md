# Hesperis-ephemeris — Julia 版

JPL DE 天体暦（`.bsp` / SPK Type 2・Type 3）を Julia から読むための  
**純 Julia ライブラリ**。外部依存なし（標準ライブラリのみ）。

---

## 対応 SPK セグメントタイプ

| タイプ | 内容 | 対応状況 |
|---|---|---|
| Type 2 | Chebyshev 多項式（位置）| ✅ 対応 |
| Type 3 | Chebyshev 多項式（位置＋速度）| ✅ 対応（full DE440 / DE441 の月秤動角セグメント） |
| Type 13 | Hermite 補間（小天体）| ❌ スコープ外 |

`de440s.bsp`（Type 2 のみ）および `de440.bsp` / `de441.bsp`（Type 3 を含む full カーネル）の両方に対応。

---

## ファイル構成

```
julia/
├── Project.toml
├── README.md
├── src/
│   ├── StellaBspReader.jl   ← モジュールルート（include + export）
│   ├── bsp.jl               ← BSP 読み込み本体（BspFile, BspSegment）
│   ├── chebyshev.jl         ← Chebyshev 多項式評価（Clenshaw algorithm）
│   └── constants.jl         ← 天文基本定数 + NAIF コード
└── test/
    └── runtests.jl          ← 単体テスト（合成 BSP）+ 結合テスト（de440s.bsp）
```

---

## 動作環境

| 環境 | 要件 |
|---|---|
| Julia | 1.9 以上 |
| 外部依存 | なし（標準ライブラリのみ） |

---

## インストール

Julia REPL で：

```julia
] add /path/to/Hesperis-Ephemeris_reader/julia
```

または `Project.toml` に `develop` で追加：

```julia
using Pkg
Pkg.develop(path="/path/to/Hesperis-Ephemeris_reader/julia")
```

---

## クイックスタート

```julia
using StellaBspReader

# 1. BSP ファイルを読み込む
bsp = load_bsp("/path/to/de440s.bsp")

jd = 2451545.0  # J2000.0 (2000-01-01 12:00 TDB)

# 2. 太陽の位置（SSB 相対, ICRS, km）
sun = compute_position(bsp, Naif.SUN, Naif.SSB, jd)
println("太陽 (km): x=$(sun[1]) y=$(sun[2]) z=$(sun[3])")

# 3. 地球の位置（SSB 相対）
earth = compute_position(bsp, Naif.EARTH, Naif.SSB, jd)
println("地球 (km): x=$(earth[1]) y=$(earth[2]) z=$(earth[3])")

# 4. 位置 + 速度
pos, vel = get_position_and_velocity(bsp, Naif.SUN, Naif.SSB, jd)
println("太陽速度 (km/day): vx=$(vel[1]) vy=$(vel[2]) vz=$(vel[3])")
```

---

## API

### `load_bsp(path::String) -> BspFile`

`.bsp` ファイルを読み込んでパースする。

### `BspFile` に対する関数

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `compute_position(bsp, target, center, jd_tdb)` | `NTuple{3,Float64}` km | セグメントチェーンを自動解決して位置を返す |
| `get_position(bsp, target, center, jd_tdb)` | `NTuple{3,Float64}` km | 直接セグメントのみ（チェーン解決なし） |
| `get_position_and_velocity(bsp, target, center, jd_tdb)` | `(NTuple{3}, NTuple{3})` | 位置（km）+ 速度（km/day） |
| `bsp.segments` | `Vector{BspSegment}` | パース済みセグメント一覧 |
| `bsp.name` | `String` | 内部ファイル名 |

### `Naif` モジュール

```julia
Naif.SSB    # Int32(0)   — 太陽系重心
Naif.SUN    # Int32(10)  — 太陽
Naif.MOON   # Int32(301) — 月
Naif.EARTH  # Int32(399) — 地球
Naif.EMB    # Int32(3)   — 地球月系重心
# … 他惑星重心は constants.jl 参照
```

---

## テスト

```bash
BSP_PATH="/path/to/de440s.bsp" julia --project=. test/runtests.jl
```

`BSP_PATH` が未設定の場合、BSP 結合テストはスキップされます（Chebyshev 単体テストは常時実行）。

---

## NAIF コード早見表

| 天体 | NAIF コード |
|---|---|
| 太陽系重心 (SSB) | 0 |
| 太陽 | 10 |
| 月 | 301 |
| 地球 | 399 |
| 地球月系重心 (EMB) | 3 |
| 水星重心 | 1 |
| 金星重心 | 2 |
| 火星重心 | 4 |
| 木星重心 | 5 |
| 土星重心 | 6 |
| 天王星重心 | 7 |
| 海王星重心 | 8 |
| 冥王星重心 | 9 |

---

## DE440s BSP の入手先

```
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
```

約 32 MB。カバー範囲: AD 1849 〜 AD 2150

---

## ライセンス

MIT License

本ライブラリは [jplephem](https://github.com/brandon-rhodes/python-jplephem)（Brandon Rhodes, MIT）の設計を参考に Julia で独自実装したものです。
