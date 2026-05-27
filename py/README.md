# Hesperis-ephemeris — Python 版（Rust + PyO3）

JPL DE 天体暦（`.bsp` / SPK Type 2・Type 3）を高速に読み込む、Rust 製ライブラリの Python バインディングです。

PyO3 を使って Rust を Python から直接呼び出す構造になっており、  
**Rust が主役・Python は呼び出し口**です。

---

## 構成

```
py/
├── Cargo.toml       # Rust クレート定義（crate-type = ["cdylib"]）
├── Cargo.lock
└── src/
    ├── lib.rs       # BSP リーダー本体 + PyO3 公開 API
    └── coord.rs     # 座標変換ロジック（ICRS ↔ 黄道、光行差補正など）
```

---

## 依存

- **Rust**: 外部クレートは `pyo3` のみ（外部天文ライブラリ依存なし）
- **Python**: `maturin`（ビルドツール）

---

## ビルド & インストール

```bash
pip install maturin
cd py
maturin develop        # 開発用（カレント環境に直接インストール）
# または
maturin build --release   # wheel ファイルを生成
```

---

## 使い方

```python
import bsp_rs

# BSP ファイルを読み込む
reader = bsp_rs.BspReader("/path/to/de440s.bsp")

# NAIF コードで天体を指定して位置を取得 [km, ICRS]
# target=399（地球）, center=10（太陽）, jd_tdb=ユリウス日（TDB）
pos = reader.compute_position(399, 10, 2451545.0)
# → [x, y, z]  単位: km

# 位置 + 速度
pos, vel = reader.compute_position_and_velocity(399, 10, 2451545.0)
# → pos [km], vel [km/day]

# バッチ計算（JD リストをまとめて渡す）
jd_list = [2451545.0 + i for i in range(365)]
positions = reader.compute_positions_batch(399, 10, jd_list)
# → [[x,y,z], ...]

# 視位置バッチ（光行時・光行差補正あり）
results = reader.compute_apparent_batch(
    naif_target=399,          # ターゲット NAIF コード
    center_naif=399,          # 観測中心 NAIF コード（399=地球）
    jd_tdb_list=jd_list,      # JD（TDB）のリスト
    use_j2000=False,          # True → J2000.0 黄道（歳差・章動なし）, False → of-date 真黄道
    aberration=True,          # True → 光偏差 + 年周光行差を適用
)
# → [(lon_deg, lat_deg, dist_km, lonspeed_deg/day, latspeed_deg/day), ...]
```

### 主な NAIF コード

| 天体 | NAIF コード |
|---|---|
| 太陽系重心 (SSB) | 0 |
| 水星バリセンター | 1 |
| 金星バリセンター | 2 |
| 地球月系重心 | 3 |
| 地球 | 399 |
| 月 | 301 |
| 火星バリセンター | 4 |
| 木星バリセンター | 5 |
| 土星バリセンター | 6 |
| 天王星バリセンター | 7 |
| 海王星バリセンター | 8 |
| 冥王星バリセンター | 9 |
| 太陽 | 10 |

---

## 公開 API

### `BspReader(bsp_path: str)`

BSP ファイルを読み込んでインスタンスを生成します。

| メソッド | 引数 | 戻り値 | 説明 |
|---|---|---|---|
| `compute_position` | `target, center, jd_tdb` | `[x,y,z]` km | 1点の ICRS 位置 |
| `compute_position_and_velocity` | `target, center, jd_tdb` | `([x,y,z], [vx,vy,vz])` | 位置＋速度 |
| `compute_positions_batch` | `target, center, jd_list` | `[[x,y,z], ...]` | 位置バッチ ※1 |
| `compute_positions_and_velocities_batch` | `target, center, jd_list` | `([[x,y,z],...], [[vx,vy,vz],...])` | 位置＋速度バッチ ※1 |
| `compute_apparent_batch` | `naif_target, center_naif, jd_list, use_j2000, aberration` | `[(lon,lat,dist,lonspd,latspd), ...]` | 視位置バッチ ※1 ※2 |
| `compute_from_center_batch` | `naif_target, center_naif, jd_list, use_j2000, aberration` | `[(lon,lat,dist,lonspd,latspd), ...]` | 任意重心視位置バッチ ※1 ※3 |
| `close()` | — | — | メモリ解放 |

> **※1 バッチ（batch）とは**
> JD（ユリウス日）のリストをまとめて1回の呼び出しで渡す方式。
> Python 側で `for` ループを回す代わりに、Rust 内のタイトループで処理するため
> Python ↔ Rust 間の呼び出しオーバーヘッドが1回で済み、大量サンプリング時に高速。
>
> **※2 `compute_apparent_batch`**
> 光行時（光が伝わる時間）・光偏差（太陽重力による光の曲がり）・年周光行差（地球公転による見かけのずれ）を補正した**視位置**を返す。
>
> **※3 `compute_from_center_batch`**
> 地球以外の任意の天体を観測中心として視位置を計算する。光偏差補正なし・年周光行差のみ適用（Skyfield の `deflectors=[]` 相当）。

---

## 対応 SPK セグメントタイプ

| タイプ | 内容 | 対応状況 |
|---|---|---|
| Type 2 | Chebyshev 多項式（位置）| ✅ 対応 |
| Type 3 | Chebyshev 多項式（位置＋速度）| ✅ 対応（full DE440 / DE441 の月秤動角セグメント） |
| Type 13 | Hermite 補間（小天体）| ❌ スコープ外（`ValueError` を送出） |

`de440s.bsp`（Type 2 のみ）および `de440.bsp` / `de441.bsp`（Type 3 を含む full カーネル）の両方に対応。

---

## ライセンス

MIT