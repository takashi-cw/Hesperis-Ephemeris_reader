# Hesperis-ephemeris — Rust 版

JPL DE 天体暦（`.bsp` / SPK Type 2・Type 3）を Rust から読むための  
**純 Rust ライブラリ**。外部依存なし（標準ライブラリのみ）。

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
rust/
├── Cargo.toml
├── README.md
├── src/
│   ├── lib.rs          ← クレートルート（pub re-export）
│   ├── bsp.rs          ← BSP 読み込み本体（BspFile, BspSegment）
│   ├── chebyshev.rs    ← Chebyshev 多項式評価（Clenshaw algorithm）
│   └── constants.rs    ← 天文基本定数 + NAIF コード
└── tests/
    └── bsp_test.rs     ← 結合テスト（de440s.bsp 使用）
```

---

## 動作環境

| 環境 | 要件 |
|---|---|
| Rust | 1.65 以上（edition 2021） |
| 外部依存 | なし（`std` のみ） |

---

## インストール

`Cargo.toml` に追加：

```toml
[dependencies]
stella-bsp-reader = { path = "/path/to/Hesperis-Ephemeris_reader/rust" }
```

---

## クイックスタート

```rust
use stella_bsp_reader::{BspFile, naif};
use std::path::Path;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1. BSP ファイルを読み込む
    let bsp = BspFile::load(Path::new("/path/to/de440s.bsp"))?;

    let jd = 2451545.0; // J2000.0 (2000-01-01 12:00 TDB)

    // 2. 太陽の位置（SSB 相対, ICRS, km）
    let sun = bsp.compute_position(naif::SUN, naif::SSB, jd)?;
    println!("太陽 (km): x={:.3} y={:.3} z={:.3}", sun[0], sun[1], sun[2]);

    // 3. 地球の位置（SSB 相対）
    let earth = bsp.compute_position(naif::EARTH, naif::SSB, jd)?;
    println!("地球 (km): x={:.3} y={:.3} z={:.3}", earth[0], earth[1], earth[2]);

    // 4. 位置 + 速度
    let (pos, vel) = bsp.get_position_and_velocity(naif::SUN, naif::SSB, jd)?;
    println!("太陽速度 (km/day): vx={:.3} vy={:.3} vz={:.3}", vel[0], vel[1], vel[2]);

    Ok(())
}
```

---

## API

### `BspFile::load(path: &Path) -> Result<BspFile, BspError>`

`.bsp` ファイルを読み込んでパースする。

### `BspFile`

| メソッド | 戻り値 | 説明 |
|---|---|---|
| `compute_position(target, center, jd_tdb)` | `[f64; 3]` km | セグメントチェーンを自動解決して位置を返す |
| `get_position(target, center, jd_tdb)` | `[f64; 3]` km | 直接セグメントのみ（チェーン解決なし） |
| `get_position_and_velocity(target, center, jd_tdb)` | `([f64;3], [f64;3])` | 位置（km）+ 速度（km/day） |
| `segments` | `&[BspSegment]` | パース済みセグメント一覧 |
| `name` | `&str` | 内部ファイル名 |

### `naif` モジュール

```rust
naif::SSB    // 0  — 太陽系重心
naif::SUN    // 10 — 太陽
naif::MOON   // 301 — 月
naif::EARTH  // 399 — 地球
naif::EMB    // 3  — 地球月系重心
// … 他惑星重心は constants.rs 参照
```

---

## テスト

```bash
BSP_PATH="/path/to/de440s.bsp" cargo test -- --nocapture
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

本ライブラリは [jplephem](https://github.com/brandon-rhodes/python-jplephem)（Brandon Rhodes, MIT）の設計を参考に Rust で独自実装したものです。
