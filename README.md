# Hesperis-ephemeris

[![MIT License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![GitHub](https://img.shields.io/badge/GitHub-takashi--cw%2FHesperis--Ephemeris__reader-black?logo=github)](https://github.com/takashi-cw/Hesperis-Ephemeris_reader)

NASA JPL DE 天体暦（`.bsp` / SPK Type 2・Type 3）を複数言語向けに実装した BSP リーダー集です。

> Multi-language readers for NASA JPL BSP (SPK Type 2 & Type 3) ephemeris files.  
> Pure implementations — no external astronomy libraries required.

---

## このプロジェクトの背景

JPL が公式配布する `.bsp` ファイルを直接読めるライブラリは、Python（[jplephem](https://github.com/brandon-rhodes/python-jplephem)）が事実上の標準実装であり、他の主要言語向けには同等のライブラリがほとんど存在しません。

```
Python  →  jplephem / Skyfield  ✅  BSP を直接読める
JS      →  —                    ❌  標準的な実装がない
Swift   →  —                    ❌  標準的な実装がない
Rust    →  —                    ❌  標準的な実装がない
Julia   →  —                    ❌  標準的な実装がない
```

`.bsp` が読めれば、座標変換・視位置計算・出没計算といった天文計算の基盤が各言語のエコシステム上で完結します。  
このリポジトリは各主要言語に BSP リーダーを提供し、**Python 以外でも本格的な天文計算ができる環境を広げること**を目的としています。

---

## 実装一覧

| ディレクトリ | 言語 | 実行環境 |
|---|---|---|
| [`js/`](./js/) | JavaScript（ES Modules） | Node.js・ブラウザ |
| [`py/`](./py/) | Rust + PyO3 | Python バインディング（高速処理） |
| [`swift/`](./swift/) | Swift（Swift Package） | macOS・iOS |
| [`rust/`](./rust/) | Rust | ネイティブ実装 |
| [`julia/`](./julia/) | Julia | ネイティブ実装 |

---

## 共通仕様

すべての実装で以下を共通とします。

| 項目 | 仕様 |
|---|---|
| 対応フォーマット | DAF/SPK Type 2・Type 3（Chebyshev 多項式） |
| エンディアン | LTL-IEEE / BIG-IEEE 両対応 |
| 出力座標系 | ICRS 直交座標 XYZ |
| 位置単位 | km |
| 速度単位 | km/day |
| 外部依存 | なし（各言語の標準ライブラリのみ） |

> **出力は ICRS（BSP 格納座標系）のまま返します。**  
> J2000.0 / of-date への座標変換は呼び出し側の責務です。

### 対応 SPK セグメントタイプ

| タイプ | 内容 | 対応状況 |
|---|---|---|
| Type 2 | Chebyshev 多項式（位置）| ✅ 全言語対応 |
| Type 3 | Chebyshev 多項式（位置＋速度）| ✅ 全言語対応（full DE440 / DE441 の月秤動角セグメント） |
| Type 13 | Hermite 補間（小天体）| ❌ スコープ外（明示的なエラーを返す） |

`de440s.bsp`（Type 2 のみ）および `de440.bsp` / `de441.bsp`（Type 3 を含む full カーネル）の両方に対応。

---

## BSP ファイルの入手

```
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
```

約 32 MB。カバー範囲: AD 1849 〜 AD 2150

---

## クイックスタート

### JavaScript（`js/`）

```js
import { loadBsp, parseBsp, NAIF } from './js/index.js';

const buffer = await loadBsp('./de440s.bsp');
const bsp    = parseBsp(buffer);
const pos    = bsp.computePosition(NAIF.SUN, NAIF.SSB, 2451545.0);
// => [x, y, z] (km, ICRS)
```

### Python（`py/`）

```python
import bsp_rs

reader = bsp_rs.BspReader("/path/to/de440s.bsp")
pos = reader.compute_position(10, 0, 2451545.0)
# => (x, y, z) km
```

### Swift（`swift/`）

```swift
import StellaBspReader
import Foundation

let bsp = try BspFile.load(url: URL(fileURLWithPath: "/path/to/de440s.bsp"))
let (x, y, z) = try bsp.computePosition(target: 10, center: 0, jdTdb: 2451545.0)
```

### Rust（`rust/`）

```rust
use stella_bsp_reader::{BspFile, naif};
use std::path::Path;

let bsp = BspFile::load(Path::new("/path/to/de440s.bsp"))?;
let pos = bsp.compute_position(naif::SUN, naif::SSB, 2451545.0)?;
// => [x, y, z] (km)
```

### Julia（`julia/`）

```julia
using StellaBspReader

bsp = load_bsp("/path/to/de440s.bsp")
pos = compute_position(bsp, Naif.SUN, Naif.SSB, 2451545.0)
# => (x, y, z) km
```

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

## ライセンス

MIT License

本リポジトリの各実装は [jplephem](https://github.com/brandon-rhodes/python-jplephem)（Brandon Rhodes, MIT License）の設計を参考に独自実装したものです。
