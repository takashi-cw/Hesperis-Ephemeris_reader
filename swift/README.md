# Hesperis-ephemeris — Swift 版

JPL DE 天体暦（`.bsp` / SPK Type 2・Type 3）を Swift から読むための  
**純 Swift ライブラリ**。外部依存なし（`Foundation` のみ）。

---

## ファイル構成

```
swift/
├── Package.swift
├── README.md
└── Sources/StellaBspReader/
    ├── BspReader.swift      ← BSP 読み込み本体（BspFile）
    ├── BspValidator.swift   ← カバー範囲検証
    ├── Chebyshev.swift      ← Chebyshev 多項式評価（Clenshaw algorithm）
    └── Constants.swift      ← 天文基本定数 + NAIF コード
```

## 対応 SPK セグメントタイプ

| タイプ | 内容 | 対応状況 |
|---|---|---|
| Type 2 | Chebyshev 多項式（位置）| ✅ 対応 |
| Type 3 | Chebyshev 多項式（位置＋速度）| ✅ 対応（full DE440 / DE441 の月秤動角セグメント） |
| Type 13 | Hermite 補間（小天体）| ❌ スコープ外 |

`de440s.bsp`（Type 2 のみ）および `de440.bsp` / `de441.bsp`（Type 3 を含む full カーネル）の両方に対応。

---

## 動作環境

| 環境 | 要件 |
|---|---|
| macOS | 14 以上（Swift 5.9+） |
| iOS | 17 以上（Swift Package として組み込み可） |
| 外部依存 | なし（`Foundation` のみ） |

**開発・検証環境:**

| 項目 | バージョン |
|---|---|
| macOS | Tahoe 26.5 |
| Swift | 6.3.2（swiftlang-6.3.2.1.108） |
| アーキテクチャ | arm64（Apple Silicon） |

---

## インストール

Swift Package Manager で `Package.swift` に追加：

```swift
// GitHub から（推奨）
dependencies: [
    .package(url: "https://github.com/takashi-cw/Hesperis-Ephemeris_reader", from: "0.1.0")
],
targets: [
    .target(name: "MyTarget", dependencies: ["StellaBspReader"])
]
```

```swift
// ローカルパスから（開発時）
dependencies: [
    .package(path: "/path/to/Hesperis-Ephemeris_reader/swift")
],
```

---

## クイックスタート

```swift
import StellaBspReader
import Foundation

// 1. BSP ファイルを読み込む
let data = try Data(contentsOf: URL(fileURLWithPath: "/path/to/de440s.bsp"))

// 2. パースして BspFile インスタンスを作成
let bsp = try BspFile(data: data)

// 3. 天体位置を取得（ICRS, km）
let J2000_TDB = 2451545.0   // 2000-01-01 12:00 TDB

// 太陽の位置（SSB 相対）
let (x, y, z) = try bsp.computePosition(target: NAIF.SUN, center: NAIF.SSB, jdTdb: J2000_TDB)
print("太陽 (km):", x, y, z)

// 月の位置（地球相対）
let moonPos = try bsp.computePosition(target: NAIF.MOON, center: NAIF.EARTH, jdTdb: J2000_TDB)
print("月 (km):", moonPos)
```

---

## API

### `BspFile(data: Data)` throws

`Foundation.Data` を解析して `BspFile` インスタンスを返す。

### `BspFile`

| メソッド / プロパティ | 説明 |
|---|---|
| `computePosition(target:center:jdTdb:)` | ICRS 位置ベクトル `(x, y, z)`（km）を返す。セグメントチェーンを自動解決 |
| `getPosition(target:center:jdTdb:)` | 直接セグメントが存在する場合のみ有効（セグメントチェーン不可） |
| `getPositionAndVelocity(target:center:jdTdb:)` | `(pos, vel)` のタプルを返す（velocity は km/day） |
| `segments` | パース済み `BspSegment` 配列 |
| `name` | 内部ファイル名 |

### `getCoverageJd(segments:naifTarget:)` → `BspCoverage`

BSP のカバー期間をユリウス日（JD TDB）で返す。

### `assertInCoverage(jdTdb:segments:)` throws

指定 JD がカバー範囲外の場合 `BspError` を throw する。

---

## NAIF コード早見表

| 天体 | NAIF コード | 定数名 |
|---|---|---|
| 太陽系重心 (SSB) | 0 | `NAIF.SSB` |
| 太陽 | 10 | `NAIF.SUN` |
| 月 | 301 | `NAIF.MOON` |
| 地球 | 399 | `NAIF.EARTH` |
| 地球月系重心 (EMB) | 3 | `NAIF.EMB` |
| 水星重心 | 1 | `NAIF.MERCURY_BARYCENTER` |
| 金星重心 | 2 | `NAIF.VENUS_BARYCENTER` |
| 火星重心 | 4 | `NAIF.MARS_BARYCENTER` |
| 木星重心 | 5 | `NAIF.JUPITER_BARYCENTER` |
| 土星重心 | 6 | `NAIF.SATURN_BARYCENTER` |
| 天王星重心 | 7 | `NAIF.URANUS_BARYCENTER` |
| 海王星重心 | 8 | `NAIF.NEPTUNE_BARYCENTER` |
| 冥王星重心 | 9 | `NAIF.PLUTO_BARYCENTER` |

---

## DE440s BSP の入手先

```
https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp
```

約 32 MB。`AD1849 〜 AD2150` をカバー。

---

## ライセンス

MIT License

本ライブラリは [jplephem](https://github.com/brandon-rhodes/python-jplephem)（Brandon Rhodes, MIT）の設計を参考に Swift で独自実装したものです。