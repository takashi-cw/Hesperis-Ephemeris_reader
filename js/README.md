# Hesperis-ephemeris — JS 版

JPL DE 天体暦（`.bsp` / SPK Type 2・Type 3）を Node.js またはブラウザから読むための  
**純 JavaScript ライブラリ**。ビルドツール・外部依存なし。

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
js/
├── index.js              ← エントリポイント（全 export を再エクスポート）
├── package.json
├── README.md
├── src/
│   ├── bsp-reader.js     ← BSP 読み込み本体（BspFile クラス、Type 2/3 対応）
│   ├── bsp-validator.js  ← カバー範囲検証
│   ├── chebyshev.js      ← Chebyshev 多項式評価（Clenshaw algorithm）
│   └── constants.js      ← 天文基本定数 + NAIF コード
└── test/
    └── bsp-reader.test.js ← 単体テスト（node --test）
```

---

## 動作環境

| 環境 | 要件 |
|---|---|
| Node.js | v18 以上（ES Modules `"type": "module"`） |
| ブラウザ | ES Modules 対応（Chrome 61+ / Firefox 60+ / Safari 10.1+） |
| ビルドツール | 不要 |

---

## インストール

npm 未公開の場合はパスで直接 import する：

```js
// Node.js（ファイル相対パス）
import { loadBsp, parseBsp, NAIF } from './path/to/js/index.js';

// ブラウザ（URL またはバンドラー経由）
import { loadBsp, parseBsp, NAIF } from '/js/index.js';
```

npm に公開後:

```bash
npm install hesperis-ephemeris
```

---

## クイックスタート

```js
import { loadBsp, parseBsp, NAIF } from 'hesperis-ephemeris';

// 1. BSP ファイルを読み込む
const buffer = await loadBsp('./de440s.bsp');  // Node.js: ファイルパス / ブラウザ: URL

// 2. パースして BspFile インスタンスを作成
const bsp = parseBsp(buffer);

// 3. 天体位置を取得（ICRS, km）
const J2000_TDB = 2451545.0;   // 2000-01-01 12:00 TDB

// 太陽の位置（SSB 相対）
const sunPos = bsp.computePosition(NAIF.SUN, NAIF.SSB, J2000_TDB);
console.log('太陽 (km):', sunPos);  // [x, y, z]

// 月の位置（地球相対）
const moonPos = bsp.computePosition(NAIF.MOON, NAIF.EARTH, J2000_TDB);
console.log('月 (km):', moonPos);
```

---

## API

### `loadBsp(pathOrUrl)` → `Promise<ArrayBuffer>`

`.bsp` ファイルを読み込む。Node.js では `fs.readFile`、ブラウザでは `fetch()` を使用。

### `parseBsp(buffer)` → `BspFile`

`ArrayBuffer` を解析して `BspFile` インスタンスを返す。

### `BspFile`

| メソッド / プロパティ | 説明 |
|---|---|
| `computePosition(target, center, jdTdb)` | ICRS 位置ベクトル `[x, y, z]`（km）を返す。セグメントチェーンを自動解決 |
| `getPosition(target, center, jdTdb)` | 直接セグメントが存在する場合のみ有効（セグメントチェーン不可） |
| `getPositionAndVelocity(target, center, jdTdb)` | `{ position, velocity }` を返す（velocity は km/day） |
| `segments` | パース済みセグメント記述子の配列 |
| `pairs` | `{ target, center, startJd, endJd }` の配列 |
| `name` | 内部ファイル名 |

### `getCoverageJd(bspFile)` → `{ startJd, endJd }`

BSP のカバー期間をユリウス日（JD TDB）で返す。

### `assertInCoverage(jdTdb, bspFile)`

指定 JD がカバー範囲外の場合 `RangeError` を throw する。

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

本ライブラリは [jplephem](https://github.com/brandon-rhodes/python-jplephem)（Brandon Rhodes, MIT）の設計を参考に JavaScript で独自実装したものです。