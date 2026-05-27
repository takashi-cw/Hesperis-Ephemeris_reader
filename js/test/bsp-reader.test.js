/**
 * bsp-reader.test.js — bsp-reader.js の単体テスト
 *
 * 実行: node --test test/bsp-reader.test.js
 *
 * テスト構成:
 *   1. 単体テスト（合成モックデータ）— .bsp ファイルなしで実行可能
 *   2. 結合テスト（実データ）        — de440s.bsp が存在する場合のみ実行
 */

import { strict as assert } from 'node:assert';
import { describe, it, before } from 'node:test';
import { existsSync } from 'node:fs';
import { fileURLToPath } from 'node:url';

import { loadBsp, parseBsp, BspFile } from '../src/bsp-reader.js';
import { AU_KM, NAIF, J2000_JD } from '../src/constants.js';

// de440s.bsp のパス（test/ からの相対パス）
const BSP_PATH = '../de440s.bsp';
const BSP_ABS  = fileURLToPath(new URL(BSP_PATH, import.meta.url));
const HAS_BSP  = existsSync(BSP_ABS);

// =========================================================================
// 1. 単体テスト（API・型の確認）
// =========================================================================
describe('bsp-reader — exports の確認', () => {
  it('loadBsp は関数', () => {
    assert.strictEqual(typeof loadBsp, 'function');
  });

  it('parseBsp は関数', () => {
    assert.strictEqual(typeof parseBsp, 'function');
  });

  it('BspFile はクラス', () => {
    assert.strictEqual(typeof BspFile, 'function');
  });
});

describe('bsp-reader — 不正な入力に対するエラー', () => {
  it('空の ArrayBuffer は例外を投げる', () => {
    const buf = new ArrayBuffer(0);
    assert.throws(() => parseBsp(buf), /TypeError|RangeError|非SPK/);
  });

  it('不正なマジックバイトは例外を投げる', () => {
    const buf = new ArrayBuffer(1024);
    const view = new DataView(buf);
    const bad = 'NOTSPK  ';
    for (let i = 0; i < 8; i++) view.setUint8(i, bad.charCodeAt(i));
    assert.throws(() => parseBsp(buf), /非SPK/);
  });
});

// =========================================================================
// 合成 BSP データ生成ヘルパー
// =========================================================================

/**
 * 最小限の DAF/SPK 合成バイナリを生成する
 *
 * Record 1 (0–1023):    ファイルヘッダー
 * Record 2 (1024–2047): サマリー（Sun/SSB, spkType=2 or 3 or 13）
 * Record 3 (2048–2167): データ
 *
 * J2000.0 での期待位置（Type 2/3 のみ）: x=100000, y=200000, z=300000 km
 */
function makeSyntheticBsp({ spkType = 2 } = {}) {
  const buf  = new ArrayBuffer(3 * 1024);
  const view = new DataView(buf);

  const writeStr = (s, offset, pad) => {
    for (let i = 0; i < pad; i++) {
      view.setUint8(offset + i, i < s.length ? s.charCodeAt(i) : 0x20);
    }
  };
  const le = true;

  // Record 1: ヘッダー
  writeStr('DAF/SPK ', 0, 8);
  view.setInt32(8,  2, le);   // ND
  view.setInt32(12, 6, le);   // NI
  writeStr('SyntheticBSP', 16, 60);
  view.setInt32(76, 2, le);   // FWARD
  view.setInt32(80, 2, le);   // BWARD
  writeStr('LTL-IEEE', 88, 8);

  // Record 2: サマリー
  // Type 2: firstAddr=257, lastAddr=271（rsize=11, n=1 → 15 doubles）
  // Type 3: firstAddr=257, lastAddr=268（rsize=8,  n=1 → 12 doubles）
  const isType3  = spkType === 3;
  const lastAddr = isType3 ? 268 : 271;
  const r2 = 1024;
  view.setFloat64(r2 + 0,  0.0, le);   // next rec
  view.setFloat64(r2 + 8,  0.0, le);
  view.setFloat64(r2 + 16, 1.0, le);   // NSUM
  const s = r2 + 24;
  view.setFloat64(s,       0.0, le);   // startJd（後で上書き）
  view.setFloat64(s + 8,   0.0, le);   // endJd（後で上書き）
  view.setInt32(s + 16,  10, le);  // target = Sun
  view.setInt32(s + 20,   0, le);  // center = SSB
  view.setInt32(s + 24,   1, le);  // frame
  view.setInt32(s + 28,   spkType, le);
  view.setInt32(s + 32,   257, le);
  view.setInt32(s + 36,   lastAddr, le);

  // Record 3: データ
  const r3 = 2048;
  view.setFloat64(r3 + 0,  0.0,     le);  // mid
  view.setFloat64(r3 + 8,  86400.0, le);  // radius

  if (isType3) {
    // [Xpos, Ypos, Zpos, Xvel, Yvel, Zvel] (ncoeff=1 each)
    view.setFloat64(r3 + 16, 100000.0, le);
    view.setFloat64(r3 + 24, 200000.0, le);
    view.setFloat64(r3 + 32, 300000.0, le);
    // Xvel, Yvel, Zvel = 0
    view.setFloat64(2112, -86400.0, le);  // init
    view.setFloat64(2120, 172800.0, le);  // intlen
    view.setFloat64(2128, 8.0,      le);  // rsize
    view.setFloat64(2136, 1.0,      le);  // n
  } else {
    // [Xpos×3, Ypos×3, Zpos×3] (ncoeff=3)
    view.setFloat64(r3 + 16, 100000.0, le);
    view.setFloat64(r3 + 40, 200000.0, le);
    view.setFloat64(r3 + 64, 300000.0, le);
    view.setFloat64(2136, -86400.0, le);  // init
    view.setFloat64(2144, 172800.0, le);  // intlen
    view.setFloat64(2152, 11.0,     le);  // rsize
    view.setFloat64(2160, 1.0,      le);  // n
  }

  // セグメントの JD 範囲（J2000.0 ± 1 day）
  const J2000 = 2451545.0;
  view.setFloat64(s,     J2000 - 1.0, le);
  view.setFloat64(s + 8, J2000 + 1.0, le);

  return buf;
}

// =========================================================================
// 合成 BSP を使った単体テスト — Type 2
// =========================================================================

describe('bsp-reader — 合成 BSP（Type 2）', () => {
  const J2000 = 2451545.0;

  it('セグメントが 1 件パースされる', () => {
    const bsp = parseBsp(makeSyntheticBsp());
    assert.strictEqual(bsp.segments.length, 1);
    assert.strictEqual(bsp.segments[0].target, 10);
    assert.strictEqual(bsp.segments[0].type, 2);
  });

  it('J2000.0 の位置が正しい（定数係数）', () => {
    const bsp = parseBsp(makeSyntheticBsp());
    const pos = bsp.getPosition(10, 0, J2000);
    assert.ok(Math.abs(pos[0] - 100000) < 1e-6, `x=${pos[0]}`);
    assert.ok(Math.abs(pos[1] - 200000) < 1e-6, `y=${pos[1]}`);
    assert.ok(Math.abs(pos[2] - 300000) < 1e-6, `z=${pos[2]}`);
  });

  it('computePosition と getPosition が一致する', () => {
    const bsp = parseBsp(makeSyntheticBsp());
    const p1 = bsp.getPosition(10, 0, J2000);
    const p2 = bsp.computePosition(10, 0, J2000);
    assert.ok(Math.abs(p1[0] - p2[0]) < 1e-9);
  });

  it('同一 target/center は [0,0,0]', () => {
    const bsp = parseBsp(makeSyntheticBsp());
    const pos = bsp.computePosition(10, 10, J2000);
    assert.deepStrictEqual(pos, [0, 0, 0]);
  });

  it('定数係数 → 速度 = 0', () => {
    const bsp = parseBsp(makeSyntheticBsp());
    const { position, velocity } = bsp.getPositionAndVelocity(10, 0, J2000);
    assert.ok(Math.abs(position[0] - 100000) < 1e-6);
    assert.ok(Math.abs(velocity[0]) < 1e-6);
    assert.ok(Math.abs(velocity[1]) < 1e-6);
    assert.ok(Math.abs(velocity[2]) < 1e-6);
  });

  it('範囲外 JD はエラーを投げる', () => {
    const bsp = parseBsp(makeSyntheticBsp());
    assert.throws(
      () => bsp.getPosition(10, 0, J2000 + 7),
      /セグメントが見つかりません|out of coverage/
    );
  });
});

// =========================================================================
// 合成 BSP を使った単体テスト — Type 3
// =========================================================================

describe('bsp-reader — 合成 BSP（Type 3）', () => {
  const J2000 = 2451545.0;

  it('Type 3 セグメントが spkType=3 としてパースされる', () => {
    const bsp = parseBsp(makeSyntheticBsp({ spkType: 3 }));
    assert.strictEqual(bsp.segments[0].type, 3);
  });

  it('Type 3 から位置が正しく取得できる', () => {
    const bsp = parseBsp(makeSyntheticBsp({ spkType: 3 }));
    const pos = bsp.getPosition(10, 0, J2000);
    assert.ok(Math.abs(pos[0] - 100000) < 1e-6, `x=${pos[0]}`);
    assert.ok(Math.abs(pos[1] - 200000) < 1e-6, `y=${pos[1]}`);
    assert.ok(Math.abs(pos[2] - 300000) < 1e-6, `z=${pos[2]}`);
  });

  it('Type 3 定数係数 → 速度 = 0', () => {
    const bsp = parseBsp(makeSyntheticBsp({ spkType: 3 }));
    const { velocity } = bsp.getPositionAndVelocity(10, 0, J2000);
    assert.ok(Math.abs(velocity[0]) < 1e-6);
    assert.ok(Math.abs(velocity[1]) < 1e-6);
    assert.ok(Math.abs(velocity[2]) < 1e-6);
  });

  it('Type 3 computePosition と getPosition が一致', () => {
    const bsp = parseBsp(makeSyntheticBsp({ spkType: 3 }));
    const p1 = bsp.getPosition(10, 0, J2000);
    const p2 = bsp.computePosition(10, 0, J2000);
    assert.ok(Math.abs(p1[0] - p2[0]) < 1e-9);
  });
});

// =========================================================================
// Type 13（スコープ外）
// =========================================================================

describe('bsp-reader — Type 13（非対応）', () => {
  const J2000 = 2451545.0;

  it('Type 13 セグメントは例外を投げる', () => {
    const bsp = parseBsp(makeSyntheticBsp({ spkType: 13 }));
    assert.throws(
      () => bsp.getPosition(10, 0, J2000),
      /未対応の SPK タイプ: 13/
    );
  });
});

// =========================================================================
// 2. 結合テスト（de440s.bsp が存在する場合のみ）
// =========================================================================

describe('bsp-reader — 実ファイル結合テスト (de440s.bsp)', { skip: !HAS_BSP }, () => {
  let bsp;

  before(async () => {
    const buf = await loadBsp(BSP_ABS);
    bsp = parseBsp(buf);
  });

  it('BspFile インスタンスが生成される', () => {
    assert.ok(bsp instanceof BspFile);
  });

  it('セグメントが 1 つ以上存在する', () => {
    assert.ok(bsp.segments.length > 0);
  });

  it('J2000.0 の Sun(10) 位置が SSB 基準で合理的な値', () => {
    const pos = bsp.getPosition(NAIF.SUN, NAIF.SSB, J2000_JD);
    assert.strictEqual(pos.length, 3);
    const dist = Math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2);
    const distAu = dist / AU_KM;
    assert.ok(distAu < 0.1, `太陽-SSB 距離 ≈ ${distAu.toFixed(4)} AU が 0.1 AU 以内`);
  });

  it('J2000.0 の Earth(399) 位置が合理的（SSB 経由合成）', () => {
    const pos = bsp.computePosition(NAIF.EARTH, NAIF.SSB, J2000_JD);
    const dist = Math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2);
    const distAu = dist / AU_KM;
    assert.ok(distAu > 0.9 && distAu < 1.1,
      `地球-SSB 距離 ≈ ${distAu.toFixed(4)} AU が 0.9〜1.1 AU の範囲`);
  });

  it('速度が取得できる: getPositionAndVelocity', () => {
    const { position, velocity } = bsp.getPositionAndVelocity(NAIF.SUN, NAIF.SSB, J2000_JD);
    assert.strictEqual(position.length, 3);
    assert.strictEqual(velocity.length, 3);
    const speed = Math.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2);
    assert.ok(speed > 0 && speed < 1e6);
  });
});

if (!HAS_BSP) {
  console.log(`⚠️  結合テストをスキップ: ${BSP_ABS} が存在しません`);
}
