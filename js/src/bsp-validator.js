/**
 * bsp-validator.js — BSP 天体暦のカバー範囲検証
 *
 * BSP ファイルのセグメント時刻範囲（J2000.0 からの秒数）を JD に変換し、
 * 入力日時が天体暦のカバー範囲内かを検証する。
 *
 * 使い方:
 *   import { assertInCoverage, getCoverageJd } from './bsp-validator.js';
 *   assertInCoverage(jdTdb, bspFile);  // 範囲外なら RangeError を throw
 *
 * ライセンス: MIT
 */

import { J2000_JD } from './constants.js';

const S_PER_DAY = 86400.0;
const NAIF_SUN  = 10;

/**
 * BSP ファイルのカバー範囲を JD で返す
 *
 * @param {import('./bsp-reader.js').BspFile} bspFile
 * @param {number} [naifTarget=10]  代表天体の NAIF コード（デフォルト: 太陽）
 * @returns {{ startJd: number, endJd: number }}
 */
export function getCoverageJd(bspFile, naifTarget = NAIF_SUN) {
  const segs = bspFile.segments.filter(s => s.target === naifTarget);

  if (segs.length === 0) {
    const all    = bspFile.segments;
    const minSec = Math.min(...all.map(s => s.startJd));
    const maxSec = Math.max(...all.map(s => s.endJd));
    return {
      startJd: J2000_JD + minSec / S_PER_DAY,
      endJd:   J2000_JD + maxSec / S_PER_DAY,
    };
  }

  const minSec = Math.min(...segs.map(s => s.startJd));
  const maxSec = Math.max(...segs.map(s => s.endJd));
  return {
    startJd: J2000_JD + minSec / S_PER_DAY,
    endJd:   J2000_JD + maxSec / S_PER_DAY,
  };
}

/**
 * カバー範囲を人間が読める文字列で返す
 *
 * @param {{ startJd: number, endJd: number }} coverage
 * @returns {string}  例: "AD1850〜AD2150"（de440s.bsp の場合）
 */
export function formatCoverageMessage(coverage) {
  const startYear = _jdToYear(coverage.startJd);
  const endYear   = _jdToYear(coverage.endJd);
  const start = startYear < 0
    ? `BC${Math.abs(Math.ceil(startYear))}`
    : `AD${Math.floor(startYear)}`;
  const end = `AD${Math.floor(endYear)}`;
  return `${start}〜${end}`;
}

/**
 * JD が BSP のカバー範囲内かを検証する
 *
 * @param {number} jdTdb   検証する JD（TDB）
 * @param {import('./bsp-reader.js').BspFile} bspFile
 * @throws {RangeError} 範囲外の場合
 */
export function assertInCoverage(jdTdb, bspFile) {
  const coverage = getCoverageJd(bspFile);
  if (jdTdb < coverage.startJd || jdTdb > coverage.endJd) {
    const range      = formatCoverageMessage(coverage);
    const inputYear  = _jdToYear(jdTdb);
    const inputLabel = inputYear < 0
      ? `BC${Math.abs(Math.ceil(inputYear))}`
      : `AD${Math.floor(inputYear)}`;
    throw new RangeError(
      `天体暦の範囲外です（${inputLabel}）。カバー範囲: ${range}。`
    );
  }
}

function _jdToYear(jd) {
  return 2000.0 + (jd - J2000_JD) / 365.25;
}
