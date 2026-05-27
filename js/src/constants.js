/**
 * constants.js — IAU / IERS / JPL 天文基本定数
 *
 * Layer 1: core（依存なし）
 * 出典:
 *   - IAU 2006 Resolution B1, B2, B3
 *   - Meeus "Astronomical Algorithms" 2nd ed.
 *   - Espenak & Meeus 2006 (NASA TechReport)
 *   - JPL DE440s / NAIF body codes
 *
 * ライセンス: MIT（数式・定数値は公知の科学的事実であり著作権なし）
 */

'use strict';

// =========================================================================
// 時刻・暦
// =========================================================================

/** J2000.0 エポック (2000-01-01T12:00:00 TT) のユリウス日 */
export const J2000_JD = 2451545.0;

/** 1 ユリウス世紀 = 36525 日 */
export const JULIAN_CENTURY = 36525.0;

/** 1 ユリウス年 = 365.25 日 */
export const JULIAN_YEAR = 365.25;

/** グレゴリオ暦施行日 (1582-10-15) の Julian Day Number */
export const GREGORIAN_CUTOVER_JDN = 2299161;

/** グレゴリオ暦施行年 */
export const GREGORIAN_CUTOVER_YEAR = 1582;

/** グレゴリオ暦施行月 */
export const GREGORIAN_CUTOVER_MONTH = 10;

/** グレゴリオ暦施行日 */
export const GREGORIAN_CUTOVER_DAY = 15;

// =========================================================================
// 時刻系オフセット
// =========================================================================

/**
 * TT - TAI の固定オフセット（秒）
 * IAU 1991 Resolution A.4
 */
export const TT_MINUS_TAI_SECONDS = 32.184;

/**
 * TDB - TT の最大振れ幅（秒）
 * Fairhead & Bretagnon 1990 の近似式の係数
 */
export const TDB_TT_AMPLITUDE = 0.001657;

/** TDB-TT の主要周期に関わる角速度係数（rad/日） */
export const TDB_TT_OMEGA = 628.3076;

// =========================================================================
// 物理定数
// =========================================================================

/** 天文単位 AU（m）— IAU 2012 */
export const AU_M = 149597870700;

/** 天文単位 AU（km） */
export const AU_KM = 149597870.700;

/** 光速（m/s）— IAU */
export const SPEED_OF_LIGHT_M_S = 299792458.0;

/** 光速（AU/日） */
export const SPEED_OF_LIGHT_AU_DAY = 173.14463267424034;

// =========================================================================
// 角度換算
// =========================================================================

/** π */
export const PI = Math.PI;

/** 2π */
export const TWO_PI = 2 * Math.PI;

/** 度 → ラジアン */
export const DEG_TO_RAD = Math.PI / 180.0;

/** ラジアン → 度 */
export const RAD_TO_DEG = 180.0 / Math.PI;

/** 秒角 → ラジアン */
export const ARCSEC_TO_RAD = Math.PI / (180.0 * 3600.0);

// =========================================================================
// JPL DE440s — NAIF ボディコード
// DE440s (.bsp) 内でセグメントを識別するターゲット番号
// 出典: JPL NAIF / SPICE Toolkit ユーザーガイド
// =========================================================================

export const NAIF = Object.freeze({
  /** 太陽系重心 (Solar System Barycenter) */
  SSB: 0,

  /** 水星重心 */
  MERCURY_BARYCENTER: 1,
  /** 金星重心 */
  VENUS_BARYCENTER: 2,
  /** 地球月系重心 */
  EMB: 3,
  /** 火星重心 */
  MARS_BARYCENTER: 4,
  /** 木星重心 */
  JUPITER_BARYCENTER: 5,
  /** 土星重心 */
  SATURN_BARYCENTER: 6,
  /** 天王星重心 */
  URANUS_BARYCENTER: 7,
  /** 海王星重心 */
  NEPTUNE_BARYCENTER: 8,
  /** 冥王星重心 */
  PLUTO_BARYCENTER: 9,

  /** 太陽 */
  SUN: 10,

  /** 月 */
  MOON: 301,
  /** 地球 */
  EARTH: 399,

  /** 水星（inertial center） */
  MERCURY: 199,
  /** 金星 */
  VENUS: 299,
  /** 火星 */
  MARS: 499,
});

/**
 * 占星術で使う天体の NAIF コードマッピング
 */
export const PLANET_NAIF = Object.freeze({
  Sun:     NAIF.SUN,
  Moon:    NAIF.MOON,
  Mercury: NAIF.MERCURY_BARYCENTER,
  Venus:   NAIF.VENUS_BARYCENTER,
  Earth:   NAIF.EARTH,
  Mars:    NAIF.MARS_BARYCENTER,
  Jupiter: NAIF.JUPITER_BARYCENTER,
  Saturn:  NAIF.SATURN_BARYCENTER,
  Uranus:  NAIF.URANUS_BARYCENTER,
  Neptune: NAIF.NEPTUNE_BARYCENTER,
  Pluto:   NAIF.PLUTO_BARYCENTER,
  EMB:     NAIF.EMB,
  SSB:     NAIF.SSB,
});
