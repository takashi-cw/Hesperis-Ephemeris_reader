// constants.rs — 天文基本定数・NAIF ボディコード
//
// 出典:
//   - IAU 2006 Resolution B1, B2, B3
//   - JPL DE440s / NAIF SPICE Toolkit ユーザーガイド
//
// ライセンス: MIT（定数値は公知の科学的事実であり著作権なし）

// MARK: - 時刻・暦

/// J2000.0 エポック（2000-01-01T12:00:00 TT）のユリウス日
pub const J2000_JD: f64 = 2_451_545.0;

/// 1 日 = 86400 秒
pub const SECS_PER_DAY: f64 = 86_400.0;

/// 1 ユリウス世紀 = 36525 日
pub const JULIAN_CENTURY: f64 = 36_525.0;

// MARK: - 物理定数

/// 天文単位（km）— IAU 2012
pub const AU_KM: f64 = 149_597_870.700;

/// 光速（m/s）— IAU
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;

/// 光速（AU/日）
pub const SPEED_OF_LIGHT_AU_DAY: f64 = 173.144_632_674_240_34;

// MARK: - NAIF ボディコード
//
// DE440s (.bsp) 内でセグメントを識別するターゲット番号
// 出典: JPL NAIF / SPICE Toolkit ユーザーガイド

pub mod naif {
    /// 太陽系重心 (Solar System Barycenter)
    pub const SSB: i32 = 0;
    /// 水星重心
    pub const MERCURY_BARYCENTER: i32 = 1;
    /// 金星重心
    pub const VENUS_BARYCENTER: i32 = 2;
    /// 地球月系重心 (Earth-Moon Barycenter)
    pub const EMB: i32 = 3;
    /// 火星重心
    pub const MARS_BARYCENTER: i32 = 4;
    /// 木星重心
    pub const JUPITER_BARYCENTER: i32 = 5;
    /// 土星重心
    pub const SATURN_BARYCENTER: i32 = 6;
    /// 天王星重心
    pub const URANUS_BARYCENTER: i32 = 7;
    /// 海王星重心
    pub const NEPTUNE_BARYCENTER: i32 = 8;
    /// 冥王星重心
    pub const PLUTO_BARYCENTER: i32 = 9;
    /// 太陽
    pub const SUN: i32 = 10;
    /// 月
    pub const MOON: i32 = 301;
    /// 地球
    pub const EARTH: i32 = 399;
    /// 水星（天体中心）
    pub const MERCURY: i32 = 199;
    /// 金星（天体中心）
    pub const VENUS: i32 = 299;
    /// 火星（天体中心）
    pub const MARS: i32 = 499;
}
