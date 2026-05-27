# constants.jl — 天文基本定数・NAIF ボディコード
#
# 出典:
#   - IAU 2006 Resolution B1, B2, B3
#   - JPL DE440s / NAIF SPICE Toolkit ユーザーガイド
#
# ライセンス: MIT（定数値は公知の科学的事実であり著作権なし）

# --- 時刻・暦 ---

"J2000.0 エポック（2000-01-01T12:00:00 TT）のユリウス日"
const J2000_JD = 2_451_545.0

"1 日 = 86400 秒"
const SECS_PER_DAY = 86_400.0

"1 ユリウス世紀 = 36525 日"
const JULIAN_CENTURY = 36_525.0

# --- 物理定数 ---

"天文単位（km）— IAU 2012"
const AU_KM = 149_597_870.700

"光速（m/s）— IAU"
const SPEED_OF_LIGHT_M_S = 299_792_458.0

"光速（AU/日）"
const SPEED_OF_LIGHT_AU_DAY = 173.144_632_674_240_34

# --- NAIF ボディコード ---
#
# DE440s (.bsp) 内でセグメントを識別するターゲット番号
# 出典: JPL NAIF / SPICE Toolkit ユーザーガイド

module Naif
    "太陽系重心 (Solar System Barycenter)"
    const SSB                = Int32(0)
    "水星重心"
    const MERCURY_BARYCENTER = Int32(1)
    "金星重心"
    const VENUS_BARYCENTER   = Int32(2)
    "地球月系重心 (Earth-Moon Barycenter)"
    const EMB                = Int32(3)
    "火星重心"
    const MARS_BARYCENTER    = Int32(4)
    "木星重心"
    const JUPITER_BARYCENTER = Int32(5)
    "土星重心"
    const SATURN_BARYCENTER  = Int32(6)
    "天王星重心"
    const URANUS_BARYCENTER  = Int32(7)
    "海王星重心"
    const NEPTUNE_BARYCENTER = Int32(8)
    "冥王星重心"
    const PLUTO_BARYCENTER   = Int32(9)
    "太陽"
    const SUN                = Int32(10)
    "月"
    const MOON               = Int32(301)
    "地球"
    const EARTH              = Int32(399)
    "水星（天体中心）"
    const MERCURY            = Int32(199)
    "金星（天体中心）"
    const VENUS              = Int32(299)
    "火星（天体中心）"
    const MARS               = Int32(499)
end
