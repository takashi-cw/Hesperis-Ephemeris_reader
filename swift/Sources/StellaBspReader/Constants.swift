// Constants.swift — IAU / IERS / JPL 天文基本定数
//
// Layer 1: core（依存なし）
// 出典:
//   - IAU 2006 Resolution B1, B2, B3
//   - Meeus "Astronomical Algorithms" 2nd ed.
//   - JPL DE440s / NAIF body codes
//
// ライセンス: MIT（数式・定数値は公知の科学的事実であり著作権なし）

// MARK: - 時刻・暦

/// J2000.0 エポック (2000-01-01T12:00:00 TT) のユリウス日
public let J2000_JD: Double = 2_451_545.0

/// 1 ユリウス世紀 = 36525 日
public let JULIAN_CENTURY: Double = 36_525.0

/// 1 ユリウス年 = 365.25 日
public let JULIAN_YEAR: Double = 365.25

/// グレゴリオ暦施行日 (1582-10-15) の Julian Day Number
public let GREGORIAN_CUTOVER_JDN: Int = 2_299_161

// MARK: - 時刻系オフセット

/// TT - TAI の固定オフセット（秒）— IAU 1991 Resolution A.4
public let TT_MINUS_TAI_SECONDS: Double = 32.184

/// TDB - TT の最大振れ幅（秒）— Fairhead & Bretagnon 1990
public let TDB_TT_AMPLITUDE: Double = 0.001657

/// TDB-TT の主要周期に関わる角速度係数（rad/Julian世紀）
/// Fairhead & Bretagnon (1990) の近似式:
///   TDB − TT ≈ 0.001657 × sin(628.3076 × T + 6.2401)  [秒]
/// ここで T は J2000.0 からの Julian 世紀数。
/// 628.3076 = 100 × 2π ≈ 地球の平均運動 × 100 [rad/世紀]
/// ※ T が「日」ではなく「世紀」である点に注意。日換算と混同すると
///    周期が 365.25日 ではなく 約15分 になる誤計算が生じる。
public let TDB_TT_OMEGA: Double = 628.3076

// MARK: - 物理定数

/// 天文単位 AU（m）— IAU 2012
public let AU_M: Double = 149_597_870_700.0

/// 天文単位 AU（km）
public let AU_KM: Double = 149_597_870.700

/// 光速（m/s）— IAU
public let SPEED_OF_LIGHT_M_S: Double = 299_792_458.0

/// 光速（AU/日）
public let SPEED_OF_LIGHT_AU_DAY: Double = 173.14463267424034

// MARK: - 角度換算

/// 度 → ラジアン
public let DEG_TO_RAD: Double = Double.pi / 180.0

/// ラジアン → 度
public let RAD_TO_DEG: Double = 180.0 / Double.pi

/// 秒角 → ラジアン
public let ARCSEC_TO_RAD: Double = Double.pi / (180.0 * 3600.0)

// MARK: - JPL DE440s / NAIF ボディコード
//
// DE440s (.bsp) 内でセグメントを識別するターゲット番号
// 出典: JPL NAIF / SPICE Toolkit ユーザーガイド

public enum NAIF {
    /// 太陽系重心 (Solar System Barycenter)
    public static let SSB: Int32 = 0
    /// 水星重心
    public static let MERCURY_BARYCENTER: Int32 = 1
    /// 金星重心
    public static let VENUS_BARYCENTER: Int32 = 2
    /// 地球月系重心
    public static let EMB: Int32 = 3
    /// 火星重心
    public static let MARS_BARYCENTER: Int32 = 4
    /// 木星重心
    public static let JUPITER_BARYCENTER: Int32 = 5
    /// 土星重心
    public static let SATURN_BARYCENTER: Int32 = 6
    /// 天王星重心
    public static let URANUS_BARYCENTER: Int32 = 7
    /// 海王星重心
    public static let NEPTUNE_BARYCENTER: Int32 = 8
    /// 冥王星重心
    public static let PLUTO_BARYCENTER: Int32 = 9
    /// 太陽
    public static let SUN: Int32 = 10
    /// 月
    public static let MOON: Int32 = 301
    /// 地球
    public static let EARTH: Int32 = 399
    /// 水星（重心ではなく天体中心）
    public static let MERCURY: Int32 = 199
    /// 金星
    public static let VENUS: Int32 = 299
    /// 火星
    public static let MARS: Int32 = 499
}

// MARK: - 太陽系不変面（SIRF）固定定数
//
// 出典: Souami & Souchay 2012, IAU WGCCRE 2015
// 参照: spacefield/ephem/solar_invariable.py _R_SIRF_J2000

/// J2000 黄道傾斜角（度）
public let EPS_J2000_DEG: Double = 23.4392911

/// 太陽系不変面の極ベクトル（黄道 J2000 基準）黄経・黄緯
/// Souami & Souchay 2012
public let INV_ECLON_DEG: Double = 107.5892
public let INV_ECLAT_DEG: Double =  88.4220

/// 太陽自転軸（ICRF 赤道座標）赤経・赤緯
/// IAU WGCCRE 2015
public let SUN_RA_DEG:  Double = 286.13
public let SUN_DEC_DEG: Double =  63.87

/// ICRS → 太陽系不変面（SIRF J2000）回転行列（行優先 3×3）
/// spacefield/ephem/solar_invariable.py _R_SIRF_J2000 と同値
/// 行0 = X軸（0°方向 L0）、行1 = Y軸、行2 = Z軸（不変面法線）
public let R_SIRF_J2000: [[Double]] = [
    [ 0.4004119179615650,  0.8487495560218467,  0.3453903402051928],
    [-0.9162974658726827,  0.3742866696641715,  0.1425076943431592],
    [-0.0083216578661234, -0.3735420726831468,  0.9275759106110604],
]

// MARK: - IAU 銀河座標回転行列
//
// 出典: IAU 1958 銀河座標系 / Liu et al. (2011) Hipparcos 改訂値, J2000 固定
// 参照: spacefield/src/spacefield/solar_invariable.py _GAL_FROM_ICRF
// ICRS 赤道座標 → 銀河座標（l, b）への回転行列（行優先 3×3）
// 行0 = 銀河 X 軸（銀河中心方向 l=0°）、行1 = Y 軸（l=90°）、行2 = Z 軸（北銀極）

public let R_GALACTIC: [[Double]] = [
    [-0.0548755604, -0.8734370902, -0.4838350155],
    [ 0.4941094279, -0.4448296300,  0.7469822445],
    [-0.8676661490, -0.1980763734,  0.4559837762],
]

// MARK: - 惑星軌道要素（J2000.0 エポック）
//
// 出典: NASA/JPL HORIZONS System（パブリックドメイン）
// 参照: spacefield/ephem/orbital_elements.py PLANET_ORBITAL_ELEMENTS

/// 惑星軌道要素（軌道傾斜角・昇交点黄経）
public struct PlanetOrbitalElements: Sendable {
    /// 軌道傾斜角（度、黄道に対する傾き）
    public let inclination: Double
    /// 昇交点黄経 Ω（度）
    public let ascendingNode: Double
}

/// NAIF コード → 軌道要素（J2000.0）
public let PLANET_ORBITAL_ELEMENTS: [Int32: PlanetOrbitalElements] = [
    NAIF.MERCURY_BARYCENTER: PlanetOrbitalElements(inclination:  7.0050, ascendingNode:  48.3313),
    NAIF.VENUS_BARYCENTER:   PlanetOrbitalElements(inclination:  3.3947, ascendingNode:  76.6799),
    NAIF.EARTH:              PlanetOrbitalElements(inclination:  0.0000, ascendingNode:   0.0000),
    NAIF.MARS_BARYCENTER:    PlanetOrbitalElements(inclination:  1.8506, ascendingNode:  49.5785),
    NAIF.JUPITER_BARYCENTER: PlanetOrbitalElements(inclination:  1.3053, ascendingNode: 100.4644),
    NAIF.SATURN_BARYCENTER:  PlanetOrbitalElements(inclination:  2.4845, ascendingNode: 113.6655),
    NAIF.URANUS_BARYCENTER:  PlanetOrbitalElements(inclination:  0.7733, ascendingNode:  74.0060),
    NAIF.NEPTUNE_BARYCENTER: PlanetOrbitalElements(inclination:  1.7700, ascendingNode: 131.7841),
    NAIF.PLUTO_BARYCENTER:   PlanetOrbitalElements(inclination: 17.1417, ascendingNode: 110.2990),
    // ⚠ MOON の ascendingNode は概算値（0.0 固定）。
    // 月の昇交点は 18.6 年周期で 0°〜360° を変動するため静的テーブルでは正確に表現できない。
    // 正確な of-date 昇交点を求めるには「地心 r × v」が必要だが、
    // eclipticToOrbitalPlaneOfDate() は日心 r × v を使うため月には非対応。
    // 現状は傾斜角のみ参考値として保持し、昇交点は未実装扱いとする。
    NAIF.MOON:               PlanetOrbitalElements(inclination:  5.1454, ascendingNode:   0.0000),
]

// MARK: - 惑星自転軸方向（IAU WGCCRE 2015/2018）
//
// 出典: IAU Working Group on Cartographic Coordinates and Rotational Elements 2015/2018
// 参照: spacefield/ephem/planet_vernal_equinox.py PLANET_ROTATION_AXIS
// ICRF J2000 赤道座標（赤経 RA / 赤緯 Dec）

/// 惑星自転軸方向（J2000 ICRF 赤道座標）
public struct PlanetRotationAxis: Sendable {
    /// 赤経（度）
    public let raJ2000: Double
    /// 赤緯（度）
    public let decJ2000: Double
}

/// NAIF コード → 自転軸方向
public let PLANET_ROTATION_AXIS: [Int32: PlanetRotationAxis] = [
    NAIF.SUN:                PlanetRotationAxis(raJ2000: 286.13,     decJ2000:  63.87),   // IAU WGCCRE 2015
    NAIF.MERCURY_BARYCENTER: PlanetRotationAxis(raJ2000: 281.01,     decJ2000:  61.45),
    NAIF.VENUS_BARYCENTER:   PlanetRotationAxis(raJ2000: 272.76,     decJ2000: -67.16),   // 逆行自転
    NAIF.MARS_BARYCENTER:    PlanetRotationAxis(raJ2000: 317.68143,  decJ2000:  52.88650),// IAU 2015
    NAIF.JUPITER_BARYCENTER: PlanetRotationAxis(raJ2000: 268.056595, decJ2000:  64.495303),// IAU 2018
    NAIF.SATURN_BARYCENTER:  PlanetRotationAxis(raJ2000:  40.589,    decJ2000:  83.538),
    NAIF.URANUS_BARYCENTER:  PlanetRotationAxis(raJ2000: 257.311,    decJ2000: -15.175),  // 逆行
    NAIF.NEPTUNE_BARYCENTER: PlanetRotationAxis(raJ2000: 299.36,     decJ2000:  43.46),
    NAIF.MOON:               PlanetRotationAxis(raJ2000: 270.00,     decJ2000:  66.54),
]

/// 惑星名 → NAIF コード のマッピング
public let PLANET_NAIF: [String: Int32] = [
    "Sun":     NAIF.SUN,
    "Moon":    NAIF.MOON,
    "Mercury": NAIF.MERCURY_BARYCENTER,
    "Venus":   NAIF.VENUS_BARYCENTER,
    "Earth":   NAIF.EARTH,
    "Mars":    NAIF.MARS_BARYCENTER,
    "Jupiter": NAIF.JUPITER_BARYCENTER,
    "Saturn":  NAIF.SATURN_BARYCENTER,
    "Uranus":  NAIF.URANUS_BARYCENTER,
    "Neptune": NAIF.NEPTUNE_BARYCENTER,
    "Pluto":   NAIF.PLUTO_BARYCENTER,
    "EMB":     NAIF.EMB,
    "SSB":     NAIF.SSB,
]
