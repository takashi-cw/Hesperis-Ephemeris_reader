//! coord.rs — 座標変換純粋関数モジュール
//!
//! coordinates.py (MIT) + apparent.py (MIT) のホットパスを Rust 移植。
//!
//! 移植元:
//!   spacefield/src/spacefield/ephem/coordinates.py（MIT）
//!   spacefield/src/spacefield/ephem/apparent.py（MIT）
//!   Stella-JS/public/src/astro/coordinates.js（MIT）
//!
//! 定数テーブル出典:
//!   IERS Conventions 2010（IAU 2000B 章動 77 項）
//!   Capitaine et al. 2003, A&A 412（IAU 2006 3角度歳差）
//!   SOFA/ERFA（定数 / パブリックドメイン相当）
//!   → MIT 配布に問題なし
//!
//! 設計方針:
//!   - 全関数が純粋関数（引数のみ使用・副作用なし）
//!   - BSP 読み取りとの分離を維持（coord.rs はデータを触らない）
//!   - lib.rs から `use crate::coord::*` で利用
//!
//! ライセンス: MIT

use std::f64::consts::PI;

// ============================================================
// 基本定数
// ============================================================

pub const J2000_JD: f64 = 2451545.0;
pub const JULIAN_CENTURY: f64 = 36525.0;

const D2R: f64 = PI / 180.0;
const A2R: f64 = PI / (180.0 * 3600.0);
const R2D: f64 = 180.0 / PI;
const TWO_PI: f64 = 2.0 * PI;

// IAU 2006 平均黄道傾斜角 J2000.0 [arcsec]
const EPS0_J2000_ARCSEC: f64 = 84381.406;

// ICRS フレームバイアス補正 ψ_bar(T=0) [arcsec]
// Capitaine Eq.37 の ψ_A(T=0)=0 と FW4 モデル間の差。これを psiA に加算して
// JPL Horizons との黄経差を 0.042" → 0.001" 台に改善する。
const PSI_BIAS_ARCSEC: f64 = -0.041775;

// 光速 [km/day]
pub const C_KM_PER_DAY: f64 = 299792.458 * 86400.0;

// 太陽シュヴァルツシルト半径 2GM_sun/c² [km]
const DEFL_CONST_KM: f64 = 2.953250;

// 年周光行差補正の有限差分幅: ±0.5 s
pub const ABERR_DT_DAYS: f64 = 0.5 / 86400.0;

// 速度計算の有限差分幅: ±30 s
pub const SPEED_DT_DAYS: f64 = 30.0 / 86400.0;

// NAIF コード（lib.rs と共有）
pub const NAIF_SSB: i32 = 0;
pub const NAIF_SUN: i32 = 10;
pub const NAIF_EARTH: i32 = 399;

// ============================================================
// IAU 2006 Capitaine 3角度歳差係数
// ============================================================
// 出典: Capitaine et al. 2003, A&A 412 Eq.37 / ERFA eraP06e
//
// 使い方:
//   T = (jd - 2451545.0) / 36525.0
//   psiA [arcsec] = PSI_A[0]*T + PSI_A[1]*T² + ... + PSI_BIAS_ARCSEC
//   omgA [arcsec] = EPS0_J2000_ARCSEC + OMG_A_POLY[0]*T + ...
//   chiA [arcsec] = CHI_A[0]*T + CHI_A[1]*T² + ...

// ψ_A (黄経一般歳差) 係数 [arcsec / T^n]
const PREC_PSI_A: [f64; 5] = [
    5038.481507, -1.0790069, -0.00114045, 0.000132851, -0.0000000951,
];

// ω_A (傾斜角一般歳差) 多項式係数 [arcsec / T^n]（基準値は EPS0_J2000_ARCSEC）
const PREC_OMG_A_POLY: [f64; 5] = [
    -0.025754, 0.0512623, -0.00772503, -0.000000467, 0.0000003337,
];

// χ_A (赤道歳差) 係数 [arcsec / T^n]
const PREC_CHI_A: [f64; 5] = [
    10.556403, -2.3814292, -0.00121197, 0.000170663, -0.0000000560,
];

// ============================================================
// IAU 2000B 章動 5 基本引数係数
// ============================================================
// 出典: IERS Conventions 2010 / Capitaine et al. 2003
// 順序: l (月平均近点角), lp (太陽平均近点角), F (月緯度引数), D (月平均離角), Ω (月昇交点黄経)
// 各行: [定数項, T, T², T³, T⁴] [arcsec]

const NUT_FUND_COEFFS: [[f64; 5]; 5] = [
    [485868.249036,  1717915923.2178,  31.8792,   0.051635, -0.00024470], // l
    [1287104.79305,   129596581.0481,  -0.5532,   0.000136, -0.00001149], // lp
    [ 335779.526232, 1739527262.8478, -12.7512,  -0.001037,  0.00000417], // F
    [1072260.70369,  1602961601.2090,  -6.3706,   0.006593, -0.00003169], // D
    [ 450160.398036,   -6962890.5431,   7.4722,   0.007702, -0.00005939], // Om
];

// ============================================================
// IAU 2000B 章動テーブル（77 項）
// ============================================================
// 出典: IERS Conventions 2010, Appendix B / Mathews, Herring, Buffett (2002)
//
// 列: (n_l, n_lp, n_F, n_D, n_Om, AA, BB, CC, DD, EE, FF)
//   n_*: 基本引数の整数係数
//   AA: dpsi sin 主係数 [×1e-7 arcsec]
//   BB: dpsi sin T 係数 [×1e-7 arcsec]
//   CC: dpsi cos 主係数 [×1e-7 arcsec]
//   DD: deps cos 主係数 [×1e-7 arcsec]
//   EE: deps cos T 係数 [×1e-7 arcsec]
//   FF: deps sin 主係数 [×1e-7 arcsec]
//
// 計算式:
//   arg = n_l*l + n_lp*lp + n_F*F + n_D*D + n_Om*Om
//   dpsi += (AA + BB*T)*sin(arg) + CC*cos(arg)
//   deps += (DD + EE*T)*cos(arg) + FF*sin(arg)
//   最終: dpsi[arcsec] = Σ * 1e-7

const NUT77: [(i32, i32, i32, i32, i32, f64, f64, f64, f64, f64, f64); 77] = [
    //  nl, nlp,  nF,  nD, nOm,          AA,      BB,      CC,          DD,      EE,      FF
    (   0,   0,   0,   0,   1, -172064161.0, -174666.0,  33386.0,  92052331.0,   9086.0,  15377.0),
    (   0,   0,   2,  -2,   2,  -13170906.0,   -1675.0, -13696.0,   5730336.0,  -3015.0,  -4587.0),
    (   0,   0,   2,   0,   2,   -2276413.0,    -234.0,   2796.0,    978459.0,   -485.0,   1374.0),
    (   0,   0,   0,   0,   2,    2074554.0,     207.0,   -698.0,   -897492.0,    470.0,   -291.0),
    (   0,   1,   0,   0,   0,    1475877.0,   -3633.0,  11817.0,     73871.0,   -184.0,  -1924.0),
    (   0,   1,   2,  -2,   2,    -516821.0,    1226.0,   -524.0,    224386.0,   -677.0,   -174.0),
    (   1,   0,   0,   0,   0,     711159.0,      73.0,   -872.0,     -6750.0,      0.0,    358.0),
    (   0,   0,   2,   0,   1,    -387298.0,    -367.0,    380.0,    200728.0,     18.0,    318.0),
    (   1,   0,   2,   0,   2,    -301461.0,     -36.0,    816.0,    129025.0,    -63.0,    367.0),
    (   0,  -1,   2,  -2,   2,     215829.0,    -494.0,    111.0,    -95929.0,    299.0,    132.0),
    (   0,   0,   2,  -2,   1,     128227.0,     137.0,    181.0,    -68982.0,     -9.0,     39.0),
    (  -1,   0,   2,   0,   2,     123457.0,      11.0,     19.0,    -53311.0,     32.0,     -4.0),
    (  -1,   0,   0,   2,   0,     156994.0,      10.0,   -168.0,     -1235.0,      0.0,     82.0),
    (   1,   0,   0,   0,   1,      63110.0,      63.0,     27.0,    -33228.0,      0.0,     -9.0),
    (  -1,   0,   0,   0,   1,     -57976.0,     -63.0,   -189.0,     31429.0,      0.0,    -75.0),
    (  -1,   0,   2,   2,   2,     -59641.0,     -11.0,    149.0,     25543.0,    -11.0,     66.0),
    (   1,   0,   2,   0,   1,     -51613.0,     -42.0,    129.0,     26366.0,      0.0,     78.0),
    (  -2,   0,   2,   0,   1,      45893.0,      50.0,     31.0,    -24236.0,    -10.0,     20.0),
    (   0,   0,   0,   2,   0,      63384.0,      11.0,   -150.0,     -1220.0,      0.0,     29.0),
    (   0,   0,   2,   2,   2,     -38571.0,      -1.0,    158.0,     16452.0,    -11.0,     68.0),
    (   0,  -2,   2,  -2,   2,      32481.0,       0.0,      0.0,    -13870.0,      0.0,      0.0),
    (  -2,   0,   0,   2,   0,     -47722.0,       0.0,    -18.0,       477.0,      0.0,    -25.0),
    (   2,   0,   2,   0,   2,     -31046.0,      -1.0,    131.0,     13238.0,    -11.0,     59.0),
    (   1,   0,   2,  -2,   2,      28593.0,       0.0,     -1.0,    -12338.0,     10.0,     -3.0),
    (  -1,   0,   2,   0,   1,      20441.0,      21.0,     10.0,    -10758.0,      0.0,     -3.0),
    (   2,   0,   0,   0,   0,      29243.0,       0.0,    -74.0,      -609.0,      0.0,     13.0),
    (   0,   0,   2,   0,   0,      25887.0,       0.0,    -66.0,      -550.0,      0.0,     11.0),
    (   0,   1,   0,   0,   1,     -14053.0,     -25.0,     79.0,      8551.0,     -2.0,    -45.0),
    (  -1,   0,   0,   2,   1,      15164.0,      10.0,     11.0,     -8001.0,      0.0,     -1.0),
    (   0,   2,   2,  -2,   2,     -15794.0,      72.0,    -16.0,      6850.0,    -42.0,     -5.0),
    (   0,   0,  -2,   2,   0,      21783.0,       0.0,     13.0,      -167.0,      0.0,     13.0),
    (   1,   0,   0,  -2,   1,     -12873.0,     -10.0,    -37.0,      6953.0,      0.0,    -14.0),
    (   0,  -1,   0,   0,   1,     -12654.0,      11.0,     63.0,      6415.0,      0.0,     26.0),
    (  -1,   0,   2,   2,   1,     -10204.0,       0.0,     25.0,      5222.0,      0.0,     15.0),
    (   0,   2,   0,   0,   0,      16707.0,     -85.0,    -10.0,       168.0,     -1.0,     10.0),
    (   1,   0,   2,   2,   2,      -7691.0,       0.0,     44.0,      3268.0,      0.0,     19.0),
    (  -2,   0,   2,   0,   0,     -11024.0,       0.0,    -14.0,       104.0,      0.0,      2.0),
    (   0,   1,   2,   0,   2,       7566.0,     -21.0,    -11.0,     -3250.0,      0.0,     -5.0),
    (   0,   0,   2,   2,   1,      -6637.0,     -11.0,     25.0,      3353.0,      0.0,     14.0),
    (   0,  -1,   2,   0,   2,      -7141.0,      21.0,      8.0,      3070.0,      0.0,      4.0),
    (   0,   0,   0,   2,   1,      -6302.0,     -11.0,      2.0,      3272.0,      0.0,      4.0),
    (   1,   0,   2,  -2,   1,       5800.0,      10.0,      2.0,     -3045.0,      0.0,     -1.0),
    (   2,   0,   2,  -2,   2,       6443.0,       0.0,     -7.0,     -2768.0,      0.0,     -4.0),
    (  -2,   0,   0,   2,   1,      -5774.0,     -11.0,    -15.0,      3041.0,      0.0,     -5.0),
    (   2,   0,   2,   0,   1,      -5350.0,       0.0,     21.0,      2695.0,      0.0,     12.0),
    (   0,  -1,   2,  -2,   1,      -4752.0,     -11.0,     -3.0,      2719.0,      0.0,     -3.0),
    (   0,   0,   0,  -2,   1,      -4940.0,     -11.0,    -21.0,      2720.0,      0.0,     -9.0),
    (  -1,  -1,   0,   2,   0,       7350.0,       0.0,     -8.0,       -51.0,      0.0,      4.0),
    (   2,   0,   0,  -2,   1,       4065.0,       0.0,      6.0,     -2206.0,      0.0,      1.0),
    (   1,   0,   0,   2,   0,       6579.0,       0.0,    -24.0,      -199.0,      0.0,      2.0),
    (   0,   1,   2,  -2,   1,       3579.0,       0.0,      5.0,     -1900.0,      0.0,      1.0),
    (   1,  -1,   0,   0,   0,       4725.0,       0.0,     -6.0,       -41.0,      0.0,      3.0),
    (  -2,   0,   2,   0,   2,      -3075.0,       0.0,     -2.0,      1313.0,      0.0,     -1.0),
    (   3,   0,   2,   0,   2,      -2904.0,       0.0,     15.0,      1233.0,      0.0,      7.0),
    (   0,  -1,   0,   2,   0,       4348.0,       0.0,    -10.0,       -81.0,      0.0,      2.0),
    (   1,  -1,   2,   0,   2,      -2878.0,       0.0,      8.0,      1232.0,      0.0,      4.0),
    (   0,   0,   0,   1,   0,      -4230.0,       0.0,      5.0,       -20.0,      0.0,     -2.0),
    (  -1,  -1,   2,   2,   2,      -2819.0,       0.0,      7.0,      1207.0,      0.0,      3.0),
    (  -1,   0,   2,   0,   0,      -4056.0,       0.0,      5.0,        40.0,      0.0,     -2.0),
    (   0,  -1,   2,   2,   2,      -2647.0,       0.0,     11.0,      1129.0,      0.0,      5.0),
    (  -2,   0,   0,   0,   1,      -2294.0,       0.0,    -10.0,      1266.0,      0.0,     -4.0),
    (   1,   1,   2,   0,   2,       2481.0,       0.0,     -7.0,     -1062.0,      0.0,     -3.0),
    (   2,   0,   0,   0,   1,       2179.0,       0.0,     -2.0,     -1129.0,      0.0,     -2.0),
    (  -1,   1,   0,   1,   0,       3276.0,       0.0,      1.0,        -9.0,      0.0,      0.0),
    (   1,   1,   0,   0,   0,      -3389.0,       0.0,      5.0,        35.0,      0.0,     -2.0),
    (   1,   0,   2,   0,   0,       3339.0,       0.0,    -13.0,      -107.0,      0.0,      1.0),
    (  -1,   0,   2,  -2,   1,      -1987.0,       0.0,     -6.0,      1073.0,      0.0,     -2.0),
    (   1,   0,   0,   0,   2,      -1981.0,       0.0,      0.0,       854.0,      0.0,      0.0),
    (  -1,   0,   0,   1,   0,       4026.0,       0.0,   -353.0,      -553.0,      0.0,   -139.0),
    (   0,   0,   2,   1,   2,       1660.0,       0.0,     -5.0,      -710.0,      0.0,     -2.0),
    (  -1,   0,   2,   4,   2,      -1521.0,       0.0,      9.0,       647.0,      0.0,      4.0),
    (  -1,   1,   0,   1,   1,       1314.0,       0.0,      0.0,      -700.0,      0.0,      0.0),
    (   0,  -2,   2,  -2,   1,      -1283.0,       0.0,      0.0,       672.0,      0.0,      0.0),
    (   1,   0,   2,   2,   1,      -1331.0,       0.0,      8.0,       663.0,      0.0,      4.0),
    (  -2,   0,   2,   2,   2,       1383.0,       0.0,     -2.0,      -594.0,      0.0,     -2.0),
    (  -1,   0,   0,   0,   2,       1405.0,       0.0,      4.0,      -610.0,      0.0,      2.0),
    (   1,   1,   2,  -2,   2,       1290.0,       0.0,      0.0,      -556.0,      0.0,      0.0),
];

// ============================================================
// 角度ユーティリティ
// ============================================================

/// 角度を [0, 360) に正規化する（純粋関数）
#[inline]
pub fn norm_angle(deg: f64) -> f64 {
    ((deg % 360.0) + 360.0) % 360.0
}

// ============================================================
// IAU 2006 平均黄道傾斜角
// ============================================================

/// IAU 2006 平均黄道傾斜角 ε_A を返す [度]（純粋関数）
///
/// Capitaine et al. 2003 / coordinates.py obliquity() と同一係数。
#[inline]
pub fn obliquity(jd: f64) -> f64 {
    let t  = (jd - J2000_JD) / JULIAN_CENTURY;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;
    let arcsec = 84381.406
        - 46.836769 * t
        - 0.0001831 * t2
        + 0.00200340 * t3
        - 0.000000576 * t4
        - 0.0000000434 * t5;
    arcsec / 3600.0
}

// ============================================================
// IAU 2000B 章動（77 項）
// ============================================================

/// IAU 2000B 5 基本引数 (l, l', F, D, Ω) を返す [rad]（内部用）
///
/// T = J2000.0 からのユリウス世紀数。coordinates.py _nut_fund_args(T) と同一。
#[inline]
fn nut_fund_args(t: f64) -> (f64, f64, f64, f64, f64) {
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let tpow = [1.0_f64, t, t2, t3, t4];

    let mut args = [0.0_f64; 5];
    for i in 0..5 {
        let arcsec: f64 = NUT_FUND_COEFFS[i]
            .iter()
            .zip(tpow.iter())
            .map(|(c, p)| c * p)
            .sum();
        let rad = (arcsec / 3600.0) * D2R;
        args[i] = ((rad % TWO_PI) + TWO_PI) % TWO_PI;
    }
    (args[0], args[1], args[2], args[3], args[4])
}

/// IAU 2000B 章動角 (ΔΨ, Δε) を返す [arcsec]（純粋関数）
///
/// IERS Conventions 2010 の 77 項シリーズを使用。
/// coordinates.py nutation_angles(jd) と同一アルゴリズム。
#[inline]
pub fn nutation_angles(jd: f64) -> (f64, f64) {
    let t = (jd - J2000_JD) / JULIAN_CENTURY;
    let (l, lp, f_, d, om) = nut_fund_args(t);

    let mut dpsi = 0.0_f64;
    let mut deps = 0.0_f64;

    for (nl, nlp, nf, nd, nom, aa, bb, cc, dd, ee, ff) in &NUT77 {
        let arg = (*nl as f64) * l
            + (*nlp as f64) * lp
            + (*nf as f64) * f_
            + (*nd as f64) * d
            + (*nom as f64) * om;
        let s = arg.sin();
        let c = arg.cos();
        dpsi += (aa + bb * t) * s + cc * c;
        deps += (dd + ee * t) * c + ff * s;
    }
    (dpsi * 1e-7, deps * 1e-7)
}

// ============================================================
// ICRS XYZ → J2000.0 黄道球面座標
// ============================================================

/// ICRS XYZ → J2000.0 黄道球面座標（ε₀ 回転のみ）（純粋関数）
///
/// coordinates.py icrs_to_j2000_ecliptic() と同一。
/// Returns: (lon_deg, lat_deg, dist)  lon/lat は度, dist は入力と同単位
pub fn icrs_to_j2000_ecliptic(x: f64, y: f64, z: f64) -> (f64, f64, f64) {
    let eps0_rad = EPS0_J2000_ARCSEC * A2R;
    let ce = eps0_rad.cos();
    let se = eps0_rad.sin();

    let xe = x;
    let ye =  ce * y + se * z;
    let ze = -se * y + ce * z;

    let dist = (xe*xe + ye*ye + ze*ze).sqrt();
    if dist == 0.0 {
        return (0.0, 0.0, 0.0);
    }
    let lat = (ze / dist).clamp(-1.0, 1.0).asin() * R2D;
    let lon = norm_angle(ye.atan2(xe) * R2D);
    (lon, lat, dist)
}

// ============================================================
// 内部ヘルパー: 歳差 + 章動 + 黄道変換
// ============================================================

/// IAU 2006 歳差 + IAU 2000B 章動 + R1(ε_true) を一括適用（内部用）
///
/// 黄道直交座標 (xl, yl, zl) を返す。
/// icrs_to_ecliptic() と _icrs_to_ecliptic_xyz_ofdate() の共通実装。
///
/// 最適化: obliquity / nutation_angles を各 1 回のみ呼ぶ。
/// （Python 版は icrs_to_equatorial_xyz_ofdate → _icrs_to_ecliptic_xyz_ofdate で各 2 回）
#[inline]
fn ecliptic_xyz_ofdate_inner(x: f64, y: f64, z: f64, jd: f64) -> (f64, f64, f64) {
    let t  = (jd - J2000_JD) / JULIAN_CENTURY;
    let t2 = t * t;
    let t3 = t2 * t;
    let t4 = t3 * t;
    let t5 = t4 * t;

    let eps0a = EPS0_J2000_ARCSEC * A2R;

    let psi_a = (PREC_PSI_A[0]*t  + PREC_PSI_A[1]*t2 + PREC_PSI_A[2]*t3
               + PREC_PSI_A[3]*t4 + PREC_PSI_A[4]*t5
               + PSI_BIAS_ARCSEC) * A2R;

    let omg_a = (EPS0_J2000_ARCSEC
               + PREC_OMG_A_POLY[0]*t  + PREC_OMG_A_POLY[1]*t2
               + PREC_OMG_A_POLY[2]*t3 + PREC_OMG_A_POLY[3]*t4
               + PREC_OMG_A_POLY[4]*t5) * A2R;

    let chi_a = (PREC_CHI_A[0]*t  + PREC_CHI_A[1]*t2 + PREC_CHI_A[2]*t3
               + PREC_CHI_A[3]*t4 + PREC_CHI_A[4]*t5) * A2R;

    let eps_a = obliquity(jd) * D2R;

    // ① R₁(ε₀): ICRS → 歳差前黄道
    let ce0 = eps0a.cos(); let se0 = eps0a.sin();
    let ay =  ce0 * y + se0 * z;
    let az = -se0 * y + ce0 * z;

    // ② R₃(−ψ_A): z 軸まわり −ψ_A 回転
    let cpsi = psi_a.cos(); let spsi = psi_a.sin();
    let bx =  cpsi * x  - spsi * ay;
    let by =  spsi * x  + cpsi * ay;

    // ③ R₁(−ω_A): x 軸まわり −ω_A 回転
    let comg = omg_a.cos(); let somg = omg_a.sin();
    let cy =  comg * by - somg * az;
    let cz =  somg * by + comg * az;

    // ④ R₃(χ_A): z 軸まわり +χ_A 回転 → of-date 平均赤道系
    let cchi = chi_a.cos(); let schi = chi_a.sin();
    let xm =  cchi * bx + schi * cy;
    let ym = -schi * bx + cchi * cy;
    let zm =  cz;

    // ⑤ 章動行列 N（1 次近似）
    let (dpsi_as, deps_as) = nutation_angles(jd);
    let dpsi = dpsi_as * A2R;
    let deps = deps_as * A2R;
    let cea = eps_a.cos(); let sea = eps_a.sin();

    let xtr = xm - dpsi * (cea * ym + sea * zm);
    let ytr = dpsi * cea * xm + ym - deps * zm;
    let ztr = dpsi * sea * xm + deps * ym + zm;

    // ⑥ R₁(ε_true): 真赤道 → 真黄道
    let eps_true = eps_a + deps;
    let se = eps_true.sin(); let ce = eps_true.cos();

    let xl = xtr;
    let yl = ytr * ce + ztr * se;
    let zl = -ytr * se + ztr * ce;

    (xl, yl, zl)
}

// ============================================================
// ICRS XYZ → of-date 真黄道球面座標
// ============================================================

/// ICRS XYZ → of-date 真黄道球面座標（純粋関数）
///
/// IAU 2006 歳差 + IAU 2000B 77 項章動 + PSI_BIAS 補正を適用。
/// coordinates.py icrs_to_ecliptic() と同一アルゴリズム。
/// 精度: JPL Horizons との黄経差 < 0.002"（2000〜2100 年）
///
/// Returns: (lon_deg, lat_deg, dist)  lon/lat は度, dist は入力と同単位
pub fn icrs_to_ecliptic(x: f64, y: f64, z: f64, jd: f64) -> (f64, f64, f64) {
    let (xl, yl, zl) = ecliptic_xyz_ofdate_inner(x, y, z, jd);

    let dist = (xl*xl + yl*yl + zl*zl).sqrt();
    if dist == 0.0 {
        return (0.0, 0.0, 0.0);
    }
    let lat = (zl / dist).clamp(-1.0, 1.0).asin() * R2D;
    let lon = norm_angle(yl.atan2(xl) * R2D);
    (lon, lat, dist)
}

// ============================================================
// 年周光行差補正（速度ベクトル法）
// ============================================================

/// 速度ベクトル法による年周光行差補正（ICRS 空間、相対論的 1 次近似）
///
/// u' = (u + β) / (1 + u·β)
/// coordinates.py apply_aberration() と同一数式。
/// 精度: < 0.001"（速度ベクトル法）
///
/// Returns: 光行差補正済み ICRS 単位ベクトル（無次元）
pub fn apply_aberration(
    ax: f64, ay: f64, az: f64,
    vx: f64, vy: f64, vz: f64,
) -> (f64, f64, f64) {
    let dist = (ax*ax + ay*ay + az*az).sqrt();
    let ux = ax / dist;
    let uy = ay / dist;
    let uz = az / dist;

    let bx = vx / C_KM_PER_DAY;
    let by = vy / C_KM_PER_DAY;
    let bz = vz / C_KM_PER_DAY;

    let udotb = ux*bx + uy*by + uz*bz;
    let inv = 1.0 / (1.0 + udotb);

    ((ux + bx) * inv, (uy + by) * inv, (uz + bz) * inv)
}

// ============================================================
// 光偏差補正（太陽重力場 / 相対論的 1 次近似）
// ============================================================

/// 太陽重力場による光偏差補正（相対論的 1 次近似）
///
/// Δê = [2GM/(c²r)] × [ê_q − ê_e(ê_q·ê_e)] / (1 + ê_q·ê_e)
/// coordinates.py apply_light_deflection() と同一数式。
///
/// Returns: 光偏差補正済み単位ベクトル（無次元）
pub fn apply_light_deflection(
    ax: f64, ay: f64, az: f64,
    sun_x: f64, sun_y: f64, sun_z: f64,
) -> (f64, f64, f64) {
    let dist_e = (ax*ax + ay*ay + az*az).sqrt();
    let ex = ax / dist_e;
    let ey = ay / dist_e;
    let ez = az / dist_e;

    let dist_q = (sun_x*sun_x + sun_y*sun_y + sun_z*sun_z).sqrt();
    let qx = sun_x / dist_q;
    let qy = sun_y / dist_q;
    let qz = sun_z / dist_q;

    let coeff = DEFL_CONST_KM / dist_q;
    let qdote = qx*ex + qy*ey + qz*ez;
    let denom = 1.0 + qdote;

    let dx = coeff * (qx - ex * qdote) / denom;
    let dy = coeff * (qy - ey * qdote) / denom;
    let dz = coeff * (qz - ez * qdote) / denom;

    (ex + dx, ey + dy, ez + dz)
}

// ============================================================
// 内部テスト（cargo test で実行）
// ============================================================

#[cfg(test)]
mod tests {
    use super::*;

    const J2000: f64 = 2451545.0;
    const TOL: f64 = 1e-6;

    #[test]
    fn test_obliquity_j2000() {
        // J2000.0 での傾斜角 ≈ 23.43927944°
        let eps = obliquity(J2000);
        assert!((eps - 23.43927944).abs() < 1e-5, "obliquity={}", eps);
    }

    #[test]
    fn test_nutation_angles_j2000() {
        // erfa_constants.py VERIFICATION 値: dpsi=-13.93166389", deps=-5.76941708"
        let (dpsi, deps) = nutation_angles(J2000);
        assert!((dpsi - (-13.93166389)).abs() < 0.001, "dpsi={}", dpsi);
        assert!((deps - (-5.76941708)).abs() < 0.001, "deps={}", deps);
    }

    #[test]
    fn test_norm_angle() {
        assert!((norm_angle(0.0) - 0.0).abs() < TOL);
        assert!((norm_angle(360.0) - 0.0).abs() < TOL);
        assert!((norm_angle(-10.0) - 350.0).abs() < TOL);
        assert!((norm_angle(370.0) - 10.0).abs() < TOL);
    }

    #[test]
    fn test_apply_aberration_unit() {
        // ゼロ速度 → 入力方向をそのまま返す（正規化済み）
        let (ax, ay, az) = (1.0, 0.0, 0.0);
        let (rx, ry, rz) = apply_aberration(ax, ay, az, 0.0, 0.0, 0.0);
        assert!((rx - 1.0).abs() < TOL, "rx={}", rx);
        assert!(ry.abs() < TOL);
        assert!(rz.abs() < TOL);
    }
}
