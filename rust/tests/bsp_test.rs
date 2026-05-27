// bsp_test.rs — 結合テスト（de440s.bsp を使った実測検証）
//
// 実行方法:
//   BSP_PATH="/Users/takashi/Stella_series/data/catalogs/de440s.bsp" cargo test
//
// BSP_PATH が未設定の場合はスキップされる。

use stella_bsp_reader::{BspFile, naif};
use std::path::Path;

/// BSP_PATH 環境変数からファイルを取得する。未設定ならテストをスキップ。
fn load_bsp() -> Option<BspFile> {
    let path = std::env::var("BSP_PATH").ok()?;
    BspFile::load(Path::new(&path)).ok()
}

#[test]
fn test_load_and_segments() {
    let Some(bsp) = load_bsp() else { return; };
    assert!(!bsp.segments.is_empty(), "セグメントが0件");
    println!("ファイル名: {}", bsp.name);
    println!("セグメント数: {}", bsp.segments.len());
    for seg in &bsp.segments {
        println!(
            "  target={:4} center={:4} type={} [{:.1}, {:.1}]",
            seg.target, seg.center, seg.spk_type, seg.start_sec, seg.end_sec
        );
    }
}

#[test]
fn test_get_position_sun_j2000() {
    // J2000.0 における太陽の SSB 基準位置
    // JPL Horizons 参照値: 太陽は SSB に近いが一致しない
    let Some(bsp) = load_bsp() else { return; };
    let jd = 2_451_545.0; // J2000.0
    let pos = bsp.compute_position(naif::SUN, naif::SSB, jd)
        .expect("太陽位置の取得失敗");
    println!("Sun @ J2000.0 (km): x={:.3} y={:.3} z={:.3}", pos[0], pos[1], pos[2]);
    // 太陽の SSB 基準位置は概ね ±数百万 km 以内
    let dist = (pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2)).sqrt();
    assert!(dist < 2_000_000.0, "太陽位置が異常: {dist:.0} km");
}

#[test]
fn test_get_position_earth_j2000() {
    // J2000.0 における地球の SSB 基準位置
    let Some(bsp) = load_bsp() else { return; };
    let jd = 2_451_545.0;
    let pos = bsp.compute_position(naif::EARTH, naif::SSB, jd)
        .expect("地球位置の取得失敗");
    println!("Earth @ J2000.0 (km): x={:.3} y={:.3} z={:.3}", pos[0], pos[1], pos[2]);
    // 地球の SSB 基準距離は概ね 1 AU（約 1.5億 km）付近
    let dist = (pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2)).sqrt();
    let au_km = 149_597_870.7;
    assert!(dist > 0.9 * au_km && dist < 1.1 * au_km,
        "地球-SSB 距離が異常: {:.3} AU", dist / au_km);
}

#[test]
fn test_same_body_returns_zero() {
    let Some(bsp) = load_bsp() else { return; };
    let pos = bsp.compute_position(naif::EARTH, naif::EARTH, 2_451_545.0)
        .expect("同一天体の取得失敗");
    assert_eq!(pos, [0.0, 0.0, 0.0]);
}

#[test]
fn test_position_and_velocity() {
    let Some(bsp) = load_bsp() else { return; };
    let jd = 2_451_545.0;
    let (pos, vel) = bsp.get_position_and_velocity(naif::SUN, naif::SSB, jd)
        .expect("太陽位置・速度の取得失敗");
    println!("Sun vel (km/day): vx={:.6} vy={:.6} vz={:.6}", vel[0], vel[1], vel[2]);
    // 速度の大きさが非ゼロであることを確認
    let speed = (vel[0].powi(2) + vel[1].powi(2) + vel[2].powi(2)).sqrt();
    assert!(speed > 0.0);
    // pos は get_position と一致すること
    let pos2 = bsp.get_position(naif::SUN, naif::SSB, jd).unwrap();
    assert_eq!(pos, pos2);
}
