// bsp_test.rs — 統合テスト（de440s.bsp を使った実測検証）
//
// 実行方法:
//   BSP_PATH="/Users/takashi/Stella_series/data/catalogs/de440s.bsp" cargo test
//
// BSP_PATH が未設定の場合はファイル依存テストをスキップする。
// ファイルなしでも実行できるテスト（load API 確認など）は常に動作する。

use stella_bsp_reader::{BspFile, naif};
use std::path::Path;

const AU_KM: f64 = 149_597_870.7;
const J2000: f64 = 2_451_545.0;

// ─────────────────────────────────────────────────────────────────────────
// ヘルパー
// ─────────────────────────────────────────────────────────────────────────

/// BSP_PATH 環境変数からファイルを取得する。
/// 未設定またはファイルが存在しない場合は None を返す（テストはスキップ）。
fn load_bsp() -> Option<BspFile> {
    let path = std::env::var("BSP_PATH").ok()?;
    BspFile::load(Path::new(&path)).ok()
}

fn dist(pos: [f64; 3]) -> f64 {
    (pos[0].powi(2) + pos[1].powi(2) + pos[2].powi(2)).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────
// seek 方式でのロード確認（BSP_PATH 不要・常に実行）
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn test_load_api_exists() {
    // load() が存在し、存在しないパスに対して Io エラーを返すことを確認する
    let result = BspFile::load(Path::new("/nonexistent/path/dummy.bsp"));
    assert!(result.is_err(), "存在しないファイルはエラーになること");
}

// ─────────────────────────────────────────────────────────────────────────
// ファイル依存テスト（BSP_PATH 環境変数が必要）
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn test_load_and_segments() {
    let Some(bsp) = load_bsp() else { return; };
    assert!(!bsp.segments.is_empty(), "セグメントが 0 件");
    println!("ファイル名: {}", bsp.name);
    println!("セグメント数: {}", bsp.segments.len());
    for seg in &bsp.segments {
        println!(
            "  target={:4} center={:4} type={} start={:.1} end={:.1}",
            seg.target, seg.center, seg.spk_type, seg.start_sec, seg.end_sec
        );
    }
    // セグメントすべてで end > start を確認
    for seg in &bsp.segments {
        assert!(seg.end_sec > seg.start_sec,
            "target={} のセグメントで end_sec <= start_sec", seg.target);
    }
}

// ─────────────────────────────────────────────────────────────────────────
// 位置精度テスト
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn test_get_position_sun_j2000() {
    let Some(bsp) = load_bsp() else { return; };
    let pos = bsp.compute_position(naif::SUN, naif::SSB, J2000)
        .expect("太陽位置の取得失敗");
    println!("Sun @ J2000.0 (km): x={:.3} y={:.3} z={:.3}", pos[0], pos[1], pos[2]);
    // 太陽の SSB 基準位置は概ね ±数百万 km 以内
    assert!(dist(pos) < 2_000_000.0, "太陽位置が異常: {:.0} km", dist(pos));
}

#[test]
fn test_get_position_earth_j2000() {
    let Some(bsp) = load_bsp() else { return; };
    let pos = bsp.compute_position(naif::EARTH, naif::SSB, J2000)
        .expect("地球位置の取得失敗");
    println!("Earth @ J2000.0 (km): x={:.3} y={:.3} z={:.3}", pos[0], pos[1], pos[2]);
    let d = dist(pos);
    assert!(d > 0.9 * AU_KM && d < 1.1 * AU_KM,
        "地球-SSB 距離が異常: {:.3} AU", d / AU_KM);
}

/// 地球-太陽間距離が ~1AU になることを確認する
#[test]
fn test_earth_sun_distance_is_approx_1au() {
    let Some(bsp) = load_bsp() else { return; };
    let pos = bsp.compute_position(naif::SUN, naif::EARTH, J2000)
        .expect("地球中心 Sun 位置の取得失敗");
    let d = dist(pos);
    println!("Earth-Sun distance @ J2000.0: {:.4} AU", d / AU_KM);
    // J2000.0 付近は近日点後1ヶ月程度: ~0.983 AU
    assert!(d > 0.97 * AU_KM && d < 1.02 * AU_KM,
        "地球-太陽距離が 1 AU から外れ過ぎ: {:.4} AU", d / AU_KM);
}

#[test]
fn test_same_body_returns_zero() {
    let Some(bsp) = load_bsp() else { return; };
    let pos = bsp.compute_position(naif::EARTH, naif::EARTH, J2000)
        .expect("同一天体の取得失敗");
    assert_eq!(pos, [0.0, 0.0, 0.0]);
}

// ─────────────────────────────────────────────────────────────────────────
// セグメント境界テスト
// ─────────────────────────────────────────────────────────────────────────

/// カバレッジ開始直後の JD で計算が通ること
#[test]
fn test_position_at_coverage_start() {
    let Some(bsp) = load_bsp() else { return; };
    // de440s.bsp の最も早いセグメント開始時刻を取得
    let min_start_sec = bsp.segments.iter()
        .map(|s| s.start_sec)
        .fold(f64::INFINITY, f64::min);
    let jd_start = 2_451_545.0 + min_start_sec / 86400.0 + 0.5;  // 境界 + 0.5 日
    let result = bsp.compute_position(naif::SUN, naif::SSB, jd_start);
    assert!(result.is_ok(), "カバレッジ開始直後に計算が通ること: {:?}", result);
}

/// カバレッジ終了直前の JD で計算が通ること
#[test]
fn test_position_at_coverage_end() {
    let Some(bsp) = load_bsp() else { return; };
    let max_end_sec = bsp.segments.iter()
        .map(|s| s.end_sec)
        .fold(f64::NEG_INFINITY, f64::max);
    let jd_end = 2_451_545.0 + max_end_sec / 86400.0 - 0.5;  // 境界 - 0.5 日
    let result = bsp.compute_position(naif::SUN, naif::SSB, jd_end);
    assert!(result.is_ok(), "カバレッジ終了直前に計算が通ること: {:?}", result);
}

// ─────────────────────────────────────────────────────────────────────────
// 範囲外クエリのエラー確認
// ─────────────────────────────────────────────────────────────────────────

/// カバレッジ範囲より前の JD はエラーを返す
#[test]
fn test_out_of_range_before_coverage() {
    let Some(bsp) = load_bsp() else { return; };
    // de440s.bsp は AD 1849 頃から開始。BC 1000 は確実に範囲外
    let jd_ancient = 1_356_001.0;  // BC 1000 頃
    let result = bsp.compute_position(naif::SUN, naif::SSB, jd_ancient);
    assert!(result.is_err(),
        "カバレッジ外（遠過去）はエラーになること。実際: {:?}", result);
}

/// カバレッジ範囲より後の JD はエラーを返す
#[test]
fn test_out_of_range_after_coverage() {
    let Some(bsp) = load_bsp() else { return; };
    // de440s.bsp は AD 2150 頃まで。AD 3000 は確実に範囲外
    let jd_future = 2_816_788.0;  // AD 3000 頃
    let result = bsp.compute_position(naif::SUN, naif::SSB, jd_future);
    assert!(result.is_err(),
        "カバレッジ外（遠未来）はエラーになること。実際: {:?}", result);
}

// ─────────────────────────────────────────────────────────────────────────
// 速度テスト
// ─────────────────────────────────────────────────────────────────────────

#[test]
fn test_position_and_velocity_sun() {
    let Some(bsp) = load_bsp() else { return; };
    let (pos, vel) = bsp.get_position_and_velocity(naif::SUN, naif::SSB, J2000)
        .expect("太陽位置・速度の取得失敗");
    println!("Sun vel (km/day): vx={:.6} vy={:.6} vz={:.6}", vel[0], vel[1], vel[2]);

    // pos は compute_position と一致すること
    let pos2 = bsp.compute_position(naif::SUN, naif::SSB, J2000).unwrap();
    assert!((pos[0] - pos2[0]).abs() < 1e-6, "位置が一致しない");

    // 太陽の SSB 基準速度: 惑星重力に引っ張られた分（実測値: ~1400 km/day 程度）
    let speed = dist(vel);
    assert!(speed > 0.0 && speed < 5_000.0,
        "太陽速度の大きさが異常: {:.3} km/day", speed);
}

#[test]
fn test_position_and_velocity_earth_relative_to_emb() {
    let Some(bsp) = load_bsp() else { return; };
    // de440s.bsp では Earth(399) は EMB(3) を center とした直接セグメントを持つ
    let (pos, vel) = bsp.get_position_and_velocity(naif::EARTH, naif::EMB, J2000)
        .expect("地球(EMB基準)の位置・速度の取得失敗");
    println!("Earth/EMB vel (km/day): vx={:.3} vy={:.3} vz={:.3}", vel[0], vel[1], vel[2]);

    // 地球は EMB の周りを ~27 日周期・~4670 km 半径で公転
    // 速度の大きさ ≈ 2π × 4670 / 27.3 ≈ 1075 km/day
    let speed = dist(vel);
    assert!(speed > 0.0 && speed < 5_000.0,
        "地球/EMB 速度が異常: {:.1} km/day", speed);

    // 地球-EMB 間距離は ~4700 km 以内
    let d = dist(pos);
    assert!(d < 10_000.0,
        "地球-EMB 距離が異常: {:.1} km（~4700 km 以内を期待）", d);
}

// ─────────────────────────────────────────────────────────────────────────
// 補間連続性テスト
// ─────────────────────────────────────────────────────────────────────────

/// 0.01 日ずつ変えた時に位置が連続すること（急激なジャンプがないこと）
#[test]
fn test_position_interpolation_continuity() {
    let Some(bsp) = load_bsp() else { return; };
    let mut prev_pos = bsp.compute_position(naif::SUN, naif::SSB, J2000).unwrap();
    for i in 1..=10 {
        let jd = J2000 + i as f64 * 0.01;
        let pos = bsp.compute_position(naif::SUN, naif::SSB, jd).unwrap();
        let jump = dist([pos[0]-prev_pos[0], pos[1]-prev_pos[1], pos[2]-prev_pos[2]]);
        // 太陽の 0.01 日の変位は高々 数千 km 以内
        assert!(jump < 10_000.0, "補間が不連続: step={i}, jump={:.1} km", jump);
        prev_pos = pos;
    }
}
