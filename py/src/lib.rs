//! bsp_rs — Rust 製 BSP リーダー（PyO3 バインド版）
//!
//! bsp_pure.py の Rust 翻訳。外部クレート依存は PyO3 のみ。
//! Python 公開 API は BspReader クラスと get_reader_rust() 関数。
//!
//! 対応 SPK セグメントタイプ:
//!   Type 2  — Chebyshev 多項式（位置のみ）        de440s.bsp / de440.bsp / de441.bsp
//!   Type 3  — Chebyshev 多項式（位置＋速度）       de440.bsp / de441.bsp の月秤動角セグメント
//!   Type 13 — Hermite 補間（小天体）               スコープ外（ValueError を送出）
//!
//! 設計:
//!   - DAF/SPK バイナリ解析は parse_daf() 純粋関数
//!   - Type 2/3 セグメント計算は compute_chebyshev() 純粋関数
//!   - Chebyshev 評価は cheby_eval / cheby_eval_with_deriv 純粋関数
//!   - BspReader はデータソース層。座標変換は行わない
//!
//! ライセンス: MIT

use pyo3::prelude::*;
use std::io::{Read, Seek, SeekFrom};
use std::sync::Mutex;

mod coord;

// ============================================================
// データソース（Memory: テスト・小ファイル / File: seek 方式）
// ============================================================

#[allow(dead_code)]
enum BspData {
    Memory(Vec<u8>),
    File(Mutex<std::fs::File>),
}

impl BspData {
    /// `offset` バイト目から `buf.len()` バイトを読み込む
    fn read_at(&self, offset: usize, buf: &mut [u8]) {
        match self {
            BspData::Memory(data) => {
                buf.copy_from_slice(&data[offset..offset + buf.len()]);
            }
            BspData::File(mutex) => {
                let mut f = mutex.lock().expect("BspData mutex not poisoned");
                f.seek(SeekFrom::Start(offset as u64)).expect("BSP seek failed");
                f.read_exact(buf).expect("BSP read failed");
            }
        }
    }
}

// ============================================================
// DAF / SPK フォーマット定数
// ============================================================

const RECORD_SIZE: usize = 1024;
const S_PER_DAY: f64 = 86400.0;
const J2000_JD: f64 = 2451545.0;
const SSB: i32 = 0;
const SPK_TYPE_2: i32 = 2;
const SPK_TYPE_3: i32 = 3;

// ============================================================
// セグメント記述子
// ============================================================

#[derive(Debug, Clone)]
struct Segment {
    /// セグメント開始時刻（J2000.0 からの秒数）
    start_sec: f64,
    /// セグメント終了時刻（J2000.0 からの秒数）
    end_sec: f64,
    target: i32,
    center: i32,
    spk_type: i32,
    first_addr: i32,
    last_addr: i32,
}

/// Type 2 → Some(3)（位置のみ）、Type 3 → Some(6)（位置＋速度）、それ以外 → None
fn spk_components(spk_type: i32) -> Option<usize> {
    match spk_type {
        SPK_TYPE_2 => Some(3),
        SPK_TYPE_3 => Some(6),
        _ => None,
    }
}

// ============================================================
// バイナリ読み取りヘルパー
// ============================================================

#[inline]
fn read_f64_le(data: &[u8], offset: usize) -> f64 {
    let b: [u8; 8] = data[offset..offset + 8].try_into().unwrap();
    f64::from_le_bytes(b)
}

#[inline]
fn read_i32_le(data: &[u8], offset: usize) -> i32 {
    let b: [u8; 4] = data[offset..offset + 4].try_into().unwrap();
    i32::from_le_bytes(b)
}

// ============================================================
// Chebyshev 評価（chebyshev.js / bsp_pure.py の Rust 翻訳）
// ============================================================

fn cheby_eval(coeffs: &[f64], x: f64) -> f64 {
    let n = coeffs.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return coeffs[0];
    }
    let mut b2 = 0.0_f64;
    let mut b1 = 0.0_f64;
    for i in (1..n).rev() {
        let b = coeffs[i] + 2.0 * x * b1 - b2;
        b2 = b1;
        b1 = b;
    }
    coeffs[0] + x * b1 - b2
}

fn cheby_eval_with_deriv(coeffs: &[f64], x: f64) -> (f64, f64) {
    let n = coeffs.len();
    if n == 0 {
        return (0.0, 0.0);
    }
    if n == 1 {
        return (coeffs[0], 0.0);
    }
    let mut b2 = 0.0_f64;
    let mut b1 = 0.0_f64;
    let mut d2 = 0.0_f64;
    let mut d1 = 0.0_f64;
    for i in (1..n).rev() {
        let b = coeffs[i] + 2.0 * x * b1 - b2;
        let d = 2.0 * b1 + 2.0 * x * d1 - d2;
        b2 = b1;
        b1 = b;
        d2 = d1;
        d1 = d;
    }
    let position = coeffs[0] + x * b1 - b2;
    let dpdx = b1 + x * d1 - d2;
    (position, dpdx)
}

// ============================================================
// DAF バイナリ解析（bsp_pure.py の _parse_daf 翻訳）
// ============================================================

#[allow(dead_code)]
fn parse_daf(data: &[u8]) -> Result<Vec<Segment>, String> {
    let locidw = std::str::from_utf8(&data[0..8]).unwrap_or("");
    if !locidw.starts_with("DAF/SPK") && !locidw.starts_with("DAF/EK") {
        return Err(format!("非 SPK ファイルです。LOCIDW='{}'", locidw));
    }

    let nd = read_i32_le(data, 8) as usize;
    let ni = read_i32_le(data, 12) as usize;
    let first_sum_rec = read_i32_le(data, 76) as usize;

    // math.ceil(ni / 2) の整数演算版
    let summary_doubles = nd + (ni + 1) / 2;
    let summary_bytes = summary_doubles * 8;

    let mut segments = Vec::new();
    let mut rec_num = first_sum_rec;

    while rec_num > 0 {
        let rec_offset = (rec_num - 1) * RECORD_SIZE;
        let next_rec = read_f64_le(data, rec_offset).round() as usize;
        let n_summaries = read_f64_le(data, rec_offset + 16).round() as usize;

        for i in 0..n_summaries {
            let off = rec_offset + 24 + i * summary_bytes;
            let int_off = off + nd * 8;
            segments.push(Segment {
                start_sec:  read_f64_le(data, off),
                end_sec:    read_f64_le(data, off + 8),
                target:     read_i32_le(data, int_off),
                center:     read_i32_le(data, int_off + 4),
                spk_type:   read_i32_le(data, int_off + 12),
                first_addr: read_i32_le(data, int_off + 16),
                last_addr:  read_i32_le(data, int_off + 20),
            });
        }

        rec_num = next_rec;
    }

    Ok(segments)
}

/// seek 方式でサマリーレコードのみを読み込む（大ファイル対応）
fn parse_summaries_from_file(
    file: &mut std::fs::File,
    nd: usize,
    ni: usize,
    first_sum_rec: usize,
) -> Result<Vec<Segment>, String> {
    let summary_doubles = nd + (ni + 1) / 2;
    let summary_bytes = summary_doubles * 8;

    let mut segments = Vec::new();
    let mut rec_num = first_sum_rec;

    while rec_num > 0 {
        let rec_offset = ((rec_num - 1) * RECORD_SIZE) as u64;
        let mut rec_buf = [0u8; RECORD_SIZE];
        file.seek(SeekFrom::Start(rec_offset))
            .map_err(|e| e.to_string())?;
        file.read_exact(&mut rec_buf)
            .map_err(|e| e.to_string())?;

        let next_rec = read_f64_le(&rec_buf, 0).round() as usize;
        let n_summaries = read_f64_le(&rec_buf, 16).round() as usize;

        for i in 0..n_summaries {
            let off = 24 + i * summary_bytes;
            let int_off = off + nd * 8;
            segments.push(Segment {
                start_sec:  read_f64_le(&rec_buf, off),
                end_sec:    read_f64_le(&rec_buf, off + 8),
                target:     read_i32_le(&rec_buf, int_off),
                center:     read_i32_le(&rec_buf, int_off + 4),
                spk_type:   read_i32_le(&rec_buf, int_off + 12),
                first_addr: read_i32_le(&rec_buf, int_off + 16),
                last_addr:  read_i32_le(&rec_buf, int_off + 20),
            });
        }

        rec_num = next_rec;
    }

    Ok(segments)
}

// ============================================================
// BspReader クラス（PyO3 公開）
// Type 2/3 計算は BspReader のプライベートメソッドとして実装（seek 対応）
// ============================================================

/// Rust 製 BSP リーダー。bsp_pure.py と同一 Python API を持つ。
///
/// compute_position(target, center, jd_tdb) -> list[float]
/// compute_position_and_velocity(target, center, jd_tdb) -> (list[float], list[float])
#[pyclass]
pub struct BspReader {
    data: BspData,
    bsp_path: String,
    segments: Vec<Segment>,
}

impl BspReader {
    /// `offset` バイト目の f64 を読む（Memory/File 両対応）
    fn read_f64_at(&self, offset: usize) -> f64 {
        let mut buf = [0u8; 8];
        self.data.read_at(offset, &mut buf);
        f64::from_le_bytes(buf)
    }

    /// Type 2（components=3）および Type 3（components=6）共通の Chebyshev 評価。
    /// jd_tdb がカバー範囲外の場合は PyValueError を返す。
    fn compute_chebyshev_seg(
        &self,
        seg: &Segment,
        jd_tdb: f64,
        with_velocity: bool,
        components: usize,
    ) -> PyResult<([f64; 3], Option<[f64; 3]>)> {
        let data_start = (seg.first_addr as usize - 1) * 8;
        let data_end = seg.last_addr as usize * 8;
        let meta_offset = data_end - 32;

        let init   = self.read_f64_at(meta_offset);
        let intlen = self.read_f64_at(meta_offset + 8);
        let rsize  = self.read_f64_at(meta_offset + 16).round() as usize;
        let n      = self.read_f64_at(meta_offset + 24).round() as i64;

        let t_sec = (jd_tdb - J2000_JD) * S_PER_DAY;

        // 範囲外チェック（サイレントクランプ禁止）
        let t_min = init;
        let t_max = init + n as f64 * intlen;
        if t_sec < t_min || t_sec > t_max {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "out of coverage: t={:.1}s (valid {:.1}s – {:.1}s) for target={}",
                t_sec, t_min, t_max, seg.target
            )));
        }

        let idx = ((t_sec - init) / intlen) as i64;
        let idx = idx.clamp(0, n - 1) as usize;

        let rec_offset = data_start + idx * rsize * 8;
        let mid    = self.read_f64_at(rec_offset);
        let radius = self.read_f64_at(rec_offset + 8);

        let x = (t_sec - mid) / radius;

        let ncoeff = (rsize - 2) / components;
        let base = rec_offset + 16;

        let cx: Vec<f64> = (0..ncoeff).map(|k| self.read_f64_at(base + k * 8)).collect();
        let cy: Vec<f64> = (0..ncoeff).map(|k| self.read_f64_at(base + ncoeff * 8 + k * 8)).collect();
        let cz: Vec<f64> = (0..ncoeff).map(|k| self.read_f64_at(base + ncoeff * 8 * 2 + k * 8)).collect();

        if with_velocity {
            if components == 6 {
                // Type 3: 速度係数を直接評価（km/s → km/day に変換）
                let vx: Vec<f64> = (0..ncoeff).map(|k| self.read_f64_at(base + ncoeff * 8 * 3 + k * 8)).collect();
                let vy: Vec<f64> = (0..ncoeff).map(|k| self.read_f64_at(base + ncoeff * 8 * 4 + k * 8)).collect();
                let vz: Vec<f64> = (0..ncoeff).map(|k| self.read_f64_at(base + ncoeff * 8 * 5 + k * 8)).collect();
                Ok((
                    [cheby_eval(&cx, x), cheby_eval(&cy, x), cheby_eval(&cz, x)],
                    Some([
                        cheby_eval(&vx, x) * S_PER_DAY,
                        cheby_eval(&vy, x) * S_PER_DAY,
                        cheby_eval(&vz, x) * S_PER_DAY,
                    ]),
                ))
            } else {
                // Type 2: 位置多項式を微分して速度を求める
                let (px, dpx) = cheby_eval_with_deriv(&cx, x);
                let (py, dpy) = cheby_eval_with_deriv(&cy, x);
                let (pz, dpz) = cheby_eval_with_deriv(&cz, x);
                let dxdt = S_PER_DAY / radius;
                Ok((
                    [px, py, pz],
                    Some([dpx * dxdt, dpy * dxdt, dpz * dxdt]),
                ))
            }
        } else {
            Ok((
                [cheby_eval(&cx, x), cheby_eval(&cy, x), cheby_eval(&cz, x)],
                None,
            ))
        }
    }
}

#[pymethods]
impl BspReader {
    #[new]
    pub fn new(bsp_path: String) -> PyResult<Self> {
        let mut file = std::fs::File::open(&bsp_path)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

        // ファイルレコード（1024 バイト）だけ読んで nd/ni/first_sum_rec を取得
        let mut header = [0u8; RECORD_SIZE];
        file.read_exact(&mut header)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;

        let locidw = std::str::from_utf8(&header[0..8]).unwrap_or("");
        if !locidw.starts_with("DAF/SPK") && !locidw.starts_with("DAF/EK") {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "非 SPK ファイルです。LOCIDW='{}'", locidw
            )));
        }
        let nd = read_i32_le(&header, 8) as usize;
        let ni = read_i32_le(&header, 12) as usize;
        let first_sum_rec = read_i32_le(&header, 76) as usize;

        let segments = parse_summaries_from_file(&mut file, nd, ni, first_sum_rec)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

        Ok(BspReader {
            data: BspData::File(Mutex::new(file)),
            bsp_path,
            segments,
        })
    }

    /// ICRS 位置ベクトルを返す [km]
    pub fn compute_position(
        &self, target: i32, center: i32, jd_tdb: f64,
    ) -> PyResult<Vec<f64>> {
        if target == center {
            return Ok(vec![0.0, 0.0, 0.0]);
        }
        if let Some(seg) = self.find_segment(target, center, jd_tdb) {
            let components = spk_components(seg.spk_type).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "unsupported SPK type: {}", seg.spk_type
                ))
            })?;
            let (pos, _) = self.compute_chebyshev_seg(seg, jd_tdb, false, components)?;
            return Ok(pos.to_vec());
        }
        let pos_t = self.pos_from_ssb(target, jd_tdb)?;
        if center == SSB {
            return Ok(pos_t);
        }
        let pos_c = self.pos_from_ssb(center, jd_tdb)?;
        Ok(vec![
            pos_t[0] - pos_c[0],
            pos_t[1] - pos_c[1],
            pos_t[2] - pos_c[2],
        ])
    }

    /// ICRS 位置・速度ベクトルを返す (pos [km], vel [km/day])
    pub fn compute_position_and_velocity(
        &self, target: i32, center: i32, jd_tdb: f64,
    ) -> PyResult<(Vec<f64>, Vec<f64>)> {
        if target == center {
            return Ok((vec![0.0, 0.0, 0.0], vec![0.0, 0.0, 0.0]));
        }
        if let Some(seg) = self.find_segment(target, center, jd_tdb) {
            let components = spk_components(seg.spk_type).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "unsupported SPK type: {}", seg.spk_type
                ))
            })?;
            let (pos, vel) = self.compute_chebyshev_seg(seg, jd_tdb, true, components)?;
            return Ok((pos.to_vec(), vel.unwrap().to_vec()));
        }
        let (pos_t, vel_t) = self.pos_vel_from_ssb(target, jd_tdb)?;
        if center == SSB {
            return Ok((pos_t, vel_t));
        }
        let (pos_c, vel_c) = self.pos_vel_from_ssb(center, jd_tdb)?;
        Ok((
            vec![pos_t[0] - pos_c[0], pos_t[1] - pos_c[1], pos_t[2] - pos_c[2]],
            vec![vel_t[0] - vel_c[0], vel_t[1] - vel_c[1], vel_t[2] - vel_c[2]],
        ))
    }

    /// JD 配列をまとめて渡して位置ベクトル配列を返す（バッチ API）
    ///
    /// Python 側の for ループを Rust 内タイトループに置き換えることで
    /// FFI オーバーヘッドを1回に抑え、大量サンプリングを高速化する。
    ///
    /// Returns: [[x,y,z], ...] — len == len(jd_list)
    pub fn compute_positions_batch(
        &self, target: i32, center: i32, jd_list: Vec<f64>,
    ) -> PyResult<Vec<Vec<f64>>> {
        jd_list.iter().map(|&jd| self.compute_position(target, center, jd)).collect()
    }

    /// JD 配列をまとめて渡して位置・速度ベクトル配列を返す（バッチ API）
    ///
    /// Returns: ([[x,y,z], ...], [[vx,vy,vz], ...]) — 各 len == len(jd_list)
    pub fn compute_positions_and_velocities_batch(
        &self, target: i32, center: i32, jd_list: Vec<f64>,
    ) -> PyResult<(Vec<Vec<f64>>, Vec<Vec<f64>>)> {
        let mut positions = Vec::with_capacity(jd_list.len());
        let mut velocities = Vec::with_capacity(jd_list.len());
        for &jd in &jd_list {
            let (pos, vel) = self.compute_position_and_velocity(target, center, jd)?;
            positions.push(pos);
            velocities.push(vel);
        }
        Ok((positions, velocities))
    }

    // ──────────────────────────────────────────────────────
    // compute_apparent_batch — 視位置バッチ（光行時+光偏差+光行差+座標変換）
    // ──────────────────────────────────────────────────────

    /// 視位置（apparent position）と経度・緯度速度をバッチ計算する（PyO3 公開）
    ///
    /// apparent.py compute_apparent() を Rust 化したバッチ版。
    /// 1 JD あたり内部で 3 回（位置 / +30s / -30s）計算して速度を求める。
    ///
    /// Args:
    ///   naif_target  : ターゲット NAIF コード
    ///   center_naif  : 観測中心 NAIF コード（399=地球）
    ///   jd_tdb_list  : JD（TDB）のリスト
    ///   use_j2000    : True → J2000.0 黄道（歳差・章動・光行差なし）
    ///   aberration   : True → 年周光行差を適用
    ///   deflection   : True → 光偏差（太陽重力場による偏向）を適用
    ///   light_time   : True → 光行時間 τ = r/c を適用（False → τ=0、瞬時位置）
    ///
    /// Returns: [(lon_deg, lat_deg, dist_km, lonspeed_deg_per_day, latspeed_deg_per_day), ...]
    ///
    /// deflection・light_time はキーワード省略可能（デフォルト true）。
    /// 2026/05 時点の既存呼び出し（5引数）との後方互換性を維持するため。
    #[pyo3(signature = (naif_target, center_naif, jd_tdb_list, use_j2000, aberration, deflection=true, light_time=true))]
    pub fn compute_apparent_batch(
        &self,
        naif_target: i32,
        center_naif: i32,
        jd_tdb_list: Vec<f64>,
        use_j2000: bool,
        aberration: bool,
        deflection: bool,
        light_time: bool,
    ) -> PyResult<Vec<(f64, f64, f64, f64, f64)>> {
        let dt = coord::SPEED_DT_DAYS;
        let two_dt = 2.0 * dt;
        let mut results = Vec::with_capacity(jd_tdb_list.len());

        for &jd in &jd_tdb_list {
            let (lon, lat, dist) =
                self.apparent_single(naif_target, center_naif, jd, use_j2000, aberration, deflection, light_time)?;
            let (lon_p, lat_p, _) =
                self.apparent_single(naif_target, center_naif, jd + dt, use_j2000, aberration, deflection, light_time)?;
            let (lon_m, lat_m, _) =
                self.apparent_single(naif_target, center_naif, jd - dt, use_j2000, aberration, deflection, light_time)?;

            let lon_p = lon_p % 360.0;
            let lon_m = lon_m % 360.0;
            let mut raw = lon_p - lon_m;
            if raw > 180.0 { raw -= 360.0; } else if raw < -180.0 { raw += 360.0; }
            let lonspeed = raw / two_dt;
            let latspeed = (lat_p - lat_m) / two_dt;

            results.push((lon % 360.0, lat, dist, lonspeed, latspeed));
        }
        Ok(results)
    }

    /// 任意重心からの視位置と速度をバッチ計算する（PyO3 公開）
    ///
    /// apparent.py compute_from_center() の Rust バッチ版。
    /// 光偏差補正なし・年周光行差のみ（Skyfield deflectors=[] 相当）。
    ///
    /// Returns: [(lon_deg, lat_deg, dist_km, lonspeed_deg_per_day, latspeed_deg_per_day), ...]
    pub fn compute_from_center_batch(
        &self,
        naif_target: i32,
        center_naif: i32,
        jd_tdb_list: Vec<f64>,
        use_j2000: bool,
        aberration: bool,
    ) -> PyResult<Vec<(f64, f64, f64, f64, f64)>> {
        if naif_target == center_naif {
            return Ok(vec![(0.0, 0.0, 0.0, 0.0, 0.0); jd_tdb_list.len()]);
        }
        let dt = coord::SPEED_DT_DAYS;
        let two_dt = 2.0 * dt;
        let mut results = Vec::with_capacity(jd_tdb_list.len());

        for &jd in &jd_tdb_list {
            let (lon, lat, dist) =
                self.from_center_single(naif_target, center_naif, jd, use_j2000, aberration)?;
            let (lon_p, lat_p, _) =
                self.from_center_single(naif_target, center_naif, jd + dt, use_j2000, aberration)?;
            let (lon_m, lat_m, _) =
                self.from_center_single(naif_target, center_naif, jd - dt, use_j2000, aberration)?;

            let lon_p = lon_p % 360.0;
            let lon_m = lon_m % 360.0;
            let mut raw = lon_p - lon_m;
            if raw > 180.0 { raw -= 360.0; } else if raw < -180.0 { raw += 360.0; }
            let lonspeed = raw / two_dt;
            let latspeed = (lat_p - lat_m) / two_dt;

            results.push((lon % 360.0, lat, dist, lonspeed, latspeed));
        }
        Ok(results)
    }

    pub fn close(&mut self) {
        self.segments.clear();
    }

    pub fn __repr__(&self) -> String {
        format!("BspReader(rust, {:?})", self.bsp_path)
    }
}

impl BspReader {
    // ──────────────────────────────────────────────────────────────────────
    // 内部ヘルパー: apparent / from_center の単一 JD 計算
    // ──────────────────────────────────────────────────────────────────────

    /// compute_apparent の 1 JD 版（Rust 内部専用）
    ///
    /// apparent.py compute_apparent() を忠実に Rust 化。
    /// トポセントリック補正は省略（stella_engine では不使用）。
    ///
    /// deflection と aberration は独立に ON/OFF できる
    /// （apparent.py【D】2026/05/30 の分離仕様に対応。ステップ番号も同スクリプトと対応させている）。
    /// light_time=false の場合は τ=0（瞬時位置。apparent.py【E】②-1 2026/08/03 対応）。
    fn apparent_single(
        &self,
        naif_target: i32,
        center_naif: i32,
        jd_tdb: f64,
        use_j2000: bool,
        aberration: bool,
        deflection: bool,
        light_time: bool,
    ) -> PyResult<(f64, f64, f64)> {
        // 1. 幾何学的距離 → 光行時間 τ（light_time=false なら τ=0）
        let tau = if light_time {
            let geo_pos = self.compute_position(naif_target, center_naif, jd_tdb)?;
            (geo_pos[0]*geo_pos[0]
                + geo_pos[1]*geo_pos[1]
                + geo_pos[2]*geo_pos[2]).sqrt() / coord::C_KM_PER_DAY
        } else {
            0.0
        };

        // 2. 実体位置: target(t−τ) − center(t)  [ICRS, km]
        let center_ssb = self.pos_from_ssb(center_naif, jd_tdb)?;
        let target_ssb = self.pos_from_ssb(naif_target, jd_tdb - tau)?;
        let ax = target_ssb[0] - center_ssb[0];
        let ay = target_ssb[1] - center_ssb[1];
        let az = target_ssb[2] - center_ssb[2];

        // 3. J2000.0 モード（ε₀ 回転のみ）
        if use_j2000 {
            return Ok(coord::icrs_to_j2000_ecliptic(ax, ay, az));
        }

        // 4. 両方 OFF → astrometric（光行時のみ）
        if !aberration && !deflection {
            return Ok(coord::icrs_to_ecliptic(ax, ay, az, jd_tdb));
        }

        // 5. 光偏差補正（太陽以外の天体・deflection=true のときのみ）
        let (mut bx, mut by, mut bz) = (ax, ay, az);
        if deflection && naif_target != coord::NAIF_SUN {
            let sun_ssb = self.pos_from_ssb(coord::NAIF_SUN, jd_tdb)?;
            let sun_x = sun_ssb[0] - center_ssb[0];
            let sun_y = sun_ssb[1] - center_ssb[1];
            let sun_z = sun_ssb[2] - center_ssb[2];
            let (dx, dy, dz) = coord::apply_light_deflection(ax, ay, az, sun_x, sun_y, sun_z);
            let dist0 = (ax*ax + ay*ay + az*az).sqrt();
            bx = dx * dist0;
            by = dy * dist0;
            bz = dz * dist0;
        }

        // 6. 年周光行差補正（観測中心の重心速度 ±0.5 s 有限差分・aberration=true のときのみ）
        if aberration {
            let e_plus  = self.pos_from_ssb(center_naif, jd_tdb + coord::ABERR_DT_DAYS)?;
            let e_minus = self.pos_from_ssb(center_naif, jd_tdb - coord::ABERR_DT_DAYS)?;
            let two_dt = 2.0 * coord::ABERR_DT_DAYS;
            let vx = (e_plus[0] - e_minus[0]) / two_dt;
            let vy = (e_plus[1] - e_minus[1]) / two_dt;
            let vz = (e_plus[2] - e_minus[2]) / two_dt;
            let (abr_x, abr_y, abr_z) = coord::apply_aberration(bx, by, bz, vx, vy, vz);
            let dist = (bx*bx + by*by + bz*bz).sqrt();
            bx = abr_x * dist;
            by = abr_y * dist;
            bz = abr_z * dist;
        }

        // 7. ICRS → of-date 真黄道
        Ok(coord::icrs_to_ecliptic(bx, by, bz, jd_tdb))
    }

    /// compute_from_center の 1 JD 版（Rust 内部専用）
    ///
    /// apparent.py compute_from_center() を Rust 化。
    /// 光偏差なし（Skyfield deflectors=[] 相当）。
    fn from_center_single(
        &self,
        naif_target: i32,
        center_naif: i32,
        jd_tdb: f64,
        use_j2000: bool,
        aberration: bool,
    ) -> PyResult<(f64, f64, f64)> {
        // 1. 幾何学的距離 → 光行時間 τ
        let geo_pos = self.compute_position(naif_target, center_naif, jd_tdb)?;
        let geo_dist = (geo_pos[0]*geo_pos[0]
                       + geo_pos[1]*geo_pos[1]
                       + geo_pos[2]*geo_pos[2]).sqrt();
        let tau = geo_dist / coord::C_KM_PER_DAY;

        // 2. 実体位置: target(t−τ) − center(t)  [ICRS, km]
        let center_ssb = self.pos_from_ssb(center_naif, jd_tdb)?;
        let target_ssb = self.pos_from_ssb(naif_target, jd_tdb - tau)?;
        let mut ax = target_ssb[0] - center_ssb[0];
        let mut ay = target_ssb[1] - center_ssb[1];
        let mut az = target_ssb[2] - center_ssb[2];

        // 3. 年周光行差補正（観測中心の重心速度 ±0.5 s 有限差分）
        if aberration {
            let e_plus  = self.pos_from_ssb(center_naif, jd_tdb + coord::ABERR_DT_DAYS)?;
            let e_minus = self.pos_from_ssb(center_naif, jd_tdb - coord::ABERR_DT_DAYS)?;
            let two_dt = 2.0 * coord::ABERR_DT_DAYS;
            let vx = (e_plus[0] - e_minus[0]) / two_dt;
            let vy = (e_plus[1] - e_minus[1]) / two_dt;
            let vz = (e_plus[2] - e_minus[2]) / two_dt;
            let (abr_x, abr_y, abr_z) = coord::apply_aberration(ax, ay, az, vx, vy, vz);
            let dist = (ax*ax + ay*ay + az*az).sqrt();
            ax = abr_x * dist;
            ay = abr_y * dist;
            az = abr_z * dist;
        }

        // 4. ICRS → 黄道
        if use_j2000 {
            Ok(coord::icrs_to_j2000_ecliptic(ax, ay, az))
        } else {
            Ok(coord::icrs_to_ecliptic(ax, ay, az, jd_tdb))
        }
    }

    fn find_segment(&self, target: i32, center: i32, jd_tdb: f64) -> Option<&Segment> {
        let t_sec = (jd_tdb - J2000_JD) * S_PER_DAY;
        self.segments.iter().find(|seg| {
            seg.target == target
                && seg.center == center
                && seg.start_sec <= t_sec
                && t_sec <= seg.end_sec
        })
    }

    fn pos_from_ssb(&self, target: i32, jd_tdb: f64) -> PyResult<Vec<f64>> {
        let t_sec = (jd_tdb - J2000_JD) * S_PER_DAY;
        if let Some(seg) = self.find_segment(target, SSB, jd_tdb) {
            let components = spk_components(seg.spk_type).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "unsupported SPK type: {}", seg.spk_type
                ))
            })?;
            let (pos, _) = self.compute_chebyshev_seg(seg, jd_tdb, false, components)?;
            return Ok(pos.to_vec());
        }
        for seg in &self.segments {
            if seg.target == target && seg.start_sec <= t_sec && t_sec <= seg.end_sec {
                let components = spk_components(seg.spk_type).ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err(format!(
                        "unsupported SPK type: {}", seg.spk_type
                    ))
                })?;
                let (pos_fc, _) = self.compute_chebyshev_seg(seg, jd_tdb, false, components)?;
                if seg.center == SSB {
                    return Ok(pos_fc.to_vec());
                }
                let pos_c = self.pos_from_ssb(seg.center, jd_tdb)?;
                return Ok(vec![
                    pos_c[0] + pos_fc[0],
                    pos_c[1] + pos_fc[1],
                    pos_c[2] + pos_fc[2],
                ]);
            }
        }
        Err(pyo3::exceptions::PyValueError::new_err(format!(
            "SSB から NAIF={} へのパスが BSP に存在しません", target
        )))
    }

    fn pos_vel_from_ssb(&self, target: i32, jd_tdb: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
        let t_sec = (jd_tdb - J2000_JD) * S_PER_DAY;
        if let Some(seg) = self.find_segment(target, SSB, jd_tdb) {
            let components = spk_components(seg.spk_type).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "unsupported SPK type: {}", seg.spk_type
                ))
            })?;
            let (pos, vel) = self.compute_chebyshev_seg(seg, jd_tdb, true, components)?;
            return Ok((pos.to_vec(), vel.unwrap().to_vec()));
        }
        for seg in &self.segments {
            if seg.target == target && seg.start_sec <= t_sec && t_sec <= seg.end_sec {
                let components = spk_components(seg.spk_type).ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err(format!(
                        "unsupported SPK type: {}", seg.spk_type
                    ))
                })?;
                let (p_fc, v_fc) = self.compute_chebyshev_seg(seg, jd_tdb, true, components)?;
                let v_fc = v_fc.unwrap();
                if seg.center == SSB {
                    return Ok((p_fc.to_vec(), v_fc.to_vec()));
                }
                let (p_c, v_c) = self.pos_vel_from_ssb(seg.center, jd_tdb)?;
                return Ok((
                    vec![p_c[0] + p_fc[0], p_c[1] + p_fc[1], p_c[2] + p_fc[2]],
                    vec![v_c[0] + v_fc[0], v_c[1] + v_fc[1], v_c[2] + v_fc[2]],
                ));
            }
        }
        Err(pyo3::exceptions::PyValueError::new_err(format!(
            "SSB から NAIF={} へのパスが BSP に存在しません", target
        )))
    }
}

// ============================================================
// モジュール定義
// ============================================================

#[pymodule]
fn bsp_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<BspReader>()?;
    Ok(())
}
