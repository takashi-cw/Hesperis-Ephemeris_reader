// bsp.rs — JPL BSP/SPK バイナリ読み取り
//
// NASA NAIF DAF/SPK フォーマットを解析し、
// 指定天体の ICRS XYZ 位置ベクトル（km）を返す。
//
// 対応フォーマット:
//   - SPK Type 2（Chebyshev 多項式：位置）    ← 惑星・月・太陽
//   - SPK Type 3（Chebyshev 多項式：位置+速度）← 月秤動角、full DE440/DE441
//   ※ Type 13（Hermite 補間：小天体）は非対応（スコープ外）
//
// 出典フォーマット仕様:
//   - NAIF SPK Required Reading (NAIF N0067)
//   - NAIF DAF Required Reading (NAIF N0067)
//   - jplephem (Brandon Rhodes, MIT License) の設計を参考に Rust で再実装
//
// ライセンス: MIT

use std::io::{Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::Mutex;

use crate::chebyshev::{chebyshev_eval3, chebyshev_eval3_with_velocity};
use crate::constants::{J2000_JD, SECS_PER_DAY};

// MARK: - 定数

const RECORD_SIZE: usize = 1024;
const SPK_TYPE_2: i32 = 2;  // Chebyshev 多項式（位置）：惑星・月・太陽
const SPK_TYPE_3: i32 = 3;  // Chebyshev 多項式（位置+速度）：月秤動角、full DE440/441

/// SPK タイプから位置計算に使う成分数を返す（Type 2: 3、Type 3: 6、非対応: None）
fn spk_components(spk_type: i32) -> Option<usize> {
    match spk_type {
        SPK_TYPE_2 => Some(3),
        SPK_TYPE_3 => Some(6),
        _          => None,
    }
}

// MARK: - エラー型

#[derive(Debug)]
pub enum BspError {
    Io(std::io::Error),
    InvalidFormat(String),
    TargetNotFound(i32),
    UnsupportedType(i32),
    OutOfCoverage { input: String, range: String },
}

impl std::fmt::Display for BspError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BspError::Io(e)              => write!(f, "IO error: {e}"),
            BspError::InvalidFormat(msg) => write!(f, "Invalid BSP format: {msg}"),
            BspError::TargetNotFound(id) => write!(f, "Target NAIF ID {id} not found"),
            BspError::UnsupportedType(t) => write!(f, "Unsupported SPK type: {t}"),
            BspError::OutOfCoverage { input, range } =>
                write!(f, "JD out of coverage: {input} (valid: {range})"),
        }
    }
}

impl std::error::Error for BspError {}

impl From<std::io::Error> for BspError {
    fn from(e: std::io::Error) -> Self {
        BspError::Io(e)
    }
}

// MARK: - セグメント

/// パース済みセグメントのメタデータ
#[derive(Debug, Clone)]
pub struct BspSegment {
    pub target:    i32,
    pub center:    i32,
    pub spk_type:  i32,
    /// セグメント開始時刻（J2000 からの秒）
    pub start_sec: f64,
    /// セグメント終了時刻（J2000 からの秒）
    pub end_sec:   f64,
    /// データ開始アドレス（1-indexed double オフセット）
    pub start_idx: usize,
    /// データ終了アドレス（1-indexed double オフセット）
    pub end_idx:   usize,
}

// MARK: - データソース

/// BSP データの保持方式
///
/// - Memory: `from_bytes` でロードした場合（テスト・小ファイル向け）
/// - File: `load` でロードした場合（seek 方式 / 大ファイル向け）
///
/// de440s.bsp（32 MB）はどちらでも動作するが、
/// de441.bsp（約 3 GB）のような大ファイルは File を必ず使うこと。
enum BspData {
    Memory(Vec<u8>),
    File(Mutex<std::fs::File>),
}

impl BspData {
    /// `offset` バイト目から `buf.len()` バイトを読み込む
    ///
    /// # Panics
    /// File モードでシーク/読み取りに失敗した場合（ディスク障害等）
    fn read_at(&self, offset: usize, buf: &mut [u8]) {
        match self {
            BspData::Memory(data) => {
                buf.copy_from_slice(&data[offset..offset + buf.len()]);
            }
            BspData::File(mutex) => {
                let mut file = mutex.lock().expect("BspData mutex not poisoned");
                file.seek(SeekFrom::Start(offset as u64))
                    .expect("BSP seek failed");
                file.read_exact(buf).expect("BSP read failed");
            }
        }
    }
}

// MARK: - BspFile

/// パース済み .bsp ファイルのラッパー
pub struct BspFile {
    /// LOCIFN（ファイル内部名）
    pub name:     String,
    /// セグメント一覧
    pub segments: Vec<BspSegment>,
    data:         BspData,
    is_le:        bool,
}

impl BspFile {
    // MARK: - ロード

    /// .bsp ファイルをパスから seek 方式でロードする（大ファイル対応）
    ///
    /// ファイル全体をメモリに読み込まず、セグメントサマリーのみを解析する。
    /// Chebyshev 係数の読み込みは `compute_chebyshev` 呼び出し時に逐次 seek して行う。
    ///
    /// # メモリ使用量
    /// - サマリー（セグメントメタデータ）: 数 KB 程度
    /// - Chebyshev 係数: 計算時に 1 レコード分（数百バイト）のみ
    pub fn load(path: &Path) -> Result<Self, BspError> {
        let mut file = std::fs::File::open(path)?;

        // ファイルヘッダー（最初の 1024 バイト）を読み込む
        let mut header = vec![0u8; RECORD_SIZE];
        file.read_exact(&mut header)?;
        let (nd, ni, first_sum_rec, is_le, name) = parse_file_record(&header)?;

        // サマリーレコードを seek 方式で読み込む（ファイル全体は読まない）
        let segments = parse_summaries_from_file(&mut file, nd, ni, first_sum_rec, is_le)?;

        Ok(Self { name, segments, data: BspData::File(Mutex::new(file)), is_le })
    }

    /// バイト列から直接初期化する（テスト・小ファイル向け）
    ///
    /// ファイル全体を `Vec<u8>` として保持する。
    /// de440s.bsp（32 MB）以下なら問題ないが、
    /// de441.bsp（約 3 GB）には `load()` を使うこと。
    pub fn from_bytes(data: Vec<u8>) -> Result<Self, BspError> {
        let (nd, ni, first_sum_rec, is_le, name) = parse_file_record(&data)?;
        let segments = parse_summaries(&data, nd, ni, first_sum_rec, is_le)?;
        Ok(Self { name, segments, data: BspData::Memory(data), is_le })
    }

    // MARK: - 公開 API

    /// ICRS 位置ベクトル [x, y, z]（km）を返す（直接セグメントのみ）
    pub fn get_position(&self, target: i32, center: i32, jd_tdb: f64) -> Result<[f64; 3], BspError> {
        let seg = self.find_segment(target, center, jd_tdb)
            .ok_or(BspError::TargetNotFound(target))?;
        let components = spk_components(seg.spk_type)
            .ok_or(BspError::UnsupportedType(seg.spk_type))?;
        Ok(self.compute_chebyshev(seg, jd_tdb, components, false)?.0)
    }

    /// ICRS 位置と速度 [x, y, z]（km, km/day）を返す（直接セグメントのみ）
    pub fn get_position_and_velocity(
        &self,
        target: i32,
        center: i32,
        jd_tdb: f64,
    ) -> Result<([f64; 3], [f64; 3]), BspError> {
        let seg = self.find_segment(target, center, jd_tdb)
            .ok_or(BspError::TargetNotFound(target))?;
        let components = spk_components(seg.spk_type)
            .ok_or(BspError::UnsupportedType(seg.spk_type))?;
        let (pos, vel) = self.compute_chebyshev(seg, jd_tdb, components, true)?;
        Ok((pos, vel.unwrap()))
    }

    /// セグメントチェーンを辿って任意の中心座標で位置を合成する
    ///
    /// DE440s のセグメント構成例:
    ///   SSB(0) → Sun(10)
    ///   SSB(0) → EMB(3) → Earth(399)
    ///   SSB(0) → EMB(3) → Moon(301)
    pub fn compute_position(&self, target: i32, center: i32, jd_tdb: f64) -> Result<[f64; 3], BspError> {
        if target == center {
            return Ok([0.0; 3]);
        }

        if self.find_segment(target, center, jd_tdb).is_some() {
            return self.get_position(target, center, jd_tdb);
        }

        let pt = self.pos_from_ssb(target, jd_tdb)?;
        let pc = if center == 0 {
            [0.0f64; 3]
        } else {
            self.pos_from_ssb(center, jd_tdb)?
        };

        Ok([pt[0] - pc[0], pt[1] - pc[1], pt[2] - pc[2]])
    }

    // MARK: - プライベート：セグメント検索

    fn find_segment(&self, target: i32, center: i32, jd_tdb: f64) -> Option<&BspSegment> {
        let t_sec = (jd_tdb - J2000_JD) * SECS_PER_DAY;
        self.segments.iter().find(|s| {
            s.target == target &&
            s.center == center &&
            t_sec >= s.start_sec &&
            t_sec <= s.end_sec
        })
    }

    /// SSB(0) からの位置を再帰的に合成する
    fn pos_from_ssb(&self, target: i32, jd_tdb: f64) -> Result<[f64; 3], BspError> {
        let t_sec = (jd_tdb - J2000_JD) * SECS_PER_DAY;

        // SSB 直接セグメントがあればそれを使う
        if let Some(seg) = self.find_segment(target, 0, jd_tdb) {
            let components = spk_components(seg.spk_type)
                .ok_or(BspError::UnsupportedType(seg.spk_type))?;
            return Ok(self.compute_chebyshev(seg, jd_tdb, components, false)?.0);
        }

        // 中間天体経由で合成
        for seg in &self.segments {
            if seg.target != target { continue; }
            if t_sec < seg.start_sec || t_sec > seg.end_sec { continue; }
            let components = spk_components(seg.spk_type)
                .ok_or(BspError::UnsupportedType(seg.spk_type))?;
            let from_center = self.compute_chebyshev(seg, jd_tdb, components, false)?.0;
            let center_from_ssb = if seg.center == 0 {
                [0.0f64; 3]
            } else {
                self.pos_from_ssb(seg.center, jd_tdb)?
            };
            return Ok([
                center_from_ssb[0] + from_center[0],
                center_from_ssb[1] + from_center[1],
                center_from_ssb[2] + from_center[2],
            ]);
        }

        Err(BspError::TargetNotFound(target))
    }

    // MARK: - プライベート：Chebyshev 計算（Type 2 / Type 3 共通）

    /// SPK Type 2 / Type 3 セグメントから位置（および速度）を計算する。
    ///
    /// - `components`: 多項式成分数（Type 2 = 3、Type 3 = 6）
    ///
    /// Type 2 と Type 3 のレコード構造は同一であり、成分数のみ異なる。
    /// `ncoeff = (rsize - 2) / components` で自動分岐。
    /// 速度は位置多項式の微分で算出（Type 3 の格納速度係数は不使用）。
    pub(crate) fn compute_chebyshev(
        &self,
        seg: &BspSegment,
        jd_tdb: f64,
        components: usize,
        with_velocity: bool,
    ) -> Result<([f64; 3], Option<[f64; 3]>), BspError> {
        let data_start = (seg.start_idx - 1) * 8;
        let data_end   = seg.end_idx * 8;

        // メタデータ（データ末尾 4 double = 32 bytes）
        let meta_offset = data_end - 32;
        let init_epoch  = self.read_f64(meta_offset);       // 最初のレコード開始時刻（秒）
        let intlen      = self.read_f64(meta_offset + 8);   // 1 レコードあたりの時刻幅（秒）
        let rsize       = self.read_f64(meta_offset + 16).round() as usize;
        let n           = self.read_f64(meta_offset + 24).round() as usize;

        let t_seconds = (jd_tdb - J2000_JD) * SECS_PER_DAY;

        // セグメントのカバー範囲チェック
        let seg_start = init_epoch;
        let seg_end   = init_epoch + (n as f64) * intlen;
        if t_seconds < seg_start || t_seconds > seg_end {
            let jd_start = seg_start / SECS_PER_DAY + J2000_JD;
            let jd_end   = seg_end   / SECS_PER_DAY + J2000_JD;
            return Err(BspError::OutOfCoverage {
                input: format!("JD {jd_tdb:.4}"),
                range: format!("JD {jd_start:.4} – {jd_end:.4}"),
            });
        }

        let idx = {
            let raw = ((t_seconds - init_epoch) / intlen) as isize;
            raw.max(0).min((n as isize) - 1) as usize  // 境界での浮動小数点丸め誤差を吸収
        };

        let rec_offset = data_start + idx * rsize * 8;

        let mid    = self.read_f64(rec_offset);
        let radius = self.read_f64(rec_offset + 8);
        let x      = (t_seconds - mid) / radius;

        let ncoeff = (rsize - 2) / components;  // Type 2: /3、Type 3: /6

        let base   = rec_offset + 16;
        let cx: Vec<f64> = (0..ncoeff).map(|i| self.read_f64(base + i * 8)).collect();
        let cy: Vec<f64> = (0..ncoeff).map(|i| self.read_f64(base + ncoeff * 8 + i * 8)).collect();
        let cz: Vec<f64> = (0..ncoeff).map(|i| self.read_f64(base + ncoeff * 8 * 2 + i * 8)).collect();

        if with_velocity {
            let interval_days = radius * 2.0 / SECS_PER_DAY;
            let (pos, vel) = chebyshev_eval3_with_velocity([&cx, &cy, &cz], x, interval_days);
            Ok((pos, Some(vel)))
        } else {
            Ok((chebyshev_eval3([&cx, &cy, &cz], x), None))
        }
    }

    // MARK: - プライベート：バイナリ読み込みヘルパー

    #[inline]
    fn read_f64(&self, offset: usize) -> f64 {
        let mut buf = [0u8; 8];
        self.data.read_at(offset, &mut buf);
        if self.is_le { f64::from_le_bytes(buf) } else { f64::from_be_bytes(buf) }
    }

    #[allow(dead_code)]
    #[inline]
    fn read_i32(&self, offset: usize) -> i32 {
        let mut buf = [0u8; 4];
        self.data.read_at(offset, &mut buf);
        if self.is_le { i32::from_le_bytes(buf) } else { i32::from_be_bytes(buf) }
    }
}

// MARK: - ファイルレコード解析

fn parse_file_record(
    data: &[u8],
) -> Result<(usize, usize, usize, bool, String), BspError> {

    if data.len() < RECORD_SIZE {
        return Err(BspError::InvalidFormat("ファイルサイズが小さすぎます".into()));
    }

    // LOCIDW（0..8）: "DAF/SPK " で始まるか確認
    let locidw = std::str::from_utf8(&data[0..8]).unwrap_or("");
    if !locidw.starts_with("DAF/SPK") && !locidw.starts_with("DAF/EK") {
        return Err(BspError::InvalidFormat(format!(
            "LOCIDW=\"{}\"", locidw.trim()
        )));
    }

    // LOCFMT（88..96）でエンディアン判定
    let locfmt = std::str::from_utf8(&data[88..96]).unwrap_or("").trim();
    let is_le = locfmt != "BIG-IEEE";

    let read_i32 = |offset: usize| -> i32 {
        let bytes: [u8; 4] = data[offset..offset + 4].try_into().unwrap();
        if is_le { i32::from_le_bytes(bytes) } else { i32::from_be_bytes(bytes) }
    };

    let nd = read_i32(8) as usize;
    let ni = read_i32(12) as usize;
    let first_sum_rec = read_i32(76) as usize;

    let locifn = std::str::from_utf8(&data[16..76])
        .unwrap_or("")
        .trim()
        .to_string();

    Ok((nd, ni, first_sum_rec, is_le, locifn))
}

// MARK: - サマリーレコード解析

fn parse_summaries(
    data: &[u8],
    nd: usize,
    ni: usize,
    first_sum_rec: usize,
    is_le: bool,
) -> Result<Vec<BspSegment>, BspError> {

    let summary_doubles = nd + (ni + 1) / 2;
    let summary_bytes   = summary_doubles * 8;

    let read_f64 = |offset: usize| -> f64 {
        let bytes: [u8; 8] = data[offset..offset + 8].try_into().unwrap();
        if is_le { f64::from_le_bytes(bytes) } else { f64::from_be_bytes(bytes) }
    };
    let read_i32 = |offset: usize| -> i32 {
        let bytes: [u8; 4] = data[offset..offset + 4].try_into().unwrap();
        if is_le { i32::from_le_bytes(bytes) } else { i32::from_be_bytes(bytes) }
    };

    let mut segments = Vec::new();
    let mut rec_num = first_sum_rec;

    while rec_num > 0 {
        let rec_offset  = (rec_num - 1) * RECORD_SIZE;
        let next_rec    = read_f64(rec_offset).round() as usize;
        let n_summaries = read_f64(rec_offset + 16).round() as usize;

        for i in 0..n_summaries {
            let base = rec_offset + 24 + i * summary_bytes;

            let start_sec = read_f64(base);
            let end_sec   = read_f64(base + 8);

            let int_base   = base + nd * 8;
            let target     = read_i32(int_base);
            let center     = read_i32(int_base + 4);
            let spk_type   = read_i32(int_base + 12);
            let first_addr = read_i32(int_base + 16) as usize;
            let last_addr  = read_i32(int_base + 20) as usize;

            segments.push(BspSegment {
                target,
                center,
                spk_type,
                start_sec,
                end_sec,
                start_idx: first_addr,
                end_idx:   last_addr,
            });
        }

        rec_num = next_rec;
    }

    Ok(segments)
}

// MARK: - サマリーレコード解析（seek 版）

/// サマリーレコードをファイルから seek 方式で読み込む
///
/// `load()` 専用。1 レコード（1024 バイト）ずつ seek して読み込むため、
/// ファイル全体をメモリに展開しない。
fn parse_summaries_from_file(
    file: &mut std::fs::File,
    nd: usize,
    ni: usize,
    first_sum_rec: usize,
    is_le: bool,
) -> Result<Vec<BspSegment>, BspError> {
    let summary_doubles = nd + (ni + 1) / 2;
    let summary_bytes   = summary_doubles * 8;

    let read_f64_s = |slice: &[u8], offset: usize| -> f64 {
        let bytes: [u8; 8] = slice[offset..offset + 8].try_into().unwrap();
        if is_le { f64::from_le_bytes(bytes) } else { f64::from_be_bytes(bytes) }
    };
    let read_i32_s = |slice: &[u8], offset: usize| -> i32 {
        let bytes: [u8; 4] = slice[offset..offset + 4].try_into().unwrap();
        if is_le { i32::from_le_bytes(bytes) } else { i32::from_be_bytes(bytes) }
    };

    let mut segments = Vec::new();
    let mut rec_num = first_sum_rec;

    while rec_num > 0 {
        let rec_offset = ((rec_num - 1) * RECORD_SIZE) as u64;
        file.seek(SeekFrom::Start(rec_offset))?;
        let mut record = vec![0u8; RECORD_SIZE];
        file.read_exact(&mut record)?;

        let next_rec    = read_f64_s(&record, 0).round() as usize;
        let n_summaries = read_f64_s(&record, 16).round() as usize;

        for i in 0..n_summaries {
            let base      = 24 + i * summary_bytes;
            let start_sec = read_f64_s(&record, base);
            let end_sec   = read_f64_s(&record, base + 8);
            let int_base  = base + nd * 8;
            let target    = read_i32_s(&record, int_base);
            let center    = read_i32_s(&record, int_base + 4);
            let spk_type  = read_i32_s(&record, int_base + 12);
            let first_addr = read_i32_s(&record, int_base + 16) as usize;
            let last_addr  = read_i32_s(&record, int_base + 20) as usize;

            segments.push(BspSegment {
                target, center, spk_type, start_sec, end_sec,
                start_idx: first_addr, end_idx: last_addr,
            });
        }

        rec_num = next_rec;
    }

    Ok(segments)
}

// MARK: - ユニットテスト

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::J2000_JD;

    // ── 合成 BSP ヘルパー ────────────────────────────────────────────────

    /// 最小限の DAF/SPK Type 2 バイナリを生成する
    ///
    /// Record 1 (0–1023):    ファイルヘッダー
    /// Record 2 (1024–2047): サマリーレコード（Sun/SSB, spkType=2）
    /// Record 3 (2048–2167): Type 2 データ（ncoeff=3、rsize=11、1 レコード）
    ///
    /// 期待位置 (J2000.0): x=100000, y=200000, z=300000 km
    fn make_synthetic_bsp() -> Vec<u8> {
        let mut bytes = vec![0u8; 3 * 1024];

        let write_str = |bytes: &mut Vec<u8>, s: &str, offset: usize, pad: usize| {
            let b = s.as_bytes();
            for i in 0..pad {
                bytes[offset + i] = if i < b.len() { b[i] } else { 0x20 };
            }
        };
        let write_i32 = |bytes: &mut Vec<u8>, v: i32, offset: usize| {
            let b = v.to_le_bytes();
            bytes[offset..offset + 4].copy_from_slice(&b);
        };
        let write_f64 = |bytes: &mut Vec<u8>, v: f64, offset: usize| {
            let b = v.to_le_bytes();
            bytes[offset..offset + 8].copy_from_slice(&b);
        };

        // Record 1: ヘッダー
        write_str(&mut bytes, "DAF/SPK ", 0, 8);
        write_i32(&mut bytes, 2, 8);    // ND
        write_i32(&mut bytes, 6, 12);   // NI
        write_str(&mut bytes, "BspTest", 16, 60);
        write_i32(&mut bytes, 2, 76);   // FWARD
        write_i32(&mut bytes, 2, 80);   // BWARD
        write_str(&mut bytes, "LTL-IEEE", 88, 8);

        // Record 2: サマリー（Type 2、firstAddr=257、lastAddr=271）
        let r2 = 1024;
        write_f64(&mut bytes, 0.0, r2 + 0);   // next rec = none
        write_f64(&mut bytes, 0.0, r2 + 8);
        write_f64(&mut bytes, 1.0, r2 + 16);  // NSUM = 1
        let s = r2 + 24;
        write_f64(&mut bytes, -86400.0, s);
        write_f64(&mut bytes,  86400.0, s + 8);
        write_i32(&mut bytes, 10, s + 16);  // target = Sun
        write_i32(&mut bytes,  0, s + 20);  // center = SSB
        write_i32(&mut bytes,  1, s + 24);  // frame
        write_i32(&mut bytes,  2, s + 28);  // type = 2
        write_i32(&mut bytes, 257, s + 32);
        write_i32(&mut bytes, 271, s + 36);

        // Record 3: Type 2 データ（byte 2048〜）
        let r3 = 2048;
        write_f64(&mut bytes,     0.0, r3 + 0);   // mid
        write_f64(&mut bytes, 86400.0, r3 + 8);   // radius
        write_f64(&mut bytes, 100000.0, r3 + 16); // coeffX[0]
        write_f64(&mut bytes,      0.0, r3 + 24);
        write_f64(&mut bytes,      0.0, r3 + 32);
        write_f64(&mut bytes, 200000.0, r3 + 40); // coeffY[0]
        write_f64(&mut bytes,      0.0, r3 + 48);
        write_f64(&mut bytes,      0.0, r3 + 56);
        write_f64(&mut bytes, 300000.0, r3 + 64); // coeffZ[0]
        write_f64(&mut bytes,      0.0, r3 + 72);
        write_f64(&mut bytes,      0.0, r3 + 80);
        // メタ: dataEnd=271*8=2168 → metaOffset=2136
        write_f64(&mut bytes,  -86400.0, 2136);  // initEpoch
        write_f64(&mut bytes,  172800.0, 2144);  // intlen
        write_f64(&mut bytes,      11.0, 2152);  // rsize
        write_f64(&mut bytes,       1.0, 2160);  // n

        bytes
    }

    /// SPK Type 3 用の最小合成 BSP データを生成する
    ///
    /// rsize=8（2 + 6×1）、ncoeff=1、1 レコード+メタ4=12 doubles → lastAddr=268
    fn make_synthetic_type3_bsp() -> Vec<u8> {
        let mut bytes = vec![0u8; 3 * 1024];

        let write_str = |bytes: &mut Vec<u8>, s: &str, offset: usize, pad: usize| {
            let b = s.as_bytes();
            for i in 0..pad {
                bytes[offset + i] = if i < b.len() { b[i] } else { 0x20 };
            }
        };
        let write_i32 = |bytes: &mut Vec<u8>, v: i32, offset: usize| {
            bytes[offset..offset + 4].copy_from_slice(&v.to_le_bytes());
        };
        let write_f64 = |bytes: &mut Vec<u8>, v: f64, offset: usize| {
            bytes[offset..offset + 8].copy_from_slice(&v.to_le_bytes());
        };

        write_str(&mut bytes, "DAF/SPK ", 0, 8);
        write_i32(&mut bytes, 2, 8);
        write_i32(&mut bytes, 6, 12);
        write_str(&mut bytes, "BspTestType3", 16, 60);
        write_i32(&mut bytes, 2, 76);
        write_i32(&mut bytes, 2, 80);
        write_str(&mut bytes, "LTL-IEEE", 88, 8);

        let r2 = 1024;
        write_f64(&mut bytes, 0.0, r2 + 0);
        write_f64(&mut bytes, 0.0, r2 + 8);
        write_f64(&mut bytes, 1.0, r2 + 16);
        let s = r2 + 24;
        write_f64(&mut bytes, -86400.0, s);
        write_f64(&mut bytes,  86400.0, s + 8);
        write_i32(&mut bytes, 10, s + 16);  // target = Sun
        write_i32(&mut bytes,  0, s + 20);  // center = SSB
        write_i32(&mut bytes,  1, s + 24);
        write_i32(&mut bytes,  3, s + 28);  // type = 3
        write_i32(&mut bytes, 257, s + 32);
        write_i32(&mut bytes, 268, s + 36);

        // Record 3: [mid, radius, Xpos, Ypos, Zpos, Xvel, Yvel, Zvel]
        let r3 = 2048;
        write_f64(&mut bytes,     0.0, r3 + 0);
        write_f64(&mut bytes, 86400.0, r3 + 8);
        write_f64(&mut bytes, 100000.0, r3 + 16); // Xpos
        write_f64(&mut bytes, 200000.0, r3 + 24); // Ypos
        write_f64(&mut bytes, 300000.0, r3 + 32); // Zpos
        // Xvel, Yvel, Zvel = 0
        // メタ: dataEnd=268*8=2144 → metaOffset=2112
        write_f64(&mut bytes, -86400.0, 2112);
        write_f64(&mut bytes, 172800.0, 2120);
        write_f64(&mut bytes,      8.0, 2128);  // rsize = 8
        write_f64(&mut bytes,      1.0, 2136);  // n = 1

        bytes
    }

    // ── Type 2 テスト ─────────────────────────────────────────────────────

    #[test]
    fn test_parse_segments() {
        let bsp = BspFile::from_bytes(make_synthetic_bsp()).unwrap();
        assert_eq!(bsp.segments.len(), 1);
        let seg = &bsp.segments[0];
        assert_eq!(seg.target, 10);
        assert_eq!(seg.center, 0);
        assert_eq!(seg.spk_type, 2);
        assert_eq!(seg.start_idx, 257);
        assert_eq!(seg.end_idx, 271);
    }

    #[test]
    fn test_get_position_j2000() {
        let bsp = BspFile::from_bytes(make_synthetic_bsp()).unwrap();
        let pos = bsp.get_position(10, 0, J2000_JD).unwrap();
        assert!((pos[0] - 100000.0).abs() < 1e-6);
        assert!((pos[1] - 200000.0).abs() < 1e-6);
        assert!((pos[2] - 300000.0).abs() < 1e-6);
    }

    #[test]
    fn test_compute_position_matches_get() {
        let bsp = BspFile::from_bytes(make_synthetic_bsp()).unwrap();
        let p1 = bsp.get_position(10, 0, J2000_JD).unwrap();
        let p2 = bsp.compute_position(10, 0, J2000_JD).unwrap();
        assert!((p1[0] - p2[0]).abs() < 1e-9);
    }

    #[test]
    fn test_same_target_center_returns_zero() {
        let bsp = BspFile::from_bytes(make_synthetic_bsp()).unwrap();
        let pos = bsp.compute_position(10, 10, J2000_JD).unwrap();
        assert_eq!(pos, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn test_target_not_found() {
        let bsp = BspFile::from_bytes(make_synthetic_bsp()).unwrap();
        assert!(matches!(
            bsp.get_position(5, 0, J2000_JD),
            Err(BspError::TargetNotFound(5))
        ));
    }

    #[test]
    fn test_invalid_format() {
        let mut bad = vec![0u8; 1024];
        bad[0] = b'X';
        assert!(matches!(
            BspFile::from_bytes(bad),
            Err(BspError::InvalidFormat(_))
        ));
    }

    #[test]
    fn test_velocity_constant_coeffs_is_zero() {
        let bsp = BspFile::from_bytes(make_synthetic_bsp()).unwrap();
        let (pos, vel) = bsp.get_position_and_velocity(10, 0, J2000_JD).unwrap();
        assert!((pos[0] - 100000.0).abs() < 1e-6);
        assert!(vel[0].abs() < 1e-6);
        assert!(vel[1].abs() < 1e-6);
        assert!(vel[2].abs() < 1e-6);
    }

    // ── out-of-range テスト ──────────────────────────────────────────────

    /// 範囲外 JD はいずれかの BspError を返す。
    /// find_segment が時刻フィルタリングするため TargetNotFound が返るが、
    /// 重要な保証は「サイレントに誤値を返さないこと」。
    #[test]
    fn test_out_of_range_returns_error() {
        let bsp = BspFile::from_bytes(make_synthetic_bsp()).unwrap();
        let jd_outside = J2000_JD + 7.0;  // セグメントは ±1 day のみ
        assert!(bsp.get_position(10, 0, jd_outside).is_err());
    }

    /// compute_chebyshev の内部範囲チェックを直接検証する。
    /// find_segment を経由せず、セグメントを直接渡して範囲外を確認する。
    #[test]
    fn test_out_of_coverage_error_from_compute() {
        let bsp = BspFile::from_bytes(make_synthetic_bsp()).unwrap();
        // セグメントを直接取得（時刻フィルタなし）
        let seg = &bsp.segments[0];
        let jd_outside = J2000_JD + 7.0;
        let result = bsp.compute_chebyshev(seg, jd_outside, 3, false);
        assert!(matches!(result, Err(BspError::OutOfCoverage { .. })));
    }

    // ── Type 3 テスト ─────────────────────────────────────────────────────

    #[test]
    fn test_type3_segment_metadata() {
        let bsp = BspFile::from_bytes(make_synthetic_type3_bsp()).unwrap();
        assert_eq!(bsp.segments[0].spk_type, 3);
    }

    #[test]
    fn test_type3_get_position() {
        let bsp = BspFile::from_bytes(make_synthetic_type3_bsp()).unwrap();
        let pos = bsp.get_position(10, 0, J2000_JD).unwrap();
        assert!((pos[0] - 100000.0).abs() < 1e-6);
        assert!((pos[1] - 200000.0).abs() < 1e-6);
        assert!((pos[2] - 300000.0).abs() < 1e-6);
    }

    #[test]
    fn test_type3_velocity_constant_coeffs_is_zero() {
        let bsp = BspFile::from_bytes(make_synthetic_type3_bsp()).unwrap();
        let (_pos, vel) = bsp.get_position_and_velocity(10, 0, J2000_JD).unwrap();
        assert!(vel[0].abs() < 1e-6);
        assert!(vel[1].abs() < 1e-6);
        assert!(vel[2].abs() < 1e-6);
    }

    #[test]
    fn test_type3_compute_position_matches_get() {
        let bsp = BspFile::from_bytes(make_synthetic_type3_bsp()).unwrap();
        let p1 = bsp.get_position(10, 0, J2000_JD).unwrap();
        let p2 = bsp.compute_position(10, 0, J2000_JD).unwrap();
        assert!((p1[0] - p2[0]).abs() < 1e-9);
    }

    // ── unsupported type (Type 13) テスト ────────────────────────────────

    #[test]
    fn test_unsupported_type13_throws() {
        // Type 3 データの spkType を 13 に改変
        let mut data = make_synthetic_type3_bsp();
        // サマリー中 spkType フィールド: byte 1024+24+28 = 1076
        data[1076] = 13; data[1077] = 0; data[1078] = 0; data[1079] = 0;
        let bsp = BspFile::from_bytes(data).unwrap();
        assert!(matches!(
            bsp.get_position(10, 0, J2000_JD),
            Err(BspError::UnsupportedType(13))
        ));
    }
}
