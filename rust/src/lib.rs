// lib.rs — stella-bsp-reader
//
// NASA JPL BSP/SPK バイナリ天体暦ファイルの Pure Rust リーダー
//
// 対応フォーマット: DAF/SPK Type 2（Chebyshev 多項式：位置・速度）
// 外部依存: なし（std のみ）
//
// ライセンス: MIT

pub mod bsp;
pub mod chebyshev;
pub mod constants;

// よく使う型・関数を re-export
pub use bsp::{BspError, BspFile, BspSegment};
pub use constants::naif;
