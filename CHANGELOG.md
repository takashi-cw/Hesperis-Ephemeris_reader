# 更新履歴 (Changelog)

Hesperis Ephemeris Reader（BSP リーダー）としての変更履歴です。
呼び出し側アプリケーション（占星術・天文計算等）固有の変更は対象外です。

---

## 2026-08-03

### 修正（Python 版 / Rust + PyO3）
- `compute_apparent_batch` で光偏差（deflection）が独立して制御できない不具合を修正
  - 従来は `aberration` 引数1個で光偏差・年周光行差をまとめて制御しており、
    光偏差のみ ON・年周光行差のみ OFF（またはその逆）という組み合わせを
    指定できなかった
  - `deflection: bool` 引数を新規追加し、光偏差・年周光行差を
    独立に ON/OFF できるよう変更
- API 変更: `compute_apparent_batch(naif_target, center_naif, jd_tdb_list, use_j2000, aberration, deflection=True)`
  - `deflection` はキーワード引数としてデフォルト値 `True` を持つため、
    既存の5引数呼び出し（`deflection` を渡さない呼び出し）はそのまま動作する
    （PyO3 `#[pyo3(signature = ...)]` によるデフォルト値対応。破壊的変更なし）

## 2026-06-16

### 追加
- `READER_GUIDE.md` を追加
  - 人間・AI 双方向けのドキュメントナビゲーションガイド
  - どのドキュメントが一次情報源か、何を仕様として扱わないかを明記

## 2026-06-14

### 追加
- `SecurityPolicy.md` を追加
  - 脆弱性・秘密情報の露出等を発見した場合の非公開報告手順を明記

## 2026-05-27

### 追加
- 全5言語実装（JavaScript, Python, Swift, Rust, Julia）で SPK Type 3
  （Chebyshev 多項式・位置＋速度、full DE440 / DE441 の月秤動角セグメント）に対応
- 各言語実装に BSP バリデータ・テストスイートを追加

### 変更
- 非対応 SPK タイプ（Type 13 等）に対して明示的なエラーを送出するよう統一
  （範囲外インデックスをクランプしてサイレントに誤った値を返す実装があったため）

### 変更（Python 版 / Rust + PyO3 のみ）
- BSP ファイルの読み込み方式を「ファイル全体をメモリに読み込む」方式から
  「ヘッダーとサマリーのみ読み込み、係数データは都度 seek して読む」方式に変更
  - メモリ使用量の目安: `de441.bsp`（フルカーネル）で 〜3GB → 数 KB に削減
- セグメント時刻の単位ミスマッチを修正
  （DAF サマリーの値は J2000.0 からの秒数であり、ユリウス日ではなかった）

## 2026-05-19

### 追加
- Rust（ネイティブ実装）を追加
- Julia 実装を追加

## 2026-05-16

### 追加
- Swift 実装（パッケージ雛形）を追加

## 2026-05-15

### 追加
- Python 実装（Rust + PyO3 バインディング）を追加
- JavaScript 実装のエントリポイント（`index.js`）を追加

## 2026-04-27

### 追加
- JavaScript 実装に定数モジュール（NAIF コード等）を追加

## 2026-03-28

### 追加
- JavaScript 実装（Chebyshev 多項式補間ロジック）を追加
