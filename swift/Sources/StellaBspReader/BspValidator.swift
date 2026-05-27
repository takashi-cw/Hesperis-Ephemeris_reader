// BspValidator.swift — BSP 天体暦のカバー範囲検証
//
// BSP ファイルのセグメント時刻範囲（J2000.0 からの秒数）を JD に変換し、
// 入力日時が天体暦のカバー範囲内かを検証する。
//
// 依存: Constants.swift（J2000_JD）
//
// ライセンス: MIT

import Foundation

// MARK: - カバー範囲

/// BSP ファイルのカバー範囲を表す値型
public struct BspCoverage: Sendable {
    public let startJd: Double
    public let endJd: Double
}

/// BSP ファイルのカバー範囲を JD で返す
///
/// - Parameters:
///   - segments: BspFile のセグメント配列
///   - naifTarget: 代表天体の NAIF コード（デフォルト: 太陽 = 10）
/// - Returns: カバー範囲（JD）
public func getCoverageJd(segments: [BspSegment], naifTarget: Int32 = NAIF.SUN) -> BspCoverage {
    let s_per_day = 86400.0

    let filtered = segments.filter { $0.target == naifTarget }
    let src = filtered.isEmpty ? segments : filtered

    guard !src.isEmpty else {
        return BspCoverage(startJd: J2000_JD, endJd: J2000_JD)
    }

    let minSec = src.map(\.startSec).min()!
    let maxSec = src.map(\.endSec).max()!

    return BspCoverage(
        startJd: J2000_JD + minSec / s_per_day,
        endJd:   J2000_JD + maxSec / s_per_day
    )
}

/// カバー範囲を人間が読める文字列で返す
///
/// - Parameter coverage: getCoverageJd の返値
/// - Returns: 例 "AD1850〜AD2150"（de440s.bsp の場合）
public func formatCoverageMessage(_ coverage: BspCoverage) -> String {
    func label(_ jd: Double) -> String {
        let year = 2000.0 + (jd - J2000_JD) / 365.25
        if year < 0 {
            return "BC\(Int(abs(ceil(year))))"
        } else {
            return "AD\(Int(year))"
        }
    }
    return "\(label(coverage.startJd))〜\(label(coverage.endJd))"
}

/// JD が BSP のカバー範囲内かを検証する
///
/// - Parameters:
///   - jdTdb: 検証する JD（TDB）
///   - segments: BspFile のセグメント配列
/// - Throws: `BspError.outOfCoverage` が範囲外の場合
public func assertInCoverage(jdTdb: Double, segments: [BspSegment]) throws {
    let coverage = getCoverageJd(segments: segments)
    guard jdTdb >= coverage.startJd && jdTdb <= coverage.endJd else {
        let range = formatCoverageMessage(coverage)
        let inputYear = 2000.0 + (jdTdb - J2000_JD) / 365.25
        let inputLabel = inputYear < 0
            ? "BC\(Int(abs(ceil(inputYear))))"
            : "AD\(Int(inputYear))"
        throw BspError.outOfCoverage(input: inputLabel, range: range)
    }
}

// MARK: - エラー型

public enum BspError: Error, Sendable {
    /// ファイルが DAF/SPK フォーマットでない
    case invalidFormat(String)
    /// 対応していない SPK タイプ
    case unsupportedType(Int32)
    /// 指定天体のセグメントが存在しない
    case targetNotFound(Int32)
    /// JD が天体暦の対象期間外
    case outOfCoverage(input: String, range: String)
}

extension BspError: LocalizedError {
    public var errorDescription: String? {
        switch self {
        case .invalidFormat(let msg):
            return "非SPKファイルです: \(msg)"
        case .unsupportedType(let t):
            return "未対応のSPKタイプです: \(t)"
        case .targetNotFound(let code):
            return "指定天体（NAIF \(code)）のセグメントが見つかりません"
        case .outOfCoverage(let input, let range):
            return "天体暦の範囲外です（\(input)）。カバー範囲: \(range)"
        }
    }
}

// MARK: - BspSegment（BspReader 向けの共有型）
// BspReader.swift で完全実装する。ここでは Validator が参照できる最小定義のみ。

/// SPK セグメントのメタデータ
public struct BspSegment: Sendable {
    /// ターゲット天体の NAIF コード
    public let target: Int32
    /// 観測者（中心天体）の NAIF コード
    public let center: Int32
    /// SPK タイプ（Type 2 = Chebyshev 位置のみ）
    public let spkType: Int32
    /// セグメント開始時刻（J2000.0 からの秒数）
    public let startSec: Double
    /// セグメント終了時刻（J2000.0 からの秒数）
    public let endSec: Double
    /// データ配列の先頭インデックス（Double 単位）
    public let startIdx: Int
    /// データ配列の末尾インデックス（Double 単位）
    public let endIdx: Int

    public init(
        target: Int32, center: Int32, spkType: Int32,
        startSec: Double, endSec: Double,
        startIdx: Int, endIdx: Int
    ) {
        self.target   = target
        self.center   = center
        self.spkType  = spkType
        self.startSec = startSec
        self.endSec   = endSec
        self.startIdx = startIdx
        self.endIdx   = endIdx
    }
}
