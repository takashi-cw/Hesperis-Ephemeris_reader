// BspReader.swift — JPL .bsp バイナリ読み込み
//
// NASA NAIF DAF/SPK フォーマットを Foundation.Data で解析し、
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
//   - jplephem (Brandon Rhodes, MIT License) の設計を参考に Swift で再実装
//
// 依存: Constants.swift / Chebyshev.swift / BspValidator.swift
//
// ライセンス: MIT

import Foundation

// MARK: - 定数（ファイルローカル）

private let recordSize  = 1024
private let sPerDay     = 86400.0
private let spkType2: Int32 = 2  // Chebyshev 多項式（位置）：惑星・月・太陽
private let spkType3: Int32 = 3  // Chebyshev 多項式（位置+速度）：月秤動角、full DE440/441

/// SPK タイプから位置計算に使う成分数を返す（Type 2: 3、Type 3: 6）
private func spkComponents(_ type: Int32) -> Int? {
    switch type {
    case spkType2: return 3
    case spkType3: return 6
    default: return nil
    }
}

// MARK: - BspFile

/// パース済み .bsp ファイルのラッパー
public struct BspFile: Sendable {

    /// LOCIFN（ファイル内部名）
    public let name: String

    /// セグメント一覧
    public let segments: [BspSegment]

    private let data: Data
    private let isLittleEndian: Bool

    // MARK: - 初期化

    /// Data から直接初期化する（テスト・インメモリ用途向け）
    public init(data: Data) throws {
        self.data = data
        let (nd, ni, firstSumRec, isLE, locifn) = try Self.parseFileRecord(data)
        self.isLittleEndian = isLE
        self.name = locifn
        self.segments = try Self.parseSummaries(data, nd: nd, ni: ni,
                                                firstSumRec: firstSumRec,
                                                isLE: isLE)
    }

    /// .bsp ファイルを URL から読み込む
    ///
    /// - Parameter url: ファイルURL（`URL(fileURLWithPath:)` で生成すること）
    /// - Returns: パース済み BspFile
    /// - Throws: `BspError.invalidFormat` または `Data` の読み込みエラー
    public static func load(url: URL) throws -> BspFile {
        let data = try Data(contentsOf: url, options: .mappedIfSafe)
        return try BspFile(data: data)
    }

    // MARK: - 公開 API

    /// ICRS 位置ベクトル（km）を返す（直接セグメントのみ）
    ///
    /// - Parameters:
    ///   - target: NAIF ターゲットコード
    ///   - center: NAIF センターコード
    ///   - jdTdb: ユリウス日（TDB）
    /// - Returns: (x, y, z) km
    public func getPosition(target: Int32, center: Int32, jdTdb: Double) throws -> (Double, Double, Double) {
        guard let seg = findSegment(target: target, center: center, jdTdb: jdTdb) else {
            throw BspError.targetNotFound(target)
        }
        guard let components = spkComponents(seg.spkType) else {
            throw BspError.unsupportedType(seg.spkType)
        }
        return try computeChebyshev(segment: seg, jdTdb: jdTdb, components: components, withVelocity: false).position
    }

    /// ICRS 位置と速度ベクトル（km, km/day）を返す（直接セグメントのみ）
    public func getPositionAndVelocity(
        target: Int32, center: Int32, jdTdb: Double
    ) throws -> (position: (Double, Double, Double), velocity: (Double, Double, Double)) {
        guard let seg = findSegment(target: target, center: center, jdTdb: jdTdb) else {
            throw BspError.targetNotFound(target)
        }
        guard let components = spkComponents(seg.spkType) else {
            throw BspError.unsupportedType(seg.spkType)
        }
        let result = try computeChebyshev(segment: seg, jdTdb: jdTdb, components: components, withVelocity: true)
        return (result.position, result.velocity!)
    }

    /// セグメントチェーンを辿って任意の中心座標で位置を合成する
    ///
    /// DE440s のセグメント構成例:
    ///   SSB(0) → Sun(10)
    ///   SSB(0) → EMB(3) → Earth(399)
    ///   SSB(0) → EMB(3) → Moon(301)
    ///
    /// - Parameters:
    ///   - target: NAIF ターゲットコード
    ///   - center: NAIF センターコード（通常 0 = SSB）
    ///   - jdTdb: ユリウス日（TDB）
    /// - Returns: (x, y, z) km
    public func computePosition(target: Int32, center: Int32, jdTdb: Double) throws -> (Double, Double, Double) {
        if target == center { return (0.0, 0.0, 0.0) }

        if findSegment(target: target, center: center, jdTdb: jdTdb) != nil {
            return try getPosition(target: target, center: center, jdTdb: jdTdb)
        }

        let ssb: Int32 = 0
        let pt = try posFromSsb(target: target, jdTdb: jdTdb)
        let pc: (Double, Double, Double) = (center == ssb)
            ? (0.0, 0.0, 0.0)
            : try posFromSsb(target: center, jdTdb: jdTdb)

        return (pt.0 - pc.0, pt.1 - pc.1, pt.2 - pc.2)
    }

    // MARK: - プライベート：セグメント検索

    private func findSegment(target: Int32, center: Int32, jdTdb: Double) -> BspSegment? {
        let tSec = (jdTdb - J2000_JD) * sPerDay
        return segments.first {
            $0.target == target &&
            $0.center == center &&
            tSec >= $0.startSec &&
            tSec <= $0.endSec
        }
    }

    /// SSB(0) からの位置を再帰的に合成する
    private func posFromSsb(target: Int32, jdTdb: Double) throws -> (Double, Double, Double) {
        let ssb: Int32 = 0
        let tSec = (jdTdb - J2000_JD) * sPerDay

        if let seg = findSegment(target: target, center: ssb, jdTdb: jdTdb) {
            guard let components = spkComponents(seg.spkType) else { throw BspError.unsupportedType(seg.spkType) }
            return try computeChebyshev(segment: seg, jdTdb: jdTdb, components: components, withVelocity: false).position
        }

        // 中間天体経由で合成
        for seg in segments where seg.target == target && tSec >= seg.startSec && tSec <= seg.endSec {
            guard let components = spkComponents(seg.spkType) else { throw BspError.unsupportedType(seg.spkType) }
            let fromCenter    = try computeChebyshev(segment: seg, jdTdb: jdTdb, components: components, withVelocity: false).position
            let centerFromSsb = seg.center == ssb ? (0.0, 0.0, 0.0) : try posFromSsb(target: seg.center, jdTdb: jdTdb)
            return (
                centerFromSsb.0 + fromCenter.0,
                centerFromSsb.1 + fromCenter.1,
                centerFromSsb.2 + fromCenter.2
            )
        }

        throw BspError.targetNotFound(target)
    }

    // MARK: - プライベート：Chebyshev 計算（Type 2 / Type 3 共通）

    /// SPK Type 2 / Type 3 セグメントから位置（および速度）を計算する。
    ///
    /// - Parameters:
    ///   - segment: 対象セグメント
    ///   - jdTdb: ユリウス日（TDB）
    ///   - components: 多項式成分数（Type 2 = 3、Type 3 = 6）
    ///   - withVelocity: true のとき速度も返す（位置係数を微分して算出）
    ///
    /// Type 2 と Type 3 のレコード構造は同一であり、成分数のみ異なる。
    /// Type 3 は [mid, radius, X_pos×n, Y_pos×n, Z_pos×n, X_vel×n, Y_vel×n, Z_vel×n] の
    /// 順に係数を持つが、本実装では位置 3 成分のみ読み取り、速度は位置多項式の微分で求める。
    private func computeChebyshev(
        segment: BspSegment,
        jdTdb: Double,
        components: Int,
        withVelocity: Bool
    ) throws -> (position: (Double, Double, Double), velocity: (Double, Double, Double)?) {

        let dataStart = (segment.startIdx - 1) * 8
        let dataEnd   = segment.endIdx * 8

        // メタデータ（データ末尾 4 double = 32 bytes）
        let metaOffset = dataEnd - 32
        let initEpoch  = readDouble(at: metaOffset)       // 最初のレコード開始時刻（秒）
        let intlen     = readDouble(at: metaOffset + 8)   // 1 レコードあたりの時刻幅（秒）
        let rsize      = Int(readDouble(at: metaOffset + 16).rounded()) // 1 レコードの Double 数
        let n          = Int(readDouble(at: metaOffset + 24).rounded()) // レコード数

        let tSeconds = (jdTdb - J2000_JD) * sPerDay

        // セグメントのカバー範囲チェック
        let segStart = initEpoch
        let segEnd   = initEpoch + Double(n) * intlen
        guard tSeconds >= segStart && tSeconds <= segEnd else {
            let jdStart = segStart / sPerDay + J2000_JD
            let jdEnd   = segEnd   / sPerDay + J2000_JD
            throw BspError.outOfCoverage(
                input: String(format: "JD %.4f", jdTdb),
                range: String(format: "JD %.4f – %.4f", jdStart, jdEnd)
            )
        }

        var idx = Int((tSeconds - initEpoch) / intlen)
        idx = max(0, min(idx, n - 1))   // 境界での浮動小数点丸め誤差を吸収

        let recOffset = dataStart + idx * rsize * 8

        let mid    = readDouble(at: recOffset)
        let radius = readDouble(at: recOffset + 8)

        let x      = (tSeconds - mid) / radius
        let ncoeff = (rsize - 2) / components  // Type 2: /3、Type 3: /6

        let coeffX = readCoeffs(at: recOffset + 16,                  count: ncoeff)
        let coeffY = readCoeffs(at: recOffset + 16 + ncoeff * 8,     count: ncoeff)
        let coeffZ = readCoeffs(at: recOffset + 16 + ncoeff * 8 * 2, count: ncoeff)

        if withVelocity {
            let intervalDays = radius * 2.0 / sPerDay
            let (pos, vel) = chebyshevEval3WithVelocity([coeffX, coeffY, coeffZ], x, intervalDays: intervalDays)
            return (pos, vel)
        } else {
            return (chebyshevEval3([coeffX, coeffY, coeffZ], x), nil)
        }
    }

    // MARK: - プライベート：バイナリ読み込みヘルパー

    private func readDouble(at offset: Int) -> Double {
        let raw: Double = data.withUnsafeBytes { ptr in
            ptr.loadUnaligned(fromByteOffset: offset, as: Double.self)
        }
        if isLittleEndian { return raw }
        return Double(bitPattern: raw.bitPattern.byteSwapped)
    }

    private func readInt32(at offset: Int) -> Int32 {
        let raw: Int32 = data.withUnsafeBytes { ptr in
            ptr.loadUnaligned(fromByteOffset: offset, as: Int32.self)
        }
        return isLittleEndian ? raw : raw.byteSwapped
    }

    private func readCoeffs(at offset: Int, count: Int) -> [Double] {
        (0..<count).map { readDouble(at: offset + $0 * 8) }
    }

    private func readString(at offset: Int, length: Int) -> String {
        let bytes = data[offset ..< offset + length]
        return String(bytes: bytes, encoding: .ascii) ?? ""
    }

    // MARK: - プライベート：ファイルレコード解析

    private static func parseFileRecord(
        _ data: Data
    ) throws -> (nd: Int, ni: Int, firstSumRec: Int, isLE: Bool, locifn: String) {

        guard data.count >= recordSize else {
            throw BspError.invalidFormat("ファイルサイズが小さすぎます")
        }

        let locidwBytes = data[0..<8]
        let locidw = String(bytes: locidwBytes, encoding: .ascii) ?? ""
        guard locidw.hasPrefix("DAF/SPK") || locidw.hasPrefix("DAF/EK") else {
            throw BspError.invalidFormat("LOCIDW=\"\(locidw.trimmingCharacters(in: .whitespaces))\"")
        }

        // ND, NI はリトルエンディアン確定前なので両方を試す
        // ファイルのエンディアンは LOCFMT フィールドで判定する
        let locfmtBytes = data[88..<96]
        let locfmt = (String(bytes: locfmtBytes, encoding: .ascii) ?? "").trimmingCharacters(in: .whitespaces)
        let isLE = (locfmt != "BIG-IEEE")

        func int32(at offset: Int) -> Int32 {
            let raw: Int32 = data.withUnsafeBytes { ptr in
                ptr.loadUnaligned(fromByteOffset: offset, as: Int32.self)
            }
            return isLE ? raw : raw.byteSwapped
        }

        let nd = Int(int32(at: 8))
        let ni = Int(int32(at: 12))
        let firstSumRec = Int(int32(at: 76))

        let locifnBytes = data[16..<76]
        let locifn = (String(bytes: locifnBytes, encoding: .ascii) ?? "").trimmingCharacters(in: .whitespaces)

        return (nd, ni, firstSumRec, isLE, locifn)
    }

    // MARK: - プライベート：サマリーレコード解析

    private static func parseSummaries(
        _ data: Data,
        nd: Int, ni: Int,
        firstSumRec: Int,
        isLE: Bool
    ) throws -> [BspSegment] {

        let summaryDoubles = nd + (ni + 1) / 2  // ceil(ni/2)
        let summaryBytes   = summaryDoubles * 8

        func dbl(at offset: Int) -> Double {
            let raw: Double = data.withUnsafeBytes { ptr in
                ptr.loadUnaligned(fromByteOffset: offset, as: Double.self)
            }
            if isLE { return raw }
            return Double(bitPattern: raw.bitPattern.byteSwapped)
        }
        func i32(at offset: Int) -> Int32 {
            let raw: Int32 = data.withUnsafeBytes { ptr in
                ptr.loadUnaligned(fromByteOffset: offset, as: Int32.self)
            }
            return isLE ? raw : raw.byteSwapped
        }

        var segments: [BspSegment] = []
        var recNum = firstSumRec

        while recNum > 0 {
            let recOffset  = (recNum - 1) * recordSize
            let nextRec    = Int(dbl(at: recOffset).rounded())
            let nSummaries = Int(dbl(at: recOffset + 16).rounded())

            for i in 0..<nSummaries {
                let base = recOffset + 24 + i * summaryBytes

                let startSec = dbl(at: base)
                let endSec   = dbl(at: base + 8)

                let intBase  = base + nd * 8
                let target   = i32(at: intBase + 0)
                let center   = i32(at: intBase + 4)
                let spkType  = i32(at: intBase + 12)
                let firstAddr = Int(i32(at: intBase + 16))
                let lastAddr  = Int(i32(at: intBase + 20))

                segments.append(BspSegment(
                    target: target, center: center, spkType: spkType,
                    startSec: startSec, endSec: endSec,
                    startIdx: firstAddr, endIdx: lastAddr
                ))
            }

            recNum = nextRec
        }

        return segments
    }
}
