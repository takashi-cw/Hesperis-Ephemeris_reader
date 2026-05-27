import Testing
import Foundation
@testable import StellaBspReader

// MARK: - Constants

@Suite("Constants")
struct ConstantsTests {

    @Test func j2000Jd() {
        #expect(J2000_JD == 2_451_545.0)
    }

    @Test func auKm() {
        #expect(AU_KM == 149_597_870.700)
    }

    @Test func speedOfLight() {
        #expect(SPEED_OF_LIGHT_M_S == 299_792_458.0)
    }

    @Test func degToRad() {
        let result = 180.0 * DEG_TO_RAD
        #expect(abs(result - Double.pi) < 1e-15)
    }

    @Test func naifCodes() {
        #expect(NAIF.SUN  == 10)
        #expect(NAIF.MOON == 301)
        #expect(NAIF.EARTH == 399)
        #expect(NAIF.JUPITER_BARYCENTER == 5)
        #expect(NAIF.PLUTO_BARYCENTER   == 9)
    }

    @Test func planetNaifMap() {
        #expect(PLANET_NAIF["Sun"]     == NAIF.SUN)
        #expect(PLANET_NAIF["Jupiter"] == NAIF.JUPITER_BARYCENTER)
        #expect(PLANET_NAIF["Pluto"]   == NAIF.PLUTO_BARYCENTER)
    }
}

// MARK: - Chebyshev

@Suite("Chebyshev")
struct ChebyshevTests {

    @Test func evalT2() {
        let coeffs = [0.0, 0.0, 1.0]
        let x = 0.5
        let expected = 2.0 * x * x - 1.0
        #expect(abs(chebyshevEval(coeffs, x) - expected) < 1e-14)
    }

    @Test func evalConstant() {
        #expect(chebyshevEval([3.0], 0.7) == 3.0)
    }

    @Test func evalEmpty() {
        #expect(chebyshevEval([], 0.5) == 0.0)
    }

    @Test func derivT1() {
        let coeffs = [0.0, 1.0]
        let (pos, dpdx) = chebyshevEvalWithDeriv(coeffs, 0.3)
        #expect(abs(pos  - 0.3) < 1e-14)
        #expect(abs(dpdx - 1.0) < 1e-14)
    }

    @Test func velocityScaling() {
        let coeffs = [0.0, 1.0]
        let (_, vel) = chebyshevEvalWithVelocity(coeffs, 0.0, intervalDays: 2.0)
        #expect(abs(vel - 1.0) < 1e-14)
    }

    @Test func eval3() {
        let cx = [1.0, 0.0]
        let cy = [0.0, 1.0]
        let cz = [0.0, 0.0, 1.0]
        let (px, py, pz) = chebyshevEval3([cx, cy, cz], 0.5)
        #expect(abs(px - 1.0)      < 1e-14)
        #expect(abs(py - 0.5)      < 1e-14)
        #expect(abs(pz - (-0.5))   < 1e-14)
    }

    @Test func normalizeTimeEdges() {
        #expect(abs(normalizeTime(jd: 2451545.5, tStart: 2451545.0, tEnd: 2451546.0) - 0.0)  < 1e-14)
        #expect(abs(normalizeTime(jd: 2451545.0, tStart: 2451545.0, tEnd: 2451546.0) - (-1.0)) < 1e-14)
        #expect(abs(normalizeTime(jd: 2451546.0, tStart: 2451545.0, tEnd: 2451546.0) - 1.0)  < 1e-14)
    }
}

// MARK: - BspValidator

@Suite("BspValidator")
struct BspValidatorTests {

    private func makeSegments() -> [BspSegment] {
        let s_per_day = 86400.0
        let startSec = -54787.5 * s_per_day
        let endSec   =  54787.5 * s_per_day
        return [
            BspSegment(target: 10, center: 0, spkType: 2,
                       startSec: startSec, endSec: endSec,
                       startIdx: 0, endIdx: 0)
        ]
    }

    @Test func coverageJd() {
        let segs = makeSegments()
        let cov  = getCoverageJd(segments: segs, naifTarget: 10)
        let startYear = 2000.0 + (cov.startJd - J2000_JD) / 365.25
        let endYear   = 2000.0 + (cov.endJd   - J2000_JD) / 365.25
        #expect(startYear < 1851.0)
        #expect(endYear   > 2149.0)
    }

    @Test func formatMessage() {
        let segs = makeSegments()
        let cov  = getCoverageJd(segments: segs)
        let msg  = formatCoverageMessage(cov)
        #expect(msg.contains("AD1"))
        #expect(msg.contains("〜"))
        #expect(msg.contains("AD2"))
    }

    @Test func assertInCoveragePass() throws {
        let segs = makeSegments()
        try assertInCoverage(jdTdb: J2000_JD, segments: segs)
    }

    @Test func assertInCoverageThrows() {
        let segs = makeSegments()
        let jd1800 = J2000_JD - 200.0 * 365.25
        #expect(throws: BspError.self) {
            try assertInCoverage(jdTdb: jd1800, segments: segs)
        }
    }

    @Test func emptySegmentsReturnsJ2000() {
        let cov = getCoverageJd(segments: [])
        #expect(cov.startJd == J2000_JD)
        #expect(cov.endJd   == J2000_JD)
    }
}

// MARK: - BspReader

private func makeSyntheticBspData() -> Data {
    var bytes = [UInt8](repeating: 0, count: 3 * 1024)

    func writeStr(_ s: String, at offset: Int, padTo length: Int) {
        let arr = Array(s.utf8)
        for i in 0..<min(arr.count, length) { bytes[offset + i] = arr[i] }
        for i in arr.count..<length { bytes[offset + i] = 0x20 }
    }
    func writeI32LE(_ v: Int32, at offset: Int) {
        var val = v.littleEndian
        withUnsafeBytes(of: &val) { src in for i in 0..<4 { bytes[offset + i] = src[i] } }
    }
    func writeDblLE(_ v: Double, at offset: Int) {
        var val = v.bitPattern.littleEndian
        withUnsafeBytes(of: &val) { src in for i in 0..<8 { bytes[offset + i] = src[i] } }
    }

    writeStr("DAF/SPK ", at: 0,  padTo: 8)
    writeI32LE(2, at: 8); writeI32LE(6, at: 12)
    writeStr("StellaBspReader TEST", at: 16, padTo: 60)
    writeI32LE(2, at: 76); writeI32LE(2, at: 80)
    writeStr("LTL-IEEE", at: 88, padTo: 8)

    let r2 = 1024
    writeDblLE(0.0, at: r2 + 0); writeDblLE(0.0, at: r2 + 8); writeDblLE(1.0, at: r2 + 16)
    let s = r2 + 24
    writeDblLE(-86400.0, at: s); writeDblLE(86400.0, at: s + 8)
    writeI32LE(10, at: s + 16); writeI32LE(0, at: s + 20)
    writeI32LE(1,  at: s + 24); writeI32LE(2, at: s + 28)
    writeI32LE(257, at: s + 32); writeI32LE(271, at: s + 36)

    let r3 = 2048
    writeDblLE(0.0, at: r3 + 0); writeDblLE(86400.0, at: r3 + 8)
    writeDblLE(100000.0, at: r3 + 16); writeDblLE(0.0, at: r3 + 24); writeDblLE(0.0, at: r3 + 32)
    writeDblLE(200000.0, at: r3 + 40); writeDblLE(0.0, at: r3 + 48); writeDblLE(0.0, at: r3 + 56)
    writeDblLE(300000.0, at: r3 + 64); writeDblLE(0.0, at: r3 + 72); writeDblLE(0.0, at: r3 + 80)
    writeDblLE(-86400.0, at: 2136); writeDblLE(172800.0, at: 2144)
    writeDblLE(11.0, at: 2152);     writeDblLE(1.0, at: 2160)

    return Data(bytes)
}

/// SPK Type 3 用の最小合成 BSP データを生成する
///
/// 構造:
///   Record 1 (bytes    0–1023): ファイルヘッダー
///   Record 2 (bytes 1024–2047): サマリーレコード（Sun/SSB セグメント、spkType=3）
///   Record 3 (bytes 2048–2143): Type 3 データ（ncoeff=1、成分数=6、1 レコード）
private func makeSyntheticType3BspData() -> Data {
    var bytes = [UInt8](repeating: 0, count: 3 * 1024)

    func writeStr(_ s: String, at offset: Int, padTo length: Int) {
        let arr = Array(s.utf8)
        for i in 0..<min(arr.count, length) { bytes[offset + i] = arr[i] }
        for i in arr.count..<length { bytes[offset + i] = 0x20 }
    }
    func writeI32LE(_ v: Int32, at offset: Int) {
        var val = v.littleEndian
        withUnsafeBytes(of: &val) { src in for i in 0..<4 { bytes[offset + i] = src[i] } }
    }
    func writeDblLE(_ v: Double, at offset: Int) {
        var val = v.bitPattern.littleEndian
        withUnsafeBytes(of: &val) { src in for i in 0..<8 { bytes[offset + i] = src[i] } }
    }

    writeStr("DAF/SPK ", at: 0,  padTo: 8)
    writeI32LE(2, at: 8); writeI32LE(6, at: 12)
    writeStr("StellaBspReader TEST Type3", at: 16, padTo: 60)
    writeI32LE(2, at: 76); writeI32LE(2, at: 80)
    writeStr("LTL-IEEE", at: 88, padTo: 8)

    // Record 2: サマリー（Type 3）
    // ncoeff=1, rsize=8 (2+6×1), 1レコード+メタ4=12 doubles → lastAddr=268
    let r2 = 1024
    writeDblLE(0.0, at: r2 + 0); writeDblLE(0.0, at: r2 + 8); writeDblLE(1.0, at: r2 + 16)
    let s = r2 + 24
    writeDblLE(-86400.0, at: s); writeDblLE(86400.0, at: s + 8)
    writeI32LE(10, at: s + 16); writeI32LE(0, at: s + 20)
    writeI32LE(1,  at: s + 24); writeI32LE(3, at: s + 28)  // type = 3
    writeI32LE(257, at: s + 32); writeI32LE(268, at: s + 36)

    // Record 3: Type 3 データ [mid, radius, Xpos, Ypos, Zpos, Xvel, Yvel, Zvel]
    let r3 = 2048
    writeDblLE(0.0, at: r3 + 0); writeDblLE(86400.0, at: r3 + 8)
    writeDblLE(100000.0, at: r3 + 16)  // coeffX_pos
    writeDblLE(200000.0, at: r3 + 24)  // coeffY_pos
    writeDblLE(300000.0, at: r3 + 32)  // coeffZ_pos
    writeDblLE(0.0, at: r3 + 40); writeDblLE(0.0, at: r3 + 48); writeDblLE(0.0, at: r3 + 56)
    // メタデータ: metaOffset = 268*8 - 32 = 2112
    writeDblLE(-86400.0, at: 2112); writeDblLE(172800.0, at: 2120)
    writeDblLE(8.0, at: 2128);      writeDblLE(1.0, at: 2136)

    return Data(bytes)
}

@Suite("BspReader")
struct BspReaderTests {

    @Test func parseFileRecord() throws {
        let bsp = try BspFile(data: makeSyntheticBspData())
        #expect(bsp.name.contains("StellaBspReader TEST"))
        #expect(bsp.segments.count == 1)
    }

    @Test func segmentMetadata() throws {
        let bsp = try BspFile(data: makeSyntheticBspData())
        let seg = try #require(bsp.segments.first)
        #expect(seg.target == 10)
        #expect(seg.center == 0)
        #expect(seg.spkType == 2)
        #expect(seg.startSec == -86400.0)
        #expect(seg.endSec   ==  86400.0)
        #expect(seg.startIdx == 257)
        #expect(seg.endIdx   == 271)
    }

    @Test func getPositionAtJ2000() throws {
        let bsp = try BspFile(data: makeSyntheticBspData())
        let (x, y, z) = try bsp.getPosition(target: 10, center: 0, jdTdb: J2000_JD)
        #expect(abs(x - 100000.0) < 1e-6)
        #expect(abs(y - 200000.0) < 1e-6)
        #expect(abs(z - 300000.0) < 1e-6)
    }

    @Test func computePositionSameAsGet() throws {
        let bsp = try BspFile(data: makeSyntheticBspData())
        let pos1 = try bsp.getPosition(target: 10, center: 0, jdTdb: J2000_JD)
        let pos2 = try bsp.computePosition(target: 10, center: 0, jdTdb: J2000_JD)
        #expect(abs(pos1.0 - pos2.0) < 1e-9)
        #expect(abs(pos1.1 - pos2.1) < 1e-9)
        #expect(abs(pos1.2 - pos2.2) < 1e-9)
    }

    @Test func sameTargetCenterReturnsZero() throws {
        let bsp = try BspFile(data: makeSyntheticBspData())
        let (x, y, z) = try bsp.computePosition(target: 10, center: 10, jdTdb: J2000_JD)
        #expect(x == 0.0 && y == 0.0 && z == 0.0)
    }

    @Test func targetNotFoundThrows() throws {
        let bsp = try BspFile(data: makeSyntheticBspData())
        #expect(throws: BspError.self) {
            try bsp.getPosition(target: 5, center: 0, jdTdb: J2000_JD)
        }
    }

    @Test func invalidFormatThrows() {
        var bad = [UInt8](repeating: 0, count: 1024)
        bad[0] = 0x58
        #expect(throws: BspError.self) {
            try BspFile(data: Data(bad))
        }
    }

    @Test func positionAndVelocity() throws {
        let bsp = try BspFile(data: makeSyntheticBspData())
        let (pos, vel) = try bsp.getPositionAndVelocity(target: 10, center: 0, jdTdb: J2000_JD)
        #expect(abs(pos.0 - 100000.0) < 1e-6)
        #expect(abs(vel.0) < 1e-6)
        #expect(abs(vel.1) < 1e-6)
        #expect(abs(vel.2) < 1e-6)
    }

    // MARK: out-of-range

    @Test func outOfRangeThrows() throws {
        let bsp = try BspFile(data: makeSyntheticBspData())
        // セグメントは J2000 ± 1 day。1 week 後は範囲外
        let jdOutside = J2000_JD + 7.0
        #expect(throws: BspError.self) {
            try bsp.getPosition(target: 10, center: 0, jdTdb: jdOutside)
        }
    }

    // MARK: Type 3

    @Test func type3SegmentMetadata() throws {
        let bsp = try BspFile(data: makeSyntheticType3BspData())
        let seg = try #require(bsp.segments.first)
        #expect(seg.spkType == 3)
        #expect(seg.target == 10)
        #expect(seg.center == 0)
    }

    @Test func type3GetPositionAtJ2000() throws {
        let bsp = try BspFile(data: makeSyntheticType3BspData())
        let (x, y, z) = try bsp.getPosition(target: 10, center: 0, jdTdb: J2000_JD)
        #expect(abs(x - 100000.0) < 1e-6)
        #expect(abs(y - 200000.0) < 1e-6)
        #expect(abs(z - 300000.0) < 1e-6)
    }

    @Test func type3PositionAndVelocityConstantCoeffs() throws {
        let bsp = try BspFile(data: makeSyntheticType3BspData())
        let (pos, vel) = try bsp.getPositionAndVelocity(target: 10, center: 0, jdTdb: J2000_JD)
        #expect(abs(pos.0 - 100000.0) < 1e-6)
        // 定数係数 → 速度 = 0
        #expect(abs(vel.0) < 1e-6)
        #expect(abs(vel.1) < 1e-6)
        #expect(abs(vel.2) < 1e-6)
    }

    @Test func type3ComputePositionMatchesGet() throws {
        let bsp = try BspFile(data: makeSyntheticType3BspData())
        let pos1 = try bsp.getPosition(target: 10, center: 0, jdTdb: J2000_JD)
        let pos2 = try bsp.computePosition(target: 10, center: 0, jdTdb: J2000_JD)
        #expect(abs(pos1.0 - pos2.0) < 1e-9)
        #expect(abs(pos1.1 - pos2.1) < 1e-9)
        #expect(abs(pos1.2 - pos2.2) < 1e-9)
    }

    // MARK: unsupportedType (Type 13)

    @Test func unsupportedType13Throws() throws {
        // Type 3 データを改変して Type 13 にする
        var data = makeSyntheticType3BspData()
        // サマリー中の spkType フィールド: byte 1024+24+28 = 1076
        data[1076] = 13; data[1077] = 0; data[1078] = 0; data[1079] = 0
        let bsp = try BspFile(data: data)
        #expect(throws: BspError.self) {
            try bsp.getPosition(target: 10, center: 0, jdTdb: J2000_JD)
        }
    }
}
