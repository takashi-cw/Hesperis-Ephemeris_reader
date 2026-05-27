# runtests.jl — StellaBspReader 単体テスト・結合テスト
#
# 実行方法:
#   julia --project=. test/runtests.jl
#
#   結合テスト（de440s.bsp が必要）:
#   BSP_PATH="/Users/takashi/Stella_series/data/catalogs/de440s.bsp" julia --project=. test/runtests.jl
#
# BSP_PATH が未設定の場合は BSP 結合テストのみスキップされる。

using Test
using Printf

include(joinpath(@__DIR__, "..", "src", "StellaBspReader.jl"))
using .StellaBspReader

const BSP_PATH = get(ENV, "BSP_PATH", "")
const SKIP_BSP = isempty(BSP_PATH)

SKIP_BSP && @warn "BSP_PATH が未設定のため BSP 結合テストをスキップします"

# =============================================================================
# 合成 BSP データ生成ヘルパー
# =============================================================================

"""
    make_synthetic_bsp(; spk_type=2) -> Vector{UInt8}

最小限の DAF/SPK 合成バイナリを生成する（テスト用）。

Record 1 (0–1023):    ファイルヘッダー
Record 2 (1024–2047): サマリー（Sun=10 / SSB=0、指定の spk_type）
Record 3 (2048–):     Chebyshev データ

セグメントは J2000.0 ± 1 日をカバー。
J2000.0 での期待位置: x=100000, y=200000, z=300000 km
"""
function make_synthetic_bsp(; spk_type::Int32=Int32(2))::Vector{UInt8}
    buf = zeros(UInt8, 3 * 1024)

    write_le_f64!(offset, v) = begin
        bytes = reinterpret(UInt8, [Float64(v)])
        buf[offset+1:offset+8] .= bytes
    end
    write_le_i32!(offset, v) = begin
        bytes = reinterpret(UInt8, [Int32(v)])
        buf[offset+1:offset+4] .= bytes
    end
    write_str!(offset, s, pad) = begin
        for i in 1:pad
            buf[offset+i] = i <= length(s) ? UInt8(s[i]) : UInt8(' ')
        end
    end

    # Record 1: ヘッダー
    write_str!(0, "DAF/SPK ", 8)
    write_le_i32!(8,  2)       # ND
    write_le_i32!(12, 6)       # NI
    write_str!(16, "SyntheticBSP", 60)
    write_le_i32!(76, 2)       # FWARD
    write_le_i32!(80, 2)       # BWARD
    write_str!(88, "LTL-IEEE", 8)

    # Record 2: サマリー
    # Type 2: lastAddr=271 (rsize=11, n=1 → 15 doubles)
    # Type 3: lastAddr=268 (rsize=8,  n=1 → 12 doubles)
    is_type3  = spk_type == Int32(3)
    last_addr = is_type3 ? 268 : 271
    r2 = 1024
    write_le_f64!(r2 + 0,  0.0)   # next rec
    write_le_f64!(r2 + 8,  0.0)
    write_le_f64!(r2 + 16, 1.0)   # NSUM

    J2000 = 2_451_545.0
    s = r2 + 24
    write_le_f64!(s,      (J2000 - 1.0 - J2000) * 86400.0)  # start_sec
    write_le_f64!(s + 8,  (J2000 + 1.0 - J2000) * 86400.0)  # end_sec
    # ※ サマリーの start_sec / end_sec は J2000 からの秒。
    # 実際には _parse_summaries がそのまま格納するので、JD 値ではなく秒値を使う。
    # ただし _find_segment は (jd - J2000_JD) * SECS_PER_DAY と比較する。
    # J2000 ± 1 日 → ±86400 秒
    write_le_i32!(s + 16, 10)     # target = Sun
    write_le_i32!(s + 20, 0)      # center = SSB
    write_le_i32!(s + 24, 1)      # frame
    write_le_i32!(s + 28, Int32(spk_type))
    write_le_i32!(s + 32, 257)    # first_addr
    write_le_i32!(s + 36, last_addr)

    # Record 3: データ
    r3 = 2048
    write_le_f64!(r3 + 0,  0.0)       # mid (秒)
    write_le_f64!(r3 + 8,  86400.0)   # radius (秒)

    if is_type3
        # [Xpos, Ypos, Zpos, Xvel, Yvel, Zvel] ncoeff=1 each
        write_le_f64!(r3 + 16, 100000.0)  # Xpos[0]
        write_le_f64!(r3 + 24, 200000.0)  # Ypos[0]
        write_le_f64!(r3 + 32, 300000.0)  # Zpos[0]
        # Xvel, Yvel, Zvel = 0.0 (定数係数なので速度=0)
        # メタ: dataEnd = 268*8 = 2144 → metaOffset = 2144 - 32 = 2112
        write_le_f64!(2112, -86400.0)  # init
        write_le_f64!(2120, 172800.0)  # intlen (2日)
        write_le_f64!(2128, 8.0)       # rsize
        write_le_f64!(2136, 1.0)       # n
    else
        # [Xpos×3, Ypos×3, Zpos×3] ncoeff=3
        write_le_f64!(r3 + 16, 100000.0)  # Xpos[0]
        # Xpos[1], [2] = 0
        write_le_f64!(r3 + 40, 200000.0)  # Ypos[0]
        # Ypos[1], [2] = 0
        write_le_f64!(r3 + 64, 300000.0)  # Zpos[0]
        # Zpos[1], [2] = 0
        # メタ: dataEnd = 271*8 = 2168 → metaOffset = 2168 - 32 = 2136
        write_le_f64!(2136, -86400.0)  # init
        write_le_f64!(2144, 172800.0)  # intlen (2日)
        write_le_f64!(2152, 11.0)      # rsize
        write_le_f64!(2160, 1.0)       # n
    end

    return buf
end

# =============================================================================
# テスト
# =============================================================================

@testset "StellaBspReader" begin

    # ── Chebyshev 単体テスト ──────────────────────────────────────────

    @testset "chebyshev_eval" begin
        @test chebyshev_eval([5.0], 0.3) ≈ 5.0
        @test chebyshev_eval([2.0, 3.0], 0.5) ≈ 2.0 + 3.0 * 0.5
        # T2(x) = 2x^2 - 1 → [0, 0, 1] → 2*0.5^2 - 1 = -0.5
        @test chebyshev_eval([0.0, 0.0, 1.0], 0.5) ≈ -0.5
        @test chebyshev_eval(Float64[], 0.5) == 0.0
    end

    @testset "chebyshev_eval_with_deriv" begin
        coeffs = [2.0, 3.0, 1.0]
        x = 0.5
        (pos, dpdx) = chebyshev_eval_with_deriv(coeffs, x)
        @test pos ≈ chebyshev_eval(coeffs, x)
        eps_val = 1e-6
        numerical = (chebyshev_eval(coeffs, x + eps_val) - chebyshev_eval(coeffs, x - eps_val)) / (2eps_val)
        @test abs(dpdx - numerical) < 1e-9
    end

    @testset "chebyshev_eval3" begin
        cx = [1.0, 0.5]; cy = [2.0, 0.0]; cz = [0.0, 1.0]
        x = 0.5
        (px, py, pz) = chebyshev_eval3(cx, cy, cz, x)
        @test px ≈ chebyshev_eval(cx, x)
        @test py ≈ chebyshev_eval(cy, x)
        @test pz ≈ chebyshev_eval(cz, x)
    end

    # ── 合成 BSP — Type 2 ─────────────────────────────────────────────

    @testset "synthetic BSP Type 2" begin
        bsp = bsp_from_bytes(make_synthetic_bsp())
        J2000 = 2_451_545.0

        @test length(bsp.segments) == 1
        @test bsp.segments[1].target   == Int32(10)
        @test bsp.segments[1].spk_type == Int32(2)

        pos = get_position(bsp, Int32(10), Int32(0), J2000)
        @test isapprox(pos[1], 100000.0; atol=1e-6)
        @test isapprox(pos[2], 200000.0; atol=1e-6)
        @test isapprox(pos[3], 300000.0; atol=1e-6)

        (pos2, vel) = get_position_and_velocity(bsp, Int32(10), Int32(0), J2000)
        @test pos2 == pos
        @test isapprox(vel[1], 0.0; atol=1e-6)
        @test isapprox(vel[2], 0.0; atol=1e-6)
        @test isapprox(vel[3], 0.0; atol=1e-6)

        cpos = compute_position(bsp, Int32(10), Int32(0), J2000)
        @test isapprox(cpos[1], pos[1]; atol=1e-9)

        zero_pos = compute_position(bsp, Int32(10), Int32(10), J2000)
        @test zero_pos == (0.0, 0.0, 0.0)
    end

    @testset "out of coverage throws BspError" begin
        bsp = bsp_from_bytes(make_synthetic_bsp())
        J2000 = 2_451_545.0
        # セグメントは J2000 ± 1 日のみ → 7 日後は範囲外 → BspError が発生すること
        @test_throws BspError get_position(bsp, Int32(10), Int32(0), J2000 + 7.0)
        @test_throws BspError get_position(bsp, Int32(10), Int32(0), J2000 - 7.0)
    end

    @testset "out of coverage error message from _compute_chebyshev" begin
        # _compute_chebyshev を直接呼び、"out of coverage" メッセージを検証する。
        # _find_segment はセグメントの summary 範囲（±1 日）で事前フィルタするため、
        # 外部 API では "target not found" になってしまう。
        # ここではセグメントを直接渡してデータ記録の範囲チェックを確認する。
        bsp = bsp_from_bytes(make_synthetic_bsp())
        J2000 = 2_451_545.0
        seg = bsp.segments[1]
        # データ記録は init_epoch=-86400 〜 init_epoch+intlen=86400 → 7 日後は範囲外
        err = try
            StellaBspReader._compute_chebyshev(bsp, seg, J2000 + 7.0, false, 3)
            nothing
        catch e
            e
        end
        @test err isa StellaBspReader.BspError
        @test occursin("out of coverage", err.msg)
    end

    # ── 合成 BSP — Type 3 ─────────────────────────────────────────────

    @testset "synthetic BSP Type 3" begin
        bsp = bsp_from_bytes(make_synthetic_bsp(spk_type=Int32(3)))
        J2000 = 2_451_545.0

        @test bsp.segments[1].spk_type == Int32(3)

        pos = get_position(bsp, Int32(10), Int32(0), J2000)
        @test isapprox(pos[1], 100000.0; atol=1e-6)
        @test isapprox(pos[2], 200000.0; atol=1e-6)
        @test isapprox(pos[3], 300000.0; atol=1e-6)

        (pos2, vel) = get_position_and_velocity(bsp, Int32(10), Int32(0), J2000)
        @test pos2 == pos
        # 定数係数なので速度 = 0
        @test isapprox(vel[1], 0.0; atol=1e-6)
        @test isapprox(vel[2], 0.0; atol=1e-6)
        @test isapprox(vel[3], 0.0; atol=1e-6)

        cpos = compute_position(bsp, Int32(10), Int32(0), J2000)
        @test isapprox(cpos[1], pos[1]; atol=1e-9)
    end

    # ── Type 13（スコープ外）──────────────────────────────────────────

    @testset "unsupported Type 13 throws BspError" begin
        bsp = bsp_from_bytes(make_synthetic_bsp(spk_type=Int32(13)))
        J2000 = 2_451_545.0
        err = try
            get_position(bsp, Int32(10), Int32(0), J2000)
            nothing
        catch e
            e
        end
        @test err isa BspError
        @test occursin("unsupported SPK type: 13", err.msg)
    end

    # ── BSP 結合テスト（BSP_PATH 必須） ──────────────────────────────

    if !SKIP_BSP
        bsp = load_bsp(BSP_PATH)

        @testset "load_and_segments" begin
            @test !isempty(bsp.segments)
            println("ファイル名: ", bsp.name)
            println("セグメント数: ", length(bsp.segments))
            for seg in bsp.segments
                @printf("  target=%4d center=%4d type=%d [%.1f, %.1f]\n",
                    seg.target, seg.center, seg.spk_type, seg.start_sec, seg.end_sec)
            end
        end

        @testset "get_position_sun_j2000" begin
            jd = 2_451_545.0
            pos = compute_position(bsp, Naif.SUN, Naif.SSB, jd)
            println("Sun @ J2000.0 (km): x=$(pos[1]) y=$(pos[2]) z=$(pos[3])")
            dist = sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2)
            @test dist < 2_000_000.0
        end

        @testset "get_position_earth_j2000" begin
            jd = 2_451_545.0
            pos = compute_position(bsp, Naif.EARTH, Naif.SSB, jd)
            println("Earth @ J2000.0 (km): x=$(pos[1]) y=$(pos[2]) z=$(pos[3])")
            dist = sqrt(pos[1]^2 + pos[2]^2 + pos[3]^2)
            @test dist > 0.9 * AU_KM && dist < 1.1 * AU_KM
        end

        @testset "same_body_returns_zero" begin
            pos = compute_position(bsp, Naif.EARTH, Naif.EARTH, 2_451_545.0)
            @test pos == (0.0, 0.0, 0.0)
        end

        @testset "position_and_velocity" begin
            jd = 2_451_545.0
            (pos, vel) = get_position_and_velocity(bsp, Naif.SUN, Naif.SSB, jd)
            println("Sun vel (km/day): vx=$(vel[1]) vy=$(vel[2]) vz=$(vel[3])")
            speed = sqrt(vel[1]^2 + vel[2]^2 + vel[3]^2)
            @test speed > 0.0
            pos2 = get_position(bsp, Naif.SUN, Naif.SSB, jd)
            @test pos == pos2
        end
    end

end
