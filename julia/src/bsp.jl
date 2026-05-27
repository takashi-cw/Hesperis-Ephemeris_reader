# bsp.jl — JPL BSP/SPK バイナリ読み取り
#
# NASA NAIF DAF/SPK フォーマットを解析し、
# 指定天体の ICRS XYZ 位置ベクトル（km）を返す。
#
# 対応 SPK セグメントタイプ:
#   Type 2  — Chebyshev 多項式（位置のみ）        ← de440s.bsp / de440.bsp / de441.bsp
#   Type 3  — Chebyshev 多項式（位置＋速度）       ← de440.bsp / de441.bsp の月秤動角セグメント
#   Type 13 — Hermite 補間（小天体）               ← スコープ外（UnsupportedSPKType を送出）
#
# 出典フォーマット仕様:
#   - NAIF SPK Required Reading (NAIF N0067)
#   - NAIF DAF Required Reading (NAIF N0067)
#   - jplephem (Brandon Rhodes, MIT License) の設計を参考に Julia で再実装
#
# ライセンス: MIT

const _RECORD_SIZE = 1024
const _SPK_TYPE_2  = Int32(2)
const _SPK_TYPE_3  = Int32(3)

# --- エラー型 ---

struct BspError <: Exception
    msg::String
end
Base.showerror(io::IO, e::BspError) = print(io, "BspError: ", e.msg)

# --- セグメント ---

"""パース済みセグメントのメタデータ"""
struct BspSegment
    target::Int32
    center::Int32
    spk_type::Int32
    "セグメント開始時刻（J2000 からの秒）"
    start_sec::Float64
    "セグメント終了時刻（J2000 からの秒）"
    end_sec::Float64
    "データ開始アドレス（1-indexed double オフセット）"
    start_idx::Int
    "データ終了アドレス（1-indexed double オフセット）"
    end_idx::Int
end

# --- BspFile ---

"""パース済み .bsp ファイルのラッパー"""
struct BspFile
    "LOCIFN（ファイル内部名）"
    name::String
    "セグメント一覧"
    segments::Vector{BspSegment}
    _data::Vector{UInt8}
    _is_le::Bool
end

# --- ロード ---

"""
    load_bsp(path::String) -> BspFile

.bsp ファイルをパスから読み込む。
"""
function load_bsp(path::String)::BspFile
    data = read(path)
    _bsp_from_bytes(data)
end

"""
    bsp_from_bytes(data::Vector{UInt8}) -> BspFile

バイト列から直接初期化する（テスト用）。
"""
function bsp_from_bytes(data::Vector{UInt8})::BspFile
    _bsp_from_bytes(data)
end

function _bsp_from_bytes(data::Vector{UInt8})::BspFile
    (nd, ni, first_sum_rec, is_le, name) = _parse_file_record(data)
    segments = _parse_summaries(data, nd, ni, first_sum_rec, is_le)
    return BspFile(name, segments, data, is_le)
end

# --- 公開 API ---

"""
    get_position(bsp, target, center, jd_tdb) -> NTuple{3,Float64}

ICRS 位置ベクトル (x, y, z) km を返す（直接セグメントのみ）。
セグメントが見つからない場合、または JD がカバー範囲外の場合は BspError を送出する。
"""
function get_position(
    bsp::BspFile,
    target::Int32,
    center::Int32,
    jd_tdb::Float64,
)::NTuple{3,Float64}
    seg = _find_segment(bsp, target, center, jd_tdb)
    isnothing(seg) && throw(BspError("target NAIF ID $target not found"))
    components = _spk_components(seg.spk_type)
    isnothing(components) && throw(BspError("unsupported SPK type: $(seg.spk_type)"))
    return _compute_chebyshev(bsp, seg, jd_tdb, false, components)[1]
end

"""
    get_position_and_velocity(bsp, target, center, jd_tdb)
        -> (position::NTuple{3,Float64}, velocity::NTuple{3,Float64})

ICRS 位置（km）と速度（km/day）を返す（直接セグメントのみ）。
"""
function get_position_and_velocity(
    bsp::BspFile,
    target::Int32,
    center::Int32,
    jd_tdb::Float64,
)::Tuple{NTuple{3,Float64},NTuple{3,Float64}}
    seg = _find_segment(bsp, target, center, jd_tdb)
    isnothing(seg) && throw(BspError("target NAIF ID $target not found"))
    components = _spk_components(seg.spk_type)
    isnothing(components) && throw(BspError("unsupported SPK type: $(seg.spk_type)"))
    return _compute_chebyshev(bsp, seg, jd_tdb, true, components)
end

"""
    compute_position(bsp, target, center, jd_tdb) -> NTuple{3,Float64}

セグメントチェーンを辿って任意の中心座標で位置を合成する。

DE440s のセグメント構成例:
  SSB(0) → Sun(10)
  SSB(0) → EMB(3) → Earth(399)
  SSB(0) → EMB(3) → Moon(301)
"""
function compute_position(
    bsp::BspFile,
    target::Int32,
    center::Int32,
    jd_tdb::Float64,
)::NTuple{3,Float64}
    target == center && return (0.0, 0.0, 0.0)

    if !isnothing(_find_segment(bsp, target, center, jd_tdb))
        return get_position(bsp, target, center, jd_tdb)
    end

    pt = _pos_from_ssb(bsp, target, jd_tdb)
    pc = center == Int32(0) ? (0.0, 0.0, 0.0) : _pos_from_ssb(bsp, center, jd_tdb)
    return (pt[1] - pc[1], pt[2] - pc[2], pt[3] - pc[3])
end

# --- プライベート：セグメント検索 ---

function _find_segment(
    bsp::BspFile,
    target::Int32,
    center::Int32,
    jd_tdb::Float64,
)::Union{BspSegment,Nothing}
    t_sec = (jd_tdb - J2000_JD) * SECS_PER_DAY
    for seg in bsp.segments
        if seg.target == target &&
           seg.center == center &&
           t_sec >= seg.start_sec &&
           t_sec <= seg.end_sec
            return seg
        end
    end
    return nothing
end

function _pos_from_ssb(bsp::BspFile, target::Int32, jd_tdb::Float64)::NTuple{3,Float64}
    t_sec = (jd_tdb - J2000_JD) * SECS_PER_DAY

    seg = _find_segment(bsp, target, Int32(0), jd_tdb)
    if !isnothing(seg)
        components = _spk_components(seg.spk_type)
        isnothing(components) && throw(BspError("unsupported SPK type: $(seg.spk_type)"))
        return _compute_chebyshev(bsp, seg, jd_tdb, false, components)[1]
    end

    # 中間天体経由で合成
    for seg in bsp.segments
        seg.target != target && continue
        t_sec < seg.start_sec || t_sec > seg.end_sec && continue
        components = _spk_components(seg.spk_type)
        isnothing(components) && throw(BspError("unsupported SPK type: $(seg.spk_type)"))
        from_center = _compute_chebyshev(bsp, seg, jd_tdb, false, components)[1]
        center_from_ssb = seg.center == Int32(0) ? (0.0, 0.0, 0.0) : _pos_from_ssb(bsp, seg.center, jd_tdb)
        return (
            center_from_ssb[1] + from_center[1],
            center_from_ssb[2] + from_center[2],
            center_from_ssb[3] + from_center[3],
        )
    end

    throw(BspError("target NAIF ID $target not found"))
end

# --- プライベート：SPK タイプヘルパー ---

"""
    _spk_components(spk_type) -> Union{Int,Nothing}

Type 2 → 3（位置のみ）、Type 3 → 6（位置＋速度）、それ以外 → nothing。
"""
function _spk_components(spk_type::Int32)::Union{Int,Nothing}
    spk_type == _SPK_TYPE_2 && return 3
    spk_type == _SPK_TYPE_3 && return 6
    return nothing
end

# --- プライベート：Chebyshev 計算（Type 2 / Type 3 共通） ---

"""
    _compute_chebyshev(bsp, seg, jd_tdb, with_velocity, components)

Type 2（components=3）および Type 3（components=6）共通の Chebyshev 評価。
jd_tdb がセグメントのカバー範囲外の場合は BspError を送出する。
"""
function _compute_chebyshev(
    bsp::BspFile,
    seg::BspSegment,
    jd_tdb::Float64,
    with_velocity::Bool,
    components::Int,
)::Tuple{NTuple{3,Float64},NTuple{3,Float64}}
    data     = bsp._data
    is_le    = bsp._is_le

    data_start = (seg.start_idx - 1) * 8
    data_end   = seg.end_idx * 8

    # メタデータ（データ末尾 4 double = 32 bytes）
    meta_offset = data_end - 32
    init_epoch  = _read_f64(data, meta_offset,      is_le)
    intlen      = _read_f64(data, meta_offset + 8,  is_le)
    rsize       = round(Int, _read_f64(data, meta_offset + 16, is_le))
    n           = round(Int, _read_f64(data, meta_offset + 24, is_le))

    t_seconds = (jd_tdb - J2000_JD) * SECS_PER_DAY

    # 範囲外チェック（サイレントクランプ禁止）
    t_min = init_epoch
    t_max = init_epoch + n * intlen
    if t_seconds < t_min || t_seconds > t_max
        throw(BspError(
            "out of coverage: t=$(t_seconds)s (valid $(t_min)s – $(t_max)s)" *
            " for target=$(seg.target)"
        ))
    end

    idx = clamp(floor(Int, (t_seconds - init_epoch) / intlen), 0, n - 1)

    rec_offset = data_start + idx * rsize * 8

    mid    = _read_f64(data, rec_offset,     is_le)
    radius = _read_f64(data, rec_offset + 8, is_le)
    x      = (t_seconds - mid) / radius

    ncoeff = (rsize - 2) ÷ components
    base   = rec_offset + 16

    cx = [_read_f64(data, base + i * 8,                  is_le) for i in 0:ncoeff-1]
    cy = [_read_f64(data, base + ncoeff * 8 + i * 8,     is_le) for i in 0:ncoeff-1]
    cz = [_read_f64(data, base + ncoeff * 8 * 2 + i * 8, is_le) for i in 0:ncoeff-1]

    if with_velocity
        if components == 6
            # Type 3: 速度係数を直接評価（km/s → km/day に変換）
            vx_coeffs = [_read_f64(data, base + ncoeff * 8 * 3 + i * 8, is_le) for i in 0:ncoeff-1]
            vy_coeffs = [_read_f64(data, base + ncoeff * 8 * 4 + i * 8, is_le) for i in 0:ncoeff-1]
            vz_coeffs = [_read_f64(data, base + ncoeff * 8 * 5 + i * 8, is_le) for i in 0:ncoeff-1]
            pos = chebyshev_eval3(cx, cy, cz, x)
            vel = (
                chebyshev_eval(vx_coeffs, x) * SECS_PER_DAY,
                chebyshev_eval(vy_coeffs, x) * SECS_PER_DAY,
                chebyshev_eval(vz_coeffs, x) * SECS_PER_DAY,
            )
            return (pos, vel)
        else
            # Type 2: 位置多項式を微分して速度を求める
            interval_days = radius * 2.0 / SECS_PER_DAY
            return chebyshev_eval3_with_velocity(cx, cy, cz, x, interval_days)
        end
    else
        return (chebyshev_eval3(cx, cy, cz, x), (0.0, 0.0, 0.0))
    end
end

# --- プライベート：バイナリ読み込みヘルパー ---

# offset は 0-indexed（Julia の配列は 1-indexed なので内部で +1 する）

@inline function _read_f64(data::Vector{UInt8}, offset::Int, is_le::Bool)::Float64
    raw = reinterpret(UInt64, data[offset+1:offset+8])[1]
    reinterpret(Float64, is_le ? raw : bswap(raw))
end

@inline function _read_i32(data::Vector{UInt8}, offset::Int, is_le::Bool)::Int32
    raw = reinterpret(UInt32, data[offset+1:offset+4])[1]
    reinterpret(Int32, is_le ? raw : bswap(raw))
end

@inline function _read_str(data::Vector{UInt8}, offset::Int, len::Int)::String
    String(data[offset+1:offset+len]) |> strip
end

# --- ファイルレコード解析 ---

function _parse_file_record(
    data::Vector{UInt8},
)::Tuple{Int,Int,Int,Bool,String}

    length(data) < _RECORD_SIZE && throw(BspError("ファイルサイズが小さすぎます"))

    locidw = String(data[1:8])
    if !startswith(locidw, "DAF/SPK") && !startswith(locidw, "DAF/EK")
        throw(BspError("LOCIDW=\"$(strip(locidw))\""))
    end

    locfmt = strip(String(data[89:96]))  # offset 88 (0-indexed) → index 89 (1-indexed)
    is_le  = locfmt != "BIG-IEEE"

    nd            = Int(_read_i32(data, 8,  is_le))
    ni            = Int(_read_i32(data, 12, is_le))
    first_sum_rec = Int(_read_i32(data, 76, is_le))
    name          = _read_str(data, 16, 60)

    return (nd, ni, first_sum_rec, is_le, name)
end

# --- サマリーレコード解析 ---

function _parse_summaries(
    data::Vector{UInt8},
    nd::Int,
    ni::Int,
    first_sum_rec::Int,
    is_le::Bool,
)::Vector{BspSegment}

    summary_doubles = nd + (ni + 1) ÷ 2
    summary_bytes   = summary_doubles * 8

    segments = BspSegment[]
    rec_num  = first_sum_rec

    while rec_num > 0
        rec_offset  = (rec_num - 1) * _RECORD_SIZE
        next_rec    = round(Int, _read_f64(data, rec_offset,      is_le))
        n_summaries = round(Int, _read_f64(data, rec_offset + 16, is_le))

        for i in 0:n_summaries-1
            base = rec_offset + 24 + i * summary_bytes

            start_sec = _read_f64(data, base,     is_le)
            end_sec   = _read_f64(data, base + 8, is_le)

            int_base   = base + nd * 8
            target     = _read_i32(data, int_base,      is_le)
            center     = _read_i32(data, int_base + 4,  is_le)
            spk_type   = _read_i32(data, int_base + 12, is_le)
            first_addr = Int(_read_i32(data, int_base + 16, is_le))
            last_addr  = Int(_read_i32(data, int_base + 20, is_le))

            push!(segments, BspSegment(
                target, center, spk_type,
                start_sec, end_sec,
                first_addr, last_addr,
            ))
        end

        rec_num = next_rec
    end

    return segments
end
