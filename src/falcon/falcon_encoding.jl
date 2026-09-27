# src/falcon/falcon_encoding.jl
# ============================================================================
# Falcon — signature compression (spec Algorithms 17/18, §3.11.2). Each
# coefficient s_i → 1 sign bit, its 7 low bits (MSB-first), then the high bits in
# unary (0^k 1, k = |s_i|>>7). Padded to a fixed bit length `slen`. Decompress
# inverts, with the spec's unique-encoding checks. Near-optimal for the Gaussian-
# distributed signature coefficients.
# ============================================================================
module FalconEncoding

export compress, decompress

# s: integer signature coefficients. slen: target length in BITS. Returns a
# Vector{UInt8} of slen÷8 bytes, or `nothing` if it overflows (caller resamples).
function compress(s::AbstractVector{<:Integer}, slen::Int)
    bits = Bool[]
    for si in s
        push!(bits, si < 0)                       # sign (1 = negative)
        a = abs(Int(si))
        for j in 6:-1:0; push!(bits, (a >> j) & 1 == 1); end   # 7 low bits, MSB first
        for _ in 1:(a >> 7); push!(bits, false); end           # unary high bits
        push!(bits, true)                                       # terminator
    end
    length(bits) > slen && return nothing
    while length(bits) < slen; push!(bits, false); end          # pad with zeros
    return _pack(bits)
end

# str: Vector{UInt8} of slen÷8 bytes. n: number of coefficients. Returns the
# coefficients, or `nothing` on an invalid/non-unique encoding.
function decompress(str::AbstractVector{UInt8}, slen::Int, n::Int)
    length(str) * 8 == slen || return nothing
    bits = _unpack(str); s = Int[]; i = 1
    for _ in 1:n
        i + 7 > slen && return nothing
        sign = bits[i]; i += 1
        low = 0; for _ in 1:7; low = (low << 1) | (bits[i] ? 1 : 0); i += 1; end
        k = 0; while i <= slen && !bits[i]; k += 1; i += 1; end
        i > slen && return nothing                 # missing terminator
        i += 1                                      # consume the terminator '1'
        val = (k << 7) + low
        sign && val == 0 && return nothing          # reject -0 (non-unique)
        push!(s, sign ? -val : val)
    end
    all(!, @view bits[i:end]) || return nothing     # trailing bits must be 0
    return s
end

function _pack(bits::Vector{Bool})
    n = length(bits) ÷ 8; out = zeros(UInt8, n)
    for b in 1:n, j in 1:8
        bits[(b-1)*8 + j] && (out[b] |= UInt8(1) << (8 - j))   # MSB-first within a byte
    end
    out
end
function _unpack(bytes::AbstractVector{UInt8})
    bits = Vector{Bool}(undef, length(bytes) * 8)
    for (b, byte) in enumerate(bytes), j in 1:8
        bits[(b-1)*8 + j] = (byte >> (8 - j)) & 1 == 1
    end
    bits
end

end # module
