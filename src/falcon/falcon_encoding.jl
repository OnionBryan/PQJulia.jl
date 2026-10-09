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
# Bits are written MSB-first straight into the output bytes; bit p (1-based) is the
# (8 - (p-1)%8)-th bit of byte ⌈p/8⌉.
function compress(s::AbstractVector{<:Integer}, slen::Int)
    out = zeros(UInt8, slen ÷ 8); nbits = 8 * length(out); p = 0
    for si in s
        a = abs(Int(si))
        p + 9 + (a >> 7) > slen && return nothing             # sign, 7 low bits, unary, terminator
        si < 0 && _setbit!(out, p + 1, nbits)
        for j in 6:-1:0
            (a >> j) & 1 == 1 && _setbit!(out, p + 8 - j, nbits)
        end
        p += 8 + (a >> 7) + 1                                  # unary zeros are already zero
        _setbit!(out, p, nbits)                                # terminator
    end
    return out                                                 # zero padding up to slen
end
@inline _setbit!(out, p, nbits) = 1 <= p <= nbits && (@inbounds out[(p + 7) >> 3] |= 0x80 >> ((p - 1) & 7))

# Bit i (1-based, MSB-first) of str.
@inline _bit(str, o, i) = @inbounds (str[o + ((i - 1) >> 3)] >> (7 - ((i - 1) & 7))) & 0x01 == 0x01

# str: Vector{UInt8} of slen÷8 bytes. n: number of coefficients. Returns the
# coefficients, or `nothing` on an invalid/non-unique encoding.
function decompress(str::AbstractVector{UInt8}, slen::Int, n::Int)
    length(str) * 8 == slen || return nothing
    o = firstindex(str); s = Vector{Int}(undef, n); i = 1
    for c in 1:n
        i + 7 > slen && return nothing
        sign = _bit(str, o, i); i += 1
        low = 0; for _ in 1:7; low = (low << 1) | Int(_bit(str, o, i)); i += 1; end
        k = 0; while i <= slen && !_bit(str, o, i); k += 1; i += 1; end
        i > slen && return nothing                 # missing terminator
        i += 1                                      # consume the terminator '1'
        val = (k << 7) + low
        sign && val == 0 && return nothing          # reject -0 (non-unique)
        s[c] = sign ? -val : val
    end
    for j in i:slen                                 # trailing bits must be 0
        _bit(str, o, j) && return nothing
    end
    return s
end

end # module
