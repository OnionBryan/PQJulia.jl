# src/falcon/falcon_chacha.jl
# ============================================================================
# Falcon — the ChaCha20 PRNG the reference signer feeds to SamplerZ (tprest
# falcon.py rng.py, which matches the round-3 C code). The 56-byte seed is 14
# little-endian words s[0..13]; ctr = s[12] | s[13] << 32. Each block is the
# RFC 7539 ChaCha20 core with state rows (CW, s[0..3], s[4..7], s[8], s[9],
# s[10] ⊻ ctr_lo, s[11] ⊻ ctr_hi). Eight consecutive blocks are word-interleaved
# into a 512-byte buffer; a request that does not fit discards the remainder.
# ============================================================================
module FalconChaCha

import ...Wipe: wipe!

export ChaCha20, randbytes

const CW = (0x61707865, 0x3320646e, 0x79622d32, 0x6b206574)

mutable struct ChaCha20
    s::NTuple{14,UInt32}
    ctr::UInt64
    buf::Vector{UInt8}
    pos::Int
end

function ChaCha20(seed::AbstractVector{UInt8})
    length(seed) == 56 || throw(ArgumentError("ChaCha20 seed must be 56 bytes"))
    s = ntuple(i -> UInt32(seed[4i-3]) | UInt32(seed[4i-2]) << 8 |
                    UInt32(seed[4i-1]) << 16 | UInt32(seed[4i]) << 24, Val(14))
    ChaCha20(s, UInt64(s[13]) | UInt64(s[14]) << 32, zeros(UInt8, 512), 513)   # empty: refill on first read
end

rotl(x::UInt32, n) = (x << n) | (x >> (32 - n))

@inline function qround(a::UInt32, b::UInt32, c::UInt32, d::UInt32)
    a += b; d = rotl(d ⊻ a, 16)
    c += d; b = rotl(b ⊻ c, 12)
    a += b; d = rotl(d ⊻ a, 8)
    c += d; b = rotl(b ⊻ c, 7)
    a, b, c, d
end

# One ChaCha20 block on the state words in registers (no heap copy of the keystream state).
function block(s::NTuple{14,UInt32}, ctr::UInt64)
    st = (CW..., s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10],
          s[11] ⊻ (ctr % UInt32), s[12] ⊻ UInt32(ctr >> 32))
    x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16 = st
    for _ in 1:10
        x1, x5, x9, x13 = qround(x1, x5, x9, x13);  x2, x6, x10, x14 = qround(x2, x6, x10, x14)
        x3, x7, x11, x15 = qround(x3, x7, x11, x15); x4, x8, x12, x16 = qround(x4, x8, x12, x16)
        x1, x6, x11, x16 = qround(x1, x6, x11, x16); x2, x7, x12, x13 = qround(x2, x7, x12, x13)
        x3, x8, x9, x14 = qround(x3, x8, x9, x14);  x4, x5, x10, x15 = qround(x4, x5, x10, x15)
    end
    map(+, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16), st)
end

# Eight blocks, word-interleaved (word w of block i is word 8(w−1) + i), little-endian, written
# over the previous 512 bytes in place.
function refill!(r::ChaCha20)
    buf = r.buf
    checkbounds(buf, 1:512)
    for i in 1:8
        blk = block(r.s, r.ctr); r.ctr += 1
        @inbounds for w in 1:16
            o = 4 * (8 * (w - 1) + i - 1)
            v = blk[w]
            buf[o+1] = v % UInt8; buf[o+2] = (v >> 8) % UInt8
            buf[o+3] = (v >> 16) % UInt8; buf[o+4] = (v >> 24) % UInt8
        end
    end
    r.pos = 1
end

function randbytes(r::ChaCha20, k::Int)
    length(r.buf) - r.pos + 1 < k && refill!(r)
    out = r.buf[r.pos:r.pos+k-1]; r.pos += k
    out
end

# randbytes(r, 1)[1] and the little-endian value of randbytes(r, 9), without the copies.
function randbyte(r::ChaCha20)
    length(r.buf) - r.pos + 1 < 1 && refill!(r)
    v = @inbounds r.buf[r.pos]; r.pos += 1
    v
end
function randu72(r::ChaCha20)
    length(r.buf) - r.pos + 1 < 9 && refill!(r)
    u = UInt128(0)
    @inbounds for i in 0:8
        u |= UInt128(r.buf[r.pos+i]) << (8i)
    end
    r.pos += 9
    u
end

wipe!(r::ChaCha20) = (r.s = ntuple(_ -> UInt32(0), Val(14)); r.ctr = 0; wipe!(r.buf); r)

end # module
