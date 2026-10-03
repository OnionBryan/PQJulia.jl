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
    ChaCha20(s, UInt64(s[13]) | UInt64(s[14]) << 32, UInt8[], 1)
end

rotl(x::UInt32, n) = (x << n) | (x >> (32 - n))

@inline function qround!(x, a, b, c, d)
    x[a] += x[b]; x[d] = rotl(x[d] ⊻ x[a], 16)
    x[c] += x[d]; x[b] = rotl(x[b] ⊻ x[c], 12)
    x[a] += x[b]; x[d] = rotl(x[d] ⊻ x[a], 8)
    x[c] += x[d]; x[b] = rotl(x[b] ⊻ x[c], 7)
end

function block!(r::ChaCha20)
    st = Vector{UInt32}(undef, 16)
    st[1:4] .= CW
    for i in 1:10; st[4+i] = r.s[i]; end
    st[15] = r.s[11] ⊻ (r.ctr % UInt32)
    st[16] = r.s[12] ⊻ UInt32(r.ctr >> 32)
    x = copy(st)
    for _ in 1:10
        qround!(x, 1, 5, 9, 13);  qround!(x, 2, 6, 10, 14)
        qround!(x, 3, 7, 11, 15); qround!(x, 4, 8, 12, 16)
        qround!(x, 1, 6, 11, 16); qround!(x, 2, 7, 12, 13)
        qround!(x, 3, 8, 9, 14);  qround!(x, 4, 5, 10, 15)
    end
    r.ctr += 1
    x .+= st
    wipe!(st)
    x
end

function refill!(r::ChaCha20)
    words = Vector{UInt32}(undef, 128)
    for i in 1:8
        blk = block!(r)
        words[i:8:end] = blk
        wipe!(blk)
    end
    wipe!(r.buf)
    r.buf = collect(reinterpret(UInt8, htol.(words))); r.pos = 1
    wipe!(words)
end

function randbytes(r::ChaCha20, k::Int)
    length(r.buf) - r.pos + 1 < k && refill!(r)
    out = r.buf[r.pos:r.pos+k-1]; r.pos += k
    out
end

wipe!(r::ChaCha20) = (r.s = ntuple(_ -> UInt32(0), Val(14)); r.ctr = 0; wipe!(r.buf); r)

end # module
