"""
    kyber.jl -- Core Kyber (CRYSTALS-Kyber) number-theoretic transform math,
    CBD sampling, poly operations, compress/decompress, encode/decode,
    and constant-time utilities.

    Port of pq-crystals/kyber/ref: reduce.c, ntt.c, cbd.c, poly.c, verify.c
    Constants from params.h, reduce.h, ntt.c

    All arithmetic follows the C reference exactly, using `% Int16` for
    truncation (never `Int16()` which throws on overflow).
"""

using SHA

# ── Constants ────────────────────────────────────────────────────────────────

const KYBER_N    = 256
const KYBER_Q    = Int16(3329)
const KYBER_MONT = Int16(-1044)   # 2^16 mod q  (centered representation)
const KYBER_QINV = Int16(-3327)   # q^{-1} mod 2^16  (signed)

# Precomputed zetas for NTT, hardcoded from pq-crystals/kyber/ref/ntt.c:39-56.
# 128 entries; Julia 1-indexed so zetas[1] = C's zetas[0] = -1044.
const KYBER_ZETAS = Int16[
    -1044,  -758,  -359, -1517,  1493,  1422,   287,   202,
     -171,   622,  1577,   182,   962, -1202, -1474,  1468,
      573, -1325,   264,   383,  -829,  1458, -1602,  -130,
     -681,  1017,   732,   608, -1542,   411,  -205, -1571,
     1223,   652,  -552,  1015, -1293,  1491,  -282, -1544,
      516,    -8,  -320,  -666, -1618, -1162,   126,  1469,
     -853,   -90,  -271,   830,   107, -1421,  -247,  -951,
     -398,   961, -1508,  -725,   448, -1065,   677, -1275,
    -1103,   430,   555,   843, -1251,   871,  1550,   105,
      422,   587,   177,  -235,  -291,  -460,  1574,  1653,
     -246,   778,  1159,  -147,  -777,  1483,  -602,  1119,
    -1590,   644,  -872,   349,   418,   329,  -156,   -75,
      817,  1097,   603,   610,  1322, -1285, -1465,   384,
    -1215,  -136,  1218, -1335,  -874,   220, -1187, -1659,
    -1185, -1530, -1278,   794, -1510,  -854,  -870,   478,
     -108,  -308,   996,   991,   958, -1460,  1522,  1628
]

# ── Montgomery Reduce ────────────────────────────────────────────────────────
# Input:  a in Int32   (product of two Int16 values * R)
# Output: a * R^{-1} mod q,  in {-(q-1), ..., q-1}
# Mirrors pq-crystals/kyber/ref/reduce.c:16-23.
#
# The C reference uses Int32 intermediates and documents the input constraint
# {-q*2^15, ..., q*2^15-1}. However, fqmul passes Int32(a)*Int32(b) where
# a,b are arbitrary Int16, producing products up to ~10^9 which exceeds that
# range. The C code relies on signed overflow wrapping on two's-complement
# machines (technically UB), which Julia does not allow.
#
# We widen to Int64 for the subtraction so no overflow occurs.  After the
# >>16 shift the exact quotient u = (a - t*q)/2^16 satisfies |u| <= 34432
# (since |a| <= 2^31 and |t*q| <= 32768*3329).  When the input IS in the
# valid range, |u| < q and rem is a no-op.  For out-of-range inputs the
# rem(u, q) folds u back into {-(q-1), ..., q-1}, giving a mathematically
# correct a * R^{-1} mod q for ALL Int32 inputs.

function montgomery_reduce(a::Int32)::Int16
    t = (a % Int16) * KYBER_QINV              # truncate to Int16, multiply (wraps)
    u = (Int64(a) - Int64(t) * Int64(KYBER_Q)) >> 16
    return rem(u, Int64(KYBER_Q)) % Int16
end

# ── fqmul ────────────────────────────────────────────────────────────────────
# Multiplication followed by Montgomery reduction.

@inline function fqmul(a::Int16, b::Int16)::Int16
    return montgomery_reduce(Int32(a) * Int32(b))
end

# ── Barrett Reduce ───────────────────────────────────────────────────────────
# Input:  a in Int16
# Output: a mod q, centered in {-(q-1)/2, ..., (q-1)/2}
# Mirrors pq-crystals/kyber/ref/reduce.c:35-42.

function barrett_reduce(a::Int16)::Int16
    v = Int16(20159)   # floor((1<<26 + q/2) / q)
    t = ((Int32(v) * Int32(a) + (Int32(1) << 25)) >> 26) % Int16
    t = (Int32(t) * Int32(KYBER_Q)) % Int16
    return a - t
end

# ── Forward NTT (in-place) ───────────────────────────────────────────────────
# Mirrors pq-crystals/kyber/ref/ntt.c:80-95.
# C starts k=1 (post-increment); Julia 1-indexed so k starts at 2.
# No Barrett reduce in the forward direction.

function kyber_ntt!(r::Vector{Int16})
    k = 2   # Julia's zetas[2] = C's zetas[1]
    len = 128
    while len >= 2
        start = 0
        while start < 256
            zeta = KYBER_ZETAS[k]
            k += 1
            for j in start:(start + len - 1)
                t = fqmul(zeta, r[j + len + 1])    # +1 for 1-indexing
                r[j + len + 1] = r[j + 1] - t
                r[j + 1]       = r[j + 1] + t
            end
            start += 2 * len
        end
        len >>= 1
    end
    return r
end

# ── Inverse NTT (in-place) ──────────────────────────────────────────────────
# Mirrors pq-crystals/kyber/ref/ntt.c:106-126.
# C starts k=127 (post-decrement); Julia 1-indexed so k starts at 128.
# Barrett reduce on the addition.  Final scaling by f = 1441 = mont^2/128.

function kyber_invntt!(r::Vector{Int16})
    f = Int16(1441)
    k = 128  # Julia's zetas[128] = C's zetas[127]
    len = 2
    while len <= 128
        start = 0
        while start < 256
            zeta = KYBER_ZETAS[k]
            k -= 1
            for j in start:(start + len - 1)
                t = r[j + 1]
                r[j + 1]       = barrett_reduce(t + r[j + len + 1])
                r[j + len + 1] = r[j + len + 1] - t
                r[j + len + 1] = fqmul(zeta, r[j + len + 1])
            end
            start += 2 * len
        end
        len <<= 1
    end
    for j in 1:256
        r[j] = fqmul(r[j], f)
    end
    return r
end

# ── Basemul ──────────────────────────────────────────────────────────────────
# Multiply two degree-1 polynomials in Z_q[X]/(X^2 - zeta).
# Mirrors pq-crystals/kyber/ref/ntt.c:139-146.
# r, a, b are length-2 views (or any AbstractVector{Int16}).

function kyber_basemul!(r::AbstractVector{Int16}, a::AbstractVector{Int16},
                        b::AbstractVector{Int16}, zeta::Int16)
    r[1]  = fqmul(a[2], b[2])
    r[1]  = fqmul(r[1], zeta)
    r[1] += fqmul(a[1], b[1])
    r[2]  = fqmul(a[1], b[2])
    r[2] += fqmul(a[2], b[1])
    return r
end

# ── Additional Constants ────────────────────────────────────────────────────

const KYBER_SYMBYTES = 32
const KYBER_POLYBYTES = 384       # 256 * 12 / 8

# ── CBD — Centered Binomial Distribution ────────────────────────────────────
# Port of pq-crystals/kyber/ref/cbd.c
# cbd2: eta=2, 128 bytes → 256 coefficients in [-2, 2]
# cbd3: eta=3, 192 bytes → 256 coefficients in [-3, 3]

"""
    kyber_cbd2!(r, buf)

Centered binomial distribution with eta=2.
`buf` must have at least 128 bytes (2*KYBER_N/4).
`r` must be a Vector{Int16} of length >= 256.
Produces coefficients in [-2, 2].
"""
function kyber_cbd2!(r::Vector{Int16}, buf::Vector{UInt8})
    for i in 0:(KYBER_N ÷ 8 - 1)   # 0:31
        t = UInt32(buf[4i + 1]) |
            (UInt32(buf[4i + 2]) << 8) |
            (UInt32(buf[4i + 3]) << 16) |
            (UInt32(buf[4i + 4]) << 24)
        d = t & 0x55555555
        d += (t >> 1) & 0x55555555
        for j in 0:7
            a = (d >> (4j)) & 0x3
            b = (d >> (4j + 2)) & 0x3
            r[8i + j + 1] = (a % Int16) - (b % Int16)
        end
    end
    return r
end

"""
    kyber_cbd3!(r, buf)

Centered binomial distribution with eta=3.
`buf` must have at least 192 bytes (3*KYBER_N/4).
`r` must be a Vector{Int16} of length >= 256.
Produces coefficients in [-3, 3].
"""
function kyber_cbd3!(r::Vector{Int16}, buf::Vector{UInt8})
    for i in 0:(KYBER_N ÷ 4 - 1)   # 0:63
        t = UInt32(buf[3i + 1]) |
            (UInt32(buf[3i + 2]) << 8) |
            (UInt32(buf[3i + 3]) << 16)
        d = t & 0x00249249
        d += (t >> 1) & 0x00249249
        d += (t >> 2) & 0x00249249
        for j in 0:3
            a = (d >> (6j)) & 0x7
            b = (d >> (6j + 3)) & 0x7
            r[4i + j + 1] = (a % Int16) - (b % Int16)
        end
    end
    return r
end

# ── PRF / XOF / Hash Wrappers ──────────────────────────────────────────────
# Using Julia's stdlib SHA (SHAKE-128, SHAKE-256, SHA3-256, SHA3-512).

"""PRF: SHAKE-256(seed || nonce, outlen)"""
function kyber_prf(seed::Vector{UInt8}, nonce::UInt8, outlen::Int)
    return SHA.shake256(vcat(seed, UInt8[nonce]), UInt64(outlen))
end

"""XOF: SHAKE-128(seed, outlen)"""
function kyber_xof(seed::Vector{UInt8}, outlen::Int)
    return SHA.shake128(seed, UInt64(outlen))
end

"""Hash H: SHA3-256"""
function kyber_hash_h(data::Vector{UInt8})
    return SHA.sha3_256(data)
end

"""Hash G: SHA3-512"""
function kyber_hash_g(data::Vector{UInt8})
    return SHA.sha3_512(data)
end

# ── Uniform Sampling (Rejection Sampling) ──────────────────────────────────
# Port of pq-crystals/kyber/ref/indcpa.c:rej_uniform
# 3 bytes → two 12-bit values, reject >= q.

"""
    kyber_rej_uniform!(r, buf) -> ctr

Rejection-sample uniform coefficients mod q from random bytes.
Returns the number of valid coefficients produced.
"""
function kyber_rej_uniform!(r::Vector{Int16}, buf::Vector{UInt8})
    len = length(r)
    buflen = length(buf)
    ctr = 0
    pos = 1  # Julia 1-indexed
    while ctr < len && pos + 2 <= buflen
        val0 = (UInt16(buf[pos]) | (UInt16(buf[pos + 1]) << 8)) & 0x0FFF
        val1 = ((UInt16(buf[pos + 1]) >> 4) | (UInt16(buf[pos + 2]) << 4)) & 0x0FFF
        pos += 3
        if val0 < UInt16(KYBER_Q)
            ctr += 1
            r[ctr] = val0 % Int16
        end
        if ctr < len && val1 < UInt16(KYBER_Q)
            ctr += 1
            r[ctr] = val1 % Int16
        end
    end
    return ctr
end

"""
    kyber_sample_uniform!(r, seed, i, j)

Sample a full polynomial (256 coeffs) uniformly mod q using SHAKE-128
seeded with seed||i||j.  Uses rejection sampling.
"""
function kyber_sample_uniform!(r::Vector{Int16}, seed::Vector{UInt8},
                               i::UInt8, j::UInt8)
    # XOF_BLOCKBYTES = 168 for SHAKE-128; generate enough
    # GEN_MATRIX_NBLOCKS ~ 3 blocks for Kyber-768
    xof_input = vcat(seed, UInt8[i, j])
    buf = kyber_xof(xof_input, 168 * 4)  # generous buffer
    ctr = kyber_rej_uniform!(r, buf)
    # In extremely rare cases we might need more bytes
    while ctr < KYBER_N
        extra = kyber_xof(xof_input, 168 * 8)
        tmp = Vector{Int16}(undef, KYBER_N - ctr)
        got = kyber_rej_uniform!(tmp, extra)
        for k in 1:min(got, KYBER_N - ctr)
            r[ctr + k] = tmp[k]
        end
        ctr += min(got, KYBER_N - ctr)
    end
    return r
end

# ── Poly Arithmetic Operations ─────────────────────────────────────────────
# Port of pq-crystals/kyber/ref/poly.c: poly_add, poly_sub, poly_reduce, poly_tomont

"""Coefficient-wise addition, no modular reduction."""
function kyber_poly_add!(r::Vector{Int16}, a::Vector{Int16}, b::Vector{Int16})
    @inbounds for i in 1:KYBER_N
        r[i] = a[i] + b[i]
    end
    return r
end

"""Coefficient-wise subtraction, no modular reduction."""
function kyber_poly_sub!(r::Vector{Int16}, a::Vector{Int16}, b::Vector{Int16})
    @inbounds for i in 1:KYBER_N
        r[i] = a[i] - b[i]
    end
    return r
end

"""Barrett-reduce every coefficient."""
function kyber_poly_reduce!(r::Vector{Int16})
    @inbounds for i in 1:KYBER_N
        r[i] = barrett_reduce(r[i])
    end
    return r
end

"""
    kyber_poly_tomont!(r)

Convert polynomial to Montgomery domain by multiplying each coefficient
by R^2 mod q = (2^32 mod q) = 1353, via fqmul.
Mirrors poly.c:poly_tomont.
"""
function kyber_poly_tomont!(r::Vector{Int16})
    f = Int16(1353)   # (1ULL << 32) % KYBER_Q
    @inbounds for i in 1:KYBER_N
        r[i] = fqmul(r[i], f)
    end
    return r
end

"""
    kyber_poly_basemul_montgomery!(r, a, b)

Full 256-coefficient pointwise multiply in NTT domain.
Uses basemul on 64 pairs of degree-1 polynomials,
with alternating +/- zetas from KYBER_ZETAS[65:128].
Mirrors poly.c:poly_basemul_montgomery.
"""
function kyber_poly_basemul_montgomery!(r::Vector{Int16},
                                        a::Vector{Int16},
                                        b::Vector{Int16})
    # C: for(i=0;i<N/4;i++) {
    #   basemul(&r[4i], &a[4i], &b[4i], zetas[64+i]);
    #   basemul(&r[4i+2], &a[4i+2], &b[4i+2], -zetas[64+i]);
    # }
    # Julia 1-indexed: zetas[64+i] in C = KYBER_ZETAS[65+i] in Julia
    for i in 0:(KYBER_N ÷ 4 - 1)   # 0:63
        z = KYBER_ZETAS[65 + i]     # C's zetas[64+i]
        base = 4i + 1               # Julia 1-indexed start
        # First pair: coeffs [base, base+1]
        kyber_basemul!(
            @view(r[base:base+1]),
            @view(a[base:base+1]),
            @view(b[base:base+1]),
            z
        )
        # Second pair: coeffs [base+2, base+3] with -zeta
        kyber_basemul!(
            @view(r[base+2:base+3]),
            @view(a[base+2:base+3]),
            @view(b[base+2:base+3]),
            -z
        )
    end
    return r
end

# ── Conditional Add Q (caddq) ──────────────────────────────────────────────
# Maps potentially negative centered representative to [0, q).
# a += (a >> 15) & q

@inline function caddq(a::Int16)::Int16
    a += (a >> 15) & KYBER_Q
    return a
end

# ── Compress / Decompress ──────────────────────────────────────────────────
# Port of pq-crystals/kyber/ref/poly.c: poly_compress / poly_decompress
# Generic d-bit versions plus the specific d=4 and d=5 optimized versions.

"""
    kyber_compress(x, d) -> UInt16

Compress: maps x in [0, q) to [0, 2^d) via round(x * 2^d / q).
Uses the same Barrett-style trick as the C reference for d=4.
"""
function kyber_compress(x::Int16, d::Int)
    # Map to positive representative
    u = caddq(x)
    # round(u * 2^d / q) = floor((u * 2^d + q/2) / q)
    return (((u % UInt32) << d) + UInt32(KYBER_Q) ÷ 2) ÷ UInt32(KYBER_Q) % UInt16 & UInt16((1 << d) - 1)
end

"""
    kyber_decompress(y, d) -> Int16

Decompress: maps y in [0, 2^d) back to round(y * q / 2^d).
"""
function kyber_decompress(y::UInt16, d::Int)::Int16
    return ((UInt32(y) * UInt32(KYBER_Q) + UInt32(1 << (d - 1))) >> d) % Int16
end

"""
    kyber_poly_compress!(r_bytes, a, d)

Compress and serialize a polynomial with d bits per coefficient.
d=4: 256 coefficients → 128 bytes
d=5: 256 coefficients → 160 bytes
Mirrors poly.c:poly_compress.
"""
function kyber_poly_compress!(r_bytes::AbstractVector{UInt8}, a::Vector{Int16}, d::Int)
    if d == 4
        # 8 coefficients → 4 bytes (4 bits each, packed in pairs)
        idx = 1
        for i in 0:(KYBER_N ÷ 8 - 1)
            t = Vector{UInt8}(undef, 8)
            for j in 0:7
                u = caddq(a[8i + j + 1])
                # Optimized Barrett-style: (u << 4) + 1665, * 80635, >> 28
                # Use `% UInt32` to match C's silent int16_t → uint32_t conversion
                # (UInt32(u) would throw InexactError on negative Int16)
                d0 = (u % UInt32) << 4
                d0 += 1665
                d0 *= 80635
                d0 >>= 28
                t[j + 1] = UInt8(d0 & 0x0f)
            end
            r_bytes[idx]     = t[1] | (t[2] << 4)
            r_bytes[idx + 1] = t[3] | (t[4] << 4)
            r_bytes[idx + 2] = t[5] | (t[6] << 4)
            r_bytes[idx + 3] = t[7] | (t[8] << 4)
            idx += 4
        end
    elseif d == 5
        # 8 coefficients → 5 bytes (5 bits each)
        idx = 1
        for i in 0:(KYBER_N ÷ 8 - 1)
            t = Vector{UInt8}(undef, 8)
            for j in 0:7
                u = caddq(a[8i + j + 1])
                # Use `% UInt32` to match C's silent int16_t → uint32_t conversion
                d0 = (u % UInt32) << 5
                d0 += 1664
                d0 *= 40318
                d0 >>= 27
                t[j + 1] = UInt8(d0 & 0x1f)
            end
            r_bytes[idx]     = t[1] | (t[2] << 5)
            r_bytes[idx + 1] = (t[2] >> 3) | (t[3] << 2) | (t[4] << 7)
            r_bytes[idx + 2] = (t[4] >> 1) | (t[5] << 4)
            r_bytes[idx + 3] = (t[5] >> 4) | (t[6] << 1) | (t[7] << 6)
            r_bytes[idx + 4] = (t[7] >> 2) | (t[8] << 3)
            idx += 5
        end
    else
        error("kyber_poly_compress!: unsupported d=$d (must be 4 or 5)")
    end
    return r_bytes
end

"""
    kyber_poly_decompress!(r, a_bytes, d)

Deserialize and decompress a polynomial from d-bit packed representation.
d=4: 128 bytes → 256 coefficients
d=5: 160 bytes → 256 coefficients
Mirrors poly.c:poly_decompress.
"""
function kyber_poly_decompress!(r::Vector{Int16}, a_bytes::AbstractVector{UInt8}, d::Int)
    if d == 4
        # 1 byte → 2 coefficients (low 4 bits, high 4 bits)
        for i in 0:(KYBER_N ÷ 2 - 1)
            r[2i + 1] = (((UInt16(a_bytes[i + 1] & 0x0f) * UInt16(KYBER_Q)) + 8) >> 4) % Int16
            r[2i + 2] = (((UInt16(a_bytes[i + 1] >> 4) * UInt16(KYBER_Q)) + 8) >> 4) % Int16
        end
    elseif d == 5
        # 5 bytes → 8 coefficients
        idx = 1
        for i in 0:(KYBER_N ÷ 8 - 1)
            t = Vector{UInt8}(undef, 8)
            t[1] = a_bytes[idx]
            t[2] = (a_bytes[idx] >> 5) | (a_bytes[idx + 1] << 3)
            t[3] = a_bytes[idx + 1] >> 2
            t[4] = (a_bytes[idx + 1] >> 7) | (a_bytes[idx + 2] << 1)
            t[5] = (a_bytes[idx + 2] >> 4) | (a_bytes[idx + 3] << 4)
            t[6] = a_bytes[idx + 3] >> 1
            t[7] = (a_bytes[idx + 3] >> 6) | (a_bytes[idx + 4] << 2)
            t[8] = a_bytes[idx + 4] >> 3
            idx += 5
            for j in 1:8
                r[8i + j] = ((UInt32(t[j] & 0x1f) * UInt32(KYBER_Q) + 16) >> 5) % Int16
            end
        end
    else
        error("kyber_poly_decompress!: unsupported d=$d (must be 4 or 5)")
    end
    return r
end

# ── Encode / Decode (12-bit Bitpacking) ────────────────────────────────────
# Port of pq-crystals/kyber/ref/poly.c: poly_tobytes / poly_frombytes

"""
    kyber_poly_tobytes!(r, a)

Serialize polynomial to bytes: 256 coefficients → 384 bytes (12 bits each).
Coefficients are mapped to positive standard representatives first.
Mirrors poly.c:poly_tobytes.
"""
function kyber_poly_tobytes!(r::AbstractVector{UInt8}, a::Vector{Int16})
    for i in 0:(KYBER_N ÷ 2 - 1)
        # Map to positive standard representatives
        # Use `% UInt16` to match C's silent int16_t → uint16_t conversion
        t0 = caddq(a[2i + 1]) % UInt16
        t1 = caddq(a[2i + 2]) % UInt16
        r[3i + 1] = (t0) % UInt8
        r[3i + 2] = ((t0 >> 8) | (t1 << 4)) % UInt8
        r[3i + 3] = (t1 >> 4) % UInt8
    end
    return r
end

"""
    kyber_poly_frombytes!(r, a)

Deserialize polynomial from bytes: 384 bytes → 256 coefficients (12 bits each).
Mirrors poly.c:poly_frombytes.
"""
function kyber_poly_frombytes!(r::Vector{Int16}, a::AbstractVector{UInt8})
    for i in 0:(KYBER_N ÷ 2 - 1)
        r[2i + 1] = ((UInt16(a[3i + 1]) | (UInt16(a[3i + 2]) << 8)) & 0x0FFF) % Int16
        r[2i + 2] = (((UInt16(a[3i + 2]) >> 4) | (UInt16(a[3i + 3]) << 4)) & 0x0FFF) % Int16
    end
    return r
end

# ── Message Encode / Decode ────────────────────────────────────────────────
# Port of pq-crystals/kyber/ref/poly.c: poly_frommsg / poly_tomsg

"""
    kyber_poly_frommsg!(r, msg)

Convert 32-byte message to polynomial.
Each bit → 0 or (q+1)/2 = 1665.
Mirrors poly.c:poly_frommsg (uses cmov_int16 for constant-time).
"""
function kyber_poly_frommsg!(r::Vector{Int16}, msg::Vector{UInt8})
    half_q = Int16((Int(KYBER_Q) + 1) ÷ 2)   # 1665
    for i in 0:(KYBER_N ÷ 8 - 1)
        for j in 0:7
            bit = (msg[i + 1] >> j) & 0x01
            # Constant-time: use mask instead of branch
            mask = -(Int16(bit))   # 0x0000 or 0xFFFF
            r[8i + j + 1] = mask & half_q
        end
    end
    return r
end

"""
    kyber_poly_tomsg!(msg, a)

Convert polynomial to 32-byte message.
Each coefficient → round(2*coeff/q) & 1.
Uses the Barrett-style trick from the C reference.
Mirrors poly.c:poly_tomsg.
"""
function kyber_poly_tomsg!(msg::Vector{UInt8}, a::Vector{Int16})
    for i in 0:(KYBER_N ÷ 8 - 1)
        msg[i + 1] = 0x00
        for j in 0:7
            # C reference does NOT call caddq here -- the Barrett-style trick
            # is designed to work on raw (possibly negative) int16_t values.
            # Use `% UInt32` to match C's silent int16_t → uint32_t assignment
            # which sign-extends then reinterprets as unsigned.
            t = a[8i + j + 1] % UInt32
            # Barrett-style: (t << 1) + 1665, * 80635, >> 28, & 1
            t <<= 1
            t += 1665
            t *= 80635
            t >>= 28
            t &= 1
            msg[i + 1] |= UInt8(t) << j
        end
    end
    return msg
end

# ── Constant-Time Utilities ────────────────────────────────────────────────
# Port of pq-crystals/kyber/ref/verify.c

"""
    kyber_verify(a, b) -> Int

Constant-time comparison of two byte arrays.
Returns 0 if equal, 1 if different.
Mirrors verify.c:verify.
"""
function kyber_verify(a::Vector{UInt8}, b::Vector{UInt8})::Int
    @assert length(a) == length(b)
    r = UInt8(0)
    @inbounds for i in 1:length(a)
        r |= a[i] ⊻ b[i]
    end
    # (-(UInt64(r))) >>> 63: maps 0→0, nonzero→1
    return Int((-UInt64(r)) >>> 63)
end

"""
    kyber_cmov!(r, x, b)

Constant-time conditional copy: if b==1, copy x into r; if b==0, do nothing.
b must be in {0, 1}. Uses XOR-mask trick for constant-time behavior.
Mirrors verify.c:cmov.
"""
function kyber_cmov!(r::Vector{UInt8}, x::Vector{UInt8}, b::UInt8)
    @assert length(r) == length(x)
    mask = -b   # 0x00 or 0xFF (two's complement wrapping for UInt8)
    @inbounds for i in 1:length(r)
        r[i] ⊻= mask & (r[i] ⊻ x[i])
    end
    return r
end

"""
    kyber_cmov_int16(r_val, v, b) -> Int16

Constant-time conditional select: returns v if b==1, r_val if b==0.
b must be in {0, 1}.
Mirrors verify.c:cmov_int16.
"""
function kyber_cmov_int16(r_val::Int16, v::Int16, b::UInt16)::Int16
    mask = (-(b % Int16)) % Int16
    return r_val ⊻ (mask & (r_val ⊻ v))
end

# ── Noise Sampling Wrappers ────────────────────────────────────────────────

"""
    kyber_poly_getnoise_eta2!(r, seed, nonce)

Sample noise polynomial with eta=2 from PRF(seed, nonce).
"""
function kyber_poly_getnoise_eta2!(r::Vector{Int16}, seed::Vector{UInt8}, nonce::UInt8)
    buf = kyber_prf(seed, nonce, 2 * KYBER_N ÷ 4)   # eta2*N/4 = 128 bytes
    kyber_cbd2!(r, buf)
    return r
end

"""
    kyber_poly_getnoise_eta1!(r, seed, nonce, eta1)

Sample noise polynomial with given eta from PRF(seed, nonce).
eta1=2 or eta1=3 depending on Kyber parameter set.
"""
function kyber_poly_getnoise_eta1!(r::Vector{Int16}, seed::Vector{UInt8},
                                   nonce::UInt8, eta1::Int)
    buf = kyber_prf(seed, nonce, eta1 * KYBER_N ÷ 4)
    if eta1 == 2
        kyber_cbd2!(r, buf)
    elseif eta1 == 3
        kyber_cbd3!(r, buf)
    else
        error("eta1 must be 2 or 3, got $eta1")
    end
    return r
end

# ══════════════════════════════════════════════════════════════════════════════
# Task 4: Kyber-768 IND-CPA encryption + CCA2 KEM wrapper
# Port of pq-crystals/kyber/ref: indcpa.c, kem.c, polyvec.c
# ══════════════════════════════════════════════════════════════════════════════
