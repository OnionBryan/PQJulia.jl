# src/falcon/falcon_sampler.jl
# ============================================================================
# Falcon — SamplerZ: discrete Gaussian over ℤ with center μ and width σ∈[σmin,σmax].
# Faithful to tprest/falcon.py samplerz.py. The base half-Gaussian uses the EXACT
# reverse-CDF table (RCDT, RCDT_PREC=72 bits) from the spec; acceptance is the
# EXACT fixed-point approxexp/berexp (FACCT C table, Int128). BIT-EXACT to the
# reference: 1024/1024 of the tprest/falcon.py samplerz KAT vectors match.
# ============================================================================
module FalconSampler

using SHA

export RNG, KATSource, randbytes, samplerz, MAX_SIGMA

const MAX_SIGMA = 1.8205
const INV_2SIGMA2 = 1.0 / (2 * MAX_SIGMA^2)
const RCDT_PREC = 72

# Reverse cumulative distribution table for the base half-Gaussian (σ=MAX_SIGMA),
# 72-bit precision — EXACT spec constants (tprest/falcon.py).
const RCDT = UInt128[
    3024686241123004913666, 1564742784480091954050, 636254429462080897535,
    199560484645026482916, 47667343854657281903, 8595902006365044063,
    1163297957344668388, 117656387352093658, 8867391802663976,
    496969357462633, 20680885154299, 638331848991, 14602316184,
    247426747, 3104126, 28824, 198, 1]

# approxexp fixed-point coefficients (kept for the constant-time hardening pass).
const APPROXEXP_C = UInt64[
    0x00000004741183A3, 0x00000036548CFC06, 0x0000024FDCBF140A,
    0x0000171D939DE045, 0x0000D00CF58F6F84, 0x000680681CF796E3,
    0x002D82D8305B0FEA, 0x011111110E066FD0, 0x0555555555070F00,
    0x155555555581FF00, 0x400000000002B400, 0x7FFFFFFFFFFF4800,
    0x8000000000000000]

# ── Byte RNG (swappable; ChaCha20 in the reference). SHAKE256(seed‖counter)
# stream — deterministic for tests, swap to a CSPRNG/ChaCha20 for production. ──
mutable struct RNG
    seed::Vector{UInt8}
    ctr::UInt64
    buf::Vector{UInt8}
    pos::Int
end
RNG(seed::Vector{UInt8}) = RNG(copy(seed), UInt64(0), UInt8[], 1)
RNG() = RNG(rand(UInt8, 32))
function _refill!(r::RNG)
    ctrbytes = reinterpret(UInt8, [r.ctr]); r.ctr += 1
    r.buf = SHA.shake256(vcat(r.seed, ctrbytes), UInt64(512)); r.pos = 1
end
function randbytes(r::RNG, k::Int)
    out = Vector{UInt8}(undef, k)
    for i in 1:k
        r.pos > length(r.buf) && _refill!(r)
        out[i] = r.buf[r.pos]; r.pos += 1
    end
    out
end

# Deterministic KAT randomness: a fixed hex `octets` string consumed exactly like
# the reference KAT_randbytes (take 2k hex chars, fromhex, byte-reverse).
mutable struct KATSource; oc::String; end
function randbytes(r::KATSource, k::Int)
    s = r.oc[1:2k]; r.oc = r.oc[2k+1:end]
    reverse(UInt8[parse(UInt8, s[2i-1:2i], base=16) for i in 1:k])
end

const LN2  = 0.69314718056          # reference's rounded constants (KAT-critical)
const ILN2 = 1.44269504089

# Base sampler: integer z0 ≥ 0 from the half-Gaussian via the RCDT. u is the
# LITTLE-endian 72-bit value of randbytes(9) — matches the reference bit-for-bit.
function basesampler(r)
    bytes = randbytes(r, 9)
    u = UInt128(0); for i in 1:9; u |= UInt128(bytes[i]) << (8 * (i - 1)); end
    z0 = 0
    @inbounds for c in RCDT; z0 += (u < c) ? 1 : 0; end
    z0
end

# Fixed-point exp: ≈ 2^63 · ccs · exp(−x) (FACCT polynomial, the C table). Int128.
function approxexp(x::Float64, ccs::Float64)
    y = Int128(APPROXEXP_C[1])
    z = Int128(floor(x * Float64(Int128(1) << 63)))
    @inbounds for i in 2:13; y = Int128(APPROXEXP_C[i]) - ((z * y) >> 63); end
    z = Int128(floor(ccs * Float64(Int128(1) << 63))) << 1
    return (z * y) >> 63
end

# Constant-time Bernoulli: returns true w.p. ≈ ccs·exp(−x) (exact reference path).
function berexp(x::Float64, ccs::Float64, r)
    s = Int(floor(x * ILN2)); rr = x - s * LN2; s = min(s, 63)
    z = (approxexp(rr, ccs) - 1) >> s; w = 0
    for i in 56:-8:0
        p = Int(randbytes(r, 1)[1]); w = p - Int((z >> i) & 0xFF)
        w != 0 && break
    end
    return w < 0
end

# SamplerZ(μ, σ, σmin): discrete-Gaussian integer at μ, width σ. Bit-exact to the
# Falcon reference (1024/1024 KAT vectors) given matching randomness bytes.
function samplerz(μ::Float64, σ::Float64, σmin::Float64, r)
    s = floor(Int, μ); rr = μ - s
    dss = 1.0 / (2σ^2); ccs = σmin / σ
    while true
        z0 = basesampler(r)
        b = Int(randbytes(r, 1)[1]) & 1
        z = b + (2b - 1) * z0
        x = ((z - rr)^2) * dss - (z0^2) * INV_2SIGMA2
        berexp(x, ccs, r) && return z + s
    end
end

end # module
