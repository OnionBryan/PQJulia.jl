# src/falcon/falcon.jl
# ============================================================================
# Falcon (FN-DSA) — keygen, sign, verify, and the spec byte encodings.
# Builds on falcon_fft.jl (FFT/NTT), falcon_sampler.jl (SamplerZ),
# falcon_chacha.jl (signing PRNG), falcon_ntrugen.jl (NTRU solver) and
# falcon_encoding.jl (signature compression). Follows the Falcon round-3
# specification v1.2 and tprest/falcon.py; the signer is byte-exact with the
# round-3 C implementation given the same randomness (test/falcon_kat.jl).
# ============================================================================
module Falcon

using LinearAlgebra, Random
import ..FalconFFT as FF
import ..FalconSampler as FS
import ..FalconChaCha as CC
import ..FalconNTRUGen as NG
import ..FalconEncoding as FE
import ..FalconCertify as FC
import ..FalconFxp as FX
import ..FalconFpr as FP
import ...Wipe: wipe!
import ...Bytes: catbytes
import ...Keccak

const q = 12289
const SALT_LEN = 40
const SEED_LEN = 56

# σ, σmin, β² (sig bound), padded signature bytes — tprest/falcon.py params.
const PARAMS = Dict(
    2    => (σ=144.81253976308423, σmin=1.1165085072329104, β2=101498,   sigbytes=44),
    4    => (σ=146.83798833523608, σmin=1.1321247692325274, β2=208714,   sigbytes=47),
    8    => (σ=148.83587593064718, σmin=1.147528535373367,  β2=428865,   sigbytes=52),
    16   => (σ=151.78340713845503, σmin=1.170254078853483,  β2=892039,   sigbytes=63),
    32   => (σ=154.6747794602761,  σmin=1.1925466358390344, β2=1852696,  sigbytes=82),
    64   => (σ=157.51308555044122, σmin=1.2144300507766141, β2=3842630,  sigbytes=122),
    128  => (σ=160.30114421975344, σmin=1.235926056771981,  β2=7959734,  sigbytes=200),
    256  => (σ=163.04153322607107, σmin=1.2570545284063217, β2=16468416, sigbytes=356),
    512  => (σ=165.7366171829776,  σmin=1.2778336969128337, β2=34034726, sigbytes=666),
    1024 => (σ=168.38857144654395, σmin=1.298280334344292,  β2=70265242, sigbytes=1280))
params(n) = haskey(PARAMS, n) ? PARAMS[n] : throw(ArgumentError("Falcon degree must be 2^k, 1 ≤ k ≤ 10"))

logn(n) = trailing_zeros(n)
# Secret-key coefficient widths for f, g (by logn) and for F, G (C reference max_fg_bits/max_FG_bits).
const FG_BITS = (8, 8, 8, 8, 8, 7, 7, 6, 6, 5)
const FFGG_BITS = 8
pk_bytes(n)  = 1 + cld(14n, 8)
sk_bytes(n)  = 1 + 2 * (n * FG_BITS[logn(n)]) ÷ 8 + (n * FFGG_BITS) ÷ 8
sig_bytes(n) = params(n).sigbytes

sysrandom(k) = rand(RandomDevice(), UInt8, k)

# ── polynomial helpers in ℤ[x]/(xⁿ+1) ──────────────────────────────────────
centermod(x) = (y = mod(x, q); y > q ÷ 2 ? y - q : y)

# NTT helpers (invertibility, h = g/f, s0 = c − s1·h), negacyclic (ψ-weighted).
# Per-degree tables (public constants), built once for n = 1, 2, …, 1024: ψ^i, n⁻¹ψ⁻ⁱ, ω^j and
# ω⁻ʲ mod q with ω = ψ². The transforms are the exact integer maps of FalconFFT.ntt_ct! with
# its twiddles read from the table instead of recomputed by powermod.
struct NTTTables
    ψpow::Vector{Int}; ψinv::Vector{Int}; ωpow::Vector{Int}; ωinv::Vector{Int}
end
function NTTTables(n::Int)
    ψ = FF.psi_2n(n); ψi = invmod(ψ, q); ω = powermod(ψ, 2, q); ωi = invmod(ω, q); ni = invmod(n, q)
    NTTTables([powermod(ψ, i, q) for i in 0:n-1], [mod(ni * powermod(ψi, i, q), q) for i in 0:n-1],
              [powermod(ω, j, q) for j in 0:n-1], [powermod(ωi, j, q) for j in 0:n-1])
end
const NTT_TABLES = [NTTTables(1 << k) for k in 0:10]
ntt_tables(n) = (ispow2(n) && n <= 1024) ? NTT_TABLES[trailing_zeros(n) + 1] :
                throw(ArgumentError("NTT degree must be 2^k, 0 ≤ k ≤ 10"))

# In-place Cooley–Tukey NTT, as FalconFFT.ntt_ct!(a, ω) with W[j+1] = ω^j mod q, on entries
# in [0, q): products stay below q² < 2³², so each butterfly is one unsigned 32-bit remainder by
# the constant q and two branch-free conditional subtractions.
@inline modq_sub(x::UInt32) = x - ifelse(x >= UInt32(q), UInt32(q), UInt32(0))
function ntt_table!(a::Vector{Int}, W::Vector{Int})
    n = length(a); FF._bitrev!(a); len = 2
    @inbounds while len <= n
        h = len >> 1; step = n ÷ len
        for i in 0:len:n-1, k in 0:h-1
            w = W[k * step + 1] % UInt32
            u = a[i+k+1] % UInt32; v = (a[i+k+h+1] % UInt32) * w % UInt32(q)
            a[i+k+1] = Int(modq_sub(u + v)); a[i+k+h+1] = Int(modq_sub(u - v + UInt32(q)))
        end
        len <<= 1
    end
    a
end

function ntt_fwd(a)
    n = length(a); T = ntt_tables(n)
    â = [Int(mod(a[i] * T.ψpow[i], q)) for i in 1:n]
    ntt_table!(â, T.ωpow)
end
ntt_inv(A) = ntt_inv!([Int(mod(x, q)) for x in A])

# The helpers below run on secret polynomials (key decoding, keygen): ntt_modq reduces mod q while
# weighting, so no mod-q copy is made, the inverse transform runs in place, and every transform of
# an operand is wiped before returning. Outputs equal the ntt_fwd(mod.(a, q)) / ntt_inv forms.
"NTT of `a mod q` in a fresh array (the caller wipes it when `a` is secret)."
function ntt_modq(a)
    n = length(a); T = ntt_tables(n)
    â = Vector{Int}(undef, n)
    for (i, x) in enumerate(a)                            # any axes: a is iterated, not indexed
        @inbounds â[i] = Int(mod(Int(mod(x, q)) * T.ψpow[i], q))
    end
    ntt_table!(â, T.ωpow)
end
"Inverse NTT of `A` (entries in [0, q)) in place."
function ntt_inv!(A::Vector{Int})
    n = length(A); T = ntt_tables(n)
    ntt_table!(A, T.ωinv)
    @inbounds for i in 1:n; A[i] = Int(mod(A[i] * T.ψinv[i], q)); end
    A
end
function is_invertible(f)
    F = ntt_modq(f); z = false
    @inbounds for x in F; z |= x == 0; end
    wipe!(F)
    !z
end
# x^(q−2) mod q: fixed public exponent, constant modulus, no Euclid on secret data.
function invq(x::Int)
    r = 1
    for k in 13:-1:0
        r = mod(r * r, q)
        ((q - 2) >> k) & 1 == 1 && (r = mod(r * x, q))
    end
    r
end
# a/b and a·b in ℤ_q[x]/(xⁿ+1). The result is a fresh array; the transform of b is wiped and that
# of a becomes the result.
function poly_div_modq(a, b)
    length(a) == length(b) || throw(DimensionMismatch("polynomials of different degrees"))
    A = ntt_modq(a); B = ntt_modq(b)
    @inbounds for i in eachindex(A); A[i] = mod(A[i] * invq(B[i]), q); end
    wipe!(B)
    ntt_inv!(A)
end
function poly_mul_modq(a, b)
    length(a) == length(b) || throw(DimensionMismatch("polynomials of different degrees"))
    A = ntt_modq(a); B = ntt_modq(b)
    @inbounds for i in eachindex(A); A[i] = mod(A[i] * B[i], q); end
    wipe!(B)
    ntt_inv!(A)
end

# f·G − g·F = q over ℤ[x]/(xⁿ+1), exact in Int64 for key-sized coefficients. The j range is split
# at the wrap x^n = −1, so the inner loops have no branch (their bounds depend on indices only).
function negamul_int(a, b)
    Base.require_one_based_indexing(a, b)
    length(a) == length(b) || throw(DimensionMismatch("polynomials of different degrees"))
    n = length(a); c = zeros(Int, n)
    @inbounds for i in 1:n
        ai = a[i]
        for j in 1:n-i+1; c[i+j-1] += ai * b[j]; end
        for j in n-i+2:n; c[i+j-1-n] -= ai * b[j]; end
    end
    c
end
# Both products are wiped; the comparison with (q, 0, …, 0) accumulates without branching.
function ntru_check_int(f, g, F, G)
    d = negamul_int(f, G); e = negamul_int(g, F)
    acc = 0
    @inbounds for i in eachindex(d)
        acc |= (d[i] - e[i]) ⊻ ifelse(i == 1, q, 0)
    end
    wipe!(d, e)
    acc == 0
end

# ── Babai reduce: shrink (F,G) keeping f·G − g·F = q (FFT, bitsize-scaled) ───
maxbits(v) = maximum(x -> (x == 0 ? 0 : ndigits(abs(x), base=2)), v)
# The Float64 arrays (f, g, F, G scaled, their FFTs, the quotient and k) are wiped on every exit;
# the BigInt values cannot be (SECURITY.md).
function reduce_FG!(f, g, F, G)
    fd = Float64.(f); gd = Float64.(g)
    ff = ComplexF64[]; gf = ComplexF64[]; den = Float64[]
    Fs = Float64[]; Gs = Float64[]; Ff = ComplexF64[]; Gf = ComplexF64[]; kf = ComplexF64[]; kc = ComplexF64[]; k = Int[]
    try
        ff = FF.fft(fd); gf = FF.fft(gd)
        den = real.(ff .* conj.(ff) .+ gf .* conj.(gf))
        while true
            sz = max(maxbits(F), maxbits(G)); base = max(maxbits(f), maxbits(g))
            sz <= base + 8 && break
            sh = sz - 52
            Fs = Float64.(F .>> sh); Gs = Float64.(G .>> sh)
            Ff = FF.fft(Fs); Gf = FF.fft(Gs)
            kf = (Ff .* conj.(ff) .+ Gf .* conj.(gf)) ./ den
            kc = FF.ifft(kf)
            k = round.(Int, real.(kc))
            wipe!(Fs, Gs, Ff, Gf, kf, kc)
            all(==(0), k) && break
            kb = BigInt.(k) .<< sh
            wipe!(k)
            F .= F .- NG.negamul(kb, f); G .= G .- NG.negamul(kb, g)
        end
    finally
        wipe!(fd, gd, ff, gf, den, Fs, Gs, Ff, Gf, kf, kc, k)
    end
    F, G
end

# ── Key generation (spec Alg. 5 NTRUGen; tprest ntrugen.py) ─────────────────
const SIGMA_FG = 1.43300980528773                  # 1.17·√(q/8192)

# f ~ D_{ℤⁿ, σ_fg}: sum of 4096/n SamplerZ(0, 1.17√(q/8192)) draws per coefficient. The raw
# draws are wiped on every exit.
function gen_poly(n, rng)
    f0 = zeros(Int, 4096)
    try
        for j in 1:4096; f0[j] = FS.samplerz(0.0, SIGMA_FG, SIGMA_FG - 0.001, rng); end
        k = 4096 ÷ n
        return [sum(@view f0[(i-1)*k+1:i*k]) for i in 1:n]
    finally
        wipe!(f0)
    end
end

# Squared Gram–Schmidt norm of [[g, −f], [G, −F]]: max(‖(f,g)‖², q²‖(ḡ, f̄)/(ff̄+gḡ)‖²).
# Every Float64 array derived from f and g is wiped on every exit.
function gs_norm2(f, g)
    fd = Float64.(f); gd = Float64.(g)
    ff = ComplexF64[]; gf = ComplexF64[]; ffgg = ComplexF64[]; tf = ComplexF64[]; tg = ComplexF64[]
    cf = ComplexF64[]; cg = ComplexF64[]; Ft = Float64[]; Gt = Float64[]
    try
        ff = FF.fft(fd); gf = FF.fft(gd)
        ffgg = ff .* conj.(ff) .+ gf .* conj.(gf)
        tf = conj.(gf) ./ ffgg; tg = conj.(ff) ./ ffgg
        cf = FF.ifft(tf); cg = FF.ifft(tg)
        Ft = real.(cf); Gt = real.(cg)
        return max(sum(abs2, f) + sum(abs2, g), q^2 * (sum(abs2, Ft) + sum(abs2, Gt)))
    finally
        wipe!(fd, gd, ff, gf, ffgg, tf, tg, cf, cg, Ft, Gt)
    end
end

fits(v, bits) = (lim = 1 << (bits - 1); all(x -> -lim < x < lim, v))

# fips206: GS bound 0.9999·1.17√q decided exactly, plus the exact leaf check (NIST FIPS 206
# status update, Perlner, Oct 2025).
# fixedpoint: no floating point (FalconFxp, ePrint 2023/290): table sampling of (f, g), fixed-point
# GS check and Babai reduction; n = 4..1024.
const FXP_GS_LIMIT_206 = Int64(round((9999 // 10000 * 117 // 100)^2 * q * big(2)^32))
# The seed is wiped once it has keyed the ChaCha20 generator, and the generator (state and
# keystream) when keygen returns or throws. Each candidate's plain-data arrays are wiped in
# keygen_candidate; only the BigInt values of the NTRU solver are left (SECURITY.md).
function keygen(n::Int; randombytes::R=sysrandom, maxtries=1000, certified=false, fips206=false, fixedpoint=false) where {R}
    params(n)
    gsf = fips206 ? 0.9999 * 1.17 : 1.17
    seed = randombytes(SEED_LEN)
    rng = try CC.ChaCha20(seed) finally wipe!(seed) end
    try
        for _ in 1:maxtries
            sk = keygen_candidate(n, rng, gsf, certified, fips206, fixedpoint)
            sk === nothing || return sk
        end
    finally
        wipe!(rng)
    end
    error("keygen failed in $maxtries tries")
end

# One candidate (f, g) drawn from rng: the secret key, or `nothing` if the candidate is rejected.
# f, g and every plain-data array derived from them are wiped on every exit; the returned key holds
# copies (secret_key).
function keygen_candidate(n, rng, gsf, certified, fips206, fixedpoint)
    f = fixedpoint ? FX.gauss_sample_poly(n, rng) : gen_poly(n, rng)
    g = Int[]
    try
        g = fixedpoint ? FX.gauss_sample_poly(n, rng) : gen_poly(n, rng)
        (fits(f, FG_BITS[logn(n)]) && fits(g, FG_BITS[logn(n)])) || return nothing
        if fixedpoint
            FX.gs_ok(f, g; limit = fips206 ? FXP_GS_LIMIT_206 : FX.GS_LIMIT) || return nothing
        else
            gs_norm2(f, g) > gsf^2 * q && return nothing
        end
        is_invertible(f) || return nothing
        fb, gb = BigInt.(f), BigInt.(g)
        FG = solve_ntru(fixedpoint, fb, gb)
        FG === nothing && return nothing
        F, G = FG
        (fixedpoint ? FX.reduce! : reduce_FG!)(fb, gb, F, G)
        NG.ntru_check(fb, gb, F, G) || return nothing
        (fits(F, FFGG_BITS) && fits(G, FFGG_BITS)) || return nothing
        if certified || fips206
            Fi, Gi = Int.(F), Int.(G)
            ok = try certify(f, g, Fi, Gi; fips206).ok finally wipe!(Fi, Gi) end
            ok || return nothing
        end
        return secret_key(f, g, F, G)
    finally
        wipe!(f, g)
    end
end

# (F, G), or `nothing` when the solver reports that this (f, g) has no solution (NTRUSolveFailure).
# Any other exception, InterruptException included, propagates.
function solve_ntru(fixedpoint, f, g)
    try
        return (fixedpoint ? FX.ntru_solve : NG.ntru_solve)(f, g)
    catch e
        e isa NG.NTRUSolveFailure || rethrow()
        return nothing
    end
end

# Copies f, g, F, G into a new key; the caller keeps ownership of (and wipes) its own arrays.
function secret_key(f, g, F, G)
    h = poly_div_modq(g, f)
    (; f=Int.(f), g=Int.(g), F=Int.(F), G=Int.(G), h, n=length(f))
end

"Exact certificate (FalconCertify): every ffLDL leaf in [σmin, σmax] and GS norm ≤ 1.17√q (0.9999·1.17√q with fips206), no floating point."
function certify(f, g, F, G; fips206=false)
    p = params(length(f))
    c = FC.certify(f, g, F, G; σ=p.σ, σmin=p.σmin, σmax=FS.MAX_SIGMA,
                   gs_factor = fips206 ? 9999 // 10000 * 117 // 100 : 117 // 100)
    merge(c, (; ok = c.leaves_ok && c.gs_ok))
end
certify(sk::NamedTuple; fips206=false) = certify(sk.f, sk.g, sk.F, sk.G; fips206)

# ── Anticirculant (negacyclic) matrix of a polynomial, and the NTRU basis ───
function anticirc(p)
    n = length(p); M = zeros(Float64, n, n)
    @inbounds for j in 1:n, i in 1:n
        k = i - j                                     # row i, col j: coeff of x^{i-j}
        M[i, j] = k >= 0 ? p[k+1] : -p[k+n+1]         # x^{-1} = -x^{n-1}
    end
    M
end
# Falcon basis [[g,-f],[G,-F]] with ROWS = the lattice basis vectors x^k·(g,-f),
# x^k·(G,-F). The basis vectors are the COLUMNS of the anticirculant blocks
# (anticirc(p)·z = negamul(p,z)), so each block is transposed to put them in rows.
# The returned basis is the caller's; the Float64 copies and the blocks are wiped.
function ntru_basis(sk)
    g, mf, G, mF = Float64.(sk.g), .-Float64.(sk.f), Float64.(sk.G), .-Float64.(sk.F)
    A = (anticirc(g), anticirc(mf), anticirc(G), anticirc(mF))
    P = map(permutedims, A)
    B = [P[1] P[2];
         P[3] P[4]]
    wipe!(g, mf, G, mF); foreach(wipe!, A); foreach(wipe!, P)
    B
end

# ── Klein / GPV sampler over the GS basis ───────────────────────────────────
# Returns z (the caller wipes it); the running center c is wiped.
function klein_sample(B, Bgs, gsnorm2, target, σ, σmin, r)
    m = size(B, 1); c = copy(target); z = zeros(Int, m)
    for i in m:-1:1
        ci = dot(c, @view Bgs[i, :]) / gsnorm2[i]
        σi = σ / sqrt(gsnorm2[i])
        zi = FS.samplerz(ci, σi, σmin, r)
        z[i] = zi
        c .-= zi .* @view B[i, :]
    end
    wipe!(c)
    z
end

# ── Hash message+salt to a point c ∈ ℤ_q^n (SHAKE256, rejection) ─────────────
function hash_to_point(msg::AbstractVector{UInt8}, salt::AbstractVector{UInt8}, n::Int)
    k = (1 << 16) ÷ q
    len = 4n + 64
    inp = catbytes(salt, msg)
    while true
        stream = Keccak.shake256(inp, UInt64(len))
        c = Int[]; i = 1
        while length(c) < n && i + 1 <= length(stream)
            elt = (Int(stream[i]) << 8) + Int(stream[i+1]); i += 2
            elt < k * q && push!(c, elt % q)
        end
        length(c) == n && return c
        len *= 2                                      # SHAKE output is prefix-stable
    end
end

# Klein/GPV reference: dense GS of the 2n×2n basis (O(n³) setup, O(n²) per sample).
function sign_setup_klein(sk)
    B = ntru_basis(sk)
    m = size(B, 1); Bgs = similar(B); gsnorm2 = zeros(m)
    for i in 1:m
        v = copy(@view B[i, :])
        for j in 1:i-1
            v .-= (dot(@view(B[i, :]), @view(Bgs[j, :])) / gsnorm2[j]) .* @view Bgs[j, :]
        end
        Bgs[i, :] .= v; gsnorm2[i] = dot(v, v)
        wipe!(v)
    end
    (; B, Bgs, gsnorm2)
end

# As sign_poly: the sampler seed and generator, the lattice coordinates z and every rejected
# attempt are wiped.
function sign_poly_klein(sk, gs, msg::AbstractVector{UInt8}; randombytes::R=sysrandom) where {R}
    n = sk.n; p = params(n)
    salt = randombytes(SALT_LEN)
    c = hash_to_point(msg, salt, n)
    target = zeros(Float64, 2n); target[1:n] .= c      # (c, 0)
    while true
        seed = randombytes(SEED_LEN)
        r = try CC.ChaCha20(seed) catch; wipe!(seed); rethrow() end
        z = try
            klein_sample(gs.B, gs.Bgs, gs.gsnorm2, target, p.σ, p.σmin, r)
        finally
            wipe!(seed); wipe!(r)
        end
        Bz = gs.B' * z
        s = target .- Bz                                # (c, 0) − z·B
        wipe!(z, Bz)
        s0 = round.(Int, @view s[1:n]); s1 = round.(Int, @view s[n+1:2n])
        wipe!(s)
        norm2 = sum(abs2, s0) + sum(abs2, s1)
        norm2 <= p.β2 && return (; salt, s1, s0, norm2)
        wipe!(s0, s1)
    end
end

# ── ffSampling (Falcon spec §3.9; tprest ffsampling.py, fft_ratio = 1) ─────
# ffLDL tree of the FFT Gram matrix of B₀ = [[g, −f], [G, −F]]; leaves hold σ/√d.
mutable struct FFLeaf                              # mutable so wipe! can clear it
    σ::Float64
end
struct FFNode
    l10::Vector{ComplexF64}
    t0::Union{FFNode,FFLeaf}
    t1::Union{FFNode,FFLeaf}
end

function ffldl(g00, g01, g11, σ)
    l10 = conj.(g01) ./ g00                        # G₁₀ = adj(G₀₁)
    d11 = g11 .- l10 .* conj.(l10) .* g00
    if length(g00) > 2
        a0, a1 = FF.split_fft(g00); b0, b1 = FF.split_fft(d11)
        node = FFNode(l10, ffldl(a0, a1, a0, σ), ffldl(b0, b1, b0, σ))
        wipe!(d11, a0, a1, b0, b1)
        return node
    end
    node = FFNode(l10, FFLeaf(σ / sqrt(real(g00[1]))), FFLeaf(σ / sqrt(real(d11[1]))))
    wipe!(d11)
    node
end
wipe!(T::FFLeaf) = (T.σ = 0.0; T)
wipe!(T::FFNode) = (wipe!(T.l10); wipe!(T.t0); wipe!(T.t1); T)

# The split halves, the subtree samples and t0b are wiped (as FalconFpr.ffsample).
function ffsample(t0, t1, T, σmin, r)
    if T isa FFLeaf
        return (ComplexF64[FS.samplerz(real(t0[1]), T.σ, σmin, r)],
                ComplexF64[FS.samplerz(real(t1[1]), T.σ, σmin, r)])
    end
    a, b = FF.split_fft(t1)
    za, zb = ffsample(a, b, T.t1, σmin, r)
    z1 = FF.merge_fft(za, zb)
    t0b = t0 .+ (t1 .- z1) .* T.l10
    c, d = FF.split_fft(t0b)
    zc, zd = ffsample(c, d, T.t0, σmin, r)
    z0 = FF.merge_fft(zc, zd)
    wipe!(a, b, za, zb, t0b, c, d, zc, zd)
    (z0, z1)
end

# Per-key signing data: FFT of B₀ and its normalized ffLDL tree (reused across signatures).
# The Float64 copies of the key and the Gram matrix are wiped, and the basis and tree too when
# the key is refused.
function sign_setup(sk)
    pg, pf, pG, pF = Float64.(sk.g), .-Float64.(sk.f), Float64.(sk.G), .-Float64.(sk.F)
    b00 = FF.fft(pg); b01 = FF.fft(pf)
    b10 = FF.fft(pG); b11 = FF.fft(pF)
    wipe!(pg, pf, pG, pF)
    g00 = b00 .* conj.(b00) .+ b01 .* conj.(b01)
    g01 = b00 .* conj.(b10) .+ b01 .* conj.(b11)
    g11 = b10 .* conj.(b10) .+ b11 .* conj.(b11)
    p = params(sk.n)
    T = ffldl(g00, g01, g11, p.σ)
    wipe!(g00, g01, g11)
    # Every leaf must lie in [σmin, σmax] or the signature distribution is wrong: refuse to sign
    # (NIST FIPS 206 status update, Perlner, Oct 2025; the round-3 spec relies on the GS-norm check alone).
    L = leaves(T, sizehint!(Float64[], sk.n))
    ok = all(σ -> p.σmin <= σ <= FS.MAX_SIGMA, L)
    wipe!(L)
    if !ok
        wipe!(T); wipe!(b00, b01, b10, b11)
        throw(ArgumentError("Falcon key has an ffLDL leaf outside [σmin, σmax]; refusing to sign"))
    end
    (; b00, b01, b10, b11, T)
end
leaves(T, out=Float64[]) = T isa FFLeaf ? push!(out, T.σ) : (leaves(T.t0, out); leaves(T.t1, out))

# (c, 0)·B₀⁻¹ in the FFT domain.
function target(gs::NamedTuple, c)
    cf = FF.fft(Float64.(c))
    cf .* gs.b11 ./ q, .-cf .* gs.b01 ./ q
end
target(gs::FP.Setup, c) = FP.target(gs, c)

# One ffSampling attempt: (s0, s1) = (c, 0) − z·B₀, so s0 + s1·h ≡ c (mod q). The samples z,
# z·B₀ and its inverse FFT are wiped (as FalconFpr.preimage).
preimage(gs::FP.Setup, c, t0, t1, σmin, r) = FP.preimage(gs, c, t0, t1, σmin, r)
function preimage(gs::NamedTuple, c, t0, t1, σmin, r)
    z0, z1 = ffsample(t0, t1, gs.T, σmin, r)
    u0 = z0 .* gs.b00 .+ z1 .* gs.b10; u1 = z0 .* gs.b01 .+ z1 .* gs.b11
    w0 = FF.ifft(u0); w1 = FF.ifft(u1)
    v0 = round.(Int, real.(w0)); v1 = round.(Int, real.(w1))
    s0, s1 = c .- v0, .-v1
    wipe!(z0, z1, u0, u1, w0, w1, v0, v1)
    s0, s1
end

# Leaf range of the integer signer's tree, re-checked before every signature (FalconFpr.check_leaves);
# the Float64 reference signer checks once, in sign_setup.
check_leaves(gs::FP.Setup, σmin) = FP.check_leaves(gs, σmin)
check_leaves(gs::NamedTuple, σmin) = gs

# ── Sign (spec Alg. 10): one salt, then resample until short and encodable ──
# Each attempt seeds a fresh ChaCha20 from `randombytes`, as the reference does. The target
# (t0, t1), every attempt's seed and generator, and each rejected (s0, s1) are wiped, also when
# signing throws.
function sign_poly(sk, gs, msg::AbstractVector{UInt8}; randombytes::R=sysrandom) where {R}
    n = sk.n; p = params(n)
    check_leaves(gs, p.σmin)
    salt = randombytes(SALT_LEN)
    c = hash_to_point(msg, salt, n)
    t0, t1 = target(gs, c)
    try
        while true
            s0, s1 = sample_preimage(gs, c, t0, t1, p.σmin, randombytes)
            norm2 = sum(abs2, s0) + sum(abs2, s1)
            if norm2 <= p.β2
                enc = FE.compress(s1, (p.sigbytes - 1 - SALT_LEN) * 8)
                enc === nothing || return (; salt, s1, s0, norm2, enc)
            end
            wipe!(s0, s1)
        end
    finally
        wipe!(t0, t1)
    end
end

"""
    sample_preimage(gs, c, t0, t1, σmin, randombytes) -> (s0, s1)

One ffSampling attempt with SamplerZ fed by a ChaCha20 generator seeded from
`randombytes(SEED_LEN)`; the seed and the generator (state and keystream) are wiped on every exit.
"""
function sample_preimage(gs, c, t0, t1, σmin, randombytes::R) where {R}
    seed = randombytes(SEED_LEN)
    r = try CC.ChaCha20(seed) catch; wipe!(seed); rethrow() end
    try
        return preimage(gs, c, t0, t1, σmin, r)
    finally
        wipe!(seed); wipe!(r)
    end
end

# ── Verify (math level) ─────────────────────────────────────────────────────
function verify_poly(pk_h, n::Int, msg::AbstractVector{UInt8}, sig)
    p = params(n)
    (length(sig.s1) == n && length(pk_h) == n) || return false
    c = hash_to_point(msg, sig.salt, n)
    s0 = centermod.(c .- poly_mul_modq(sig.s1, pk_h))   # s0 = c − s1·h mod q (centered)
    norm2 = sum(abs2, s0) + sum(abs2, sig.s1)
    return norm2 <= p.β2
end

# Public key in NTT form, ĥ = NTT(f)⁻¹·NTT(g) (FIPS 206 default), and verification from it.
pk_ntt(h) = ntt_fwd(mod.(h, q))
function verify_poly_ntt(ĥ, n::Int, msg::AbstractVector{UInt8}, sig)
    p = params(n)
    (length(sig.s1) == n && length(ĥ) == n) || return false
    c = hash_to_point(msg, sig.salt, n)
    S = ntt_fwd(mod.(sig.s1, q))
    s0 = centermod.(c .- ntt_inv([mod(S[i] * ĥ[i], q) for i in 1:n]))
    sum(abs2, s0) + sum(abs2, sig.s1) <= p.β2
end

# ── Byte encodings (spec §3.11; C reference codec.c) ────────────────────────
# Fixed-width fields, MSB-first; the final partial byte is zero-padded on the right.
function pack_bits(vals, bits)
    out = sizehint!(UInt8[], cld(length(vals) * bits, 8))   # no regrowth copies
    acc = UInt64(0); nacc = 0; mask = (UInt64(1) << bits) - 1
    for v in vals
        acc = (acc << bits) | ((Int64(v) % UInt64) & mask); nacc += bits
        while nacc >= 8
            nacc -= 8; push!(out, UInt8((acc >> nacc) & 0xff))
        end
        acc &= (UInt64(1) << nacc) - 1
    end
    nacc > 0 && push!(out, UInt8((acc << (8 - nacc)) & 0xff))
    out
end
# Returns the unsigned fields, or `nothing` on a length mismatch or non-zero padding.
function unpack_bits(bytes, n, bits)
    length(bytes) == cld(n * bits, 8) || return nothing
    vals = zeros(Int, n); acc = UInt64(0); nacc = 0; k = 0; mask = (UInt64(1) << bits) - 1
    try                                                     # a read of `bytes` can throw
        for b in bytes
            acc = (acc << 8) | b; nacc += 8
            while nacc >= bits && k < n
                nacc -= bits; k += 1; vals[k] = Int((acc >> nacc) & mask)
            end
            acc &= (UInt64(1) << nacc) - 1
        end
    catch
        wipe!(vals); rethrow()
    end
    (k == n && acc == 0) ? vals : (wipe!(vals); nothing)   # secret-key fields: no partial copy left
end

encode_pk(h, n) = catbytes(UInt8(logn(n)), pack_bits(h, 14))

function decode_pk(pk::AbstractVector{UInt8}, n)
    (length(pk) == pk_bytes(n) && pk[1] == logn(n)) || return nothing
    h = unpack_bits(@view(pk[2:end]), n, 14)
    (h === nothing || any(>=(q), h)) ? nothing : h
end

function encode_sk(sk)
    n = sk.n; b = FG_BITS[logn(n)]
    pf, pg, pF = pack_bits(sk.f, b), pack_bits(sk.g, b), pack_bits(sk.F, FFGG_BITS)
    out = catbytes(UInt8(0x50 + logn(n)), pf, pg, pF)
    wipe!(pf, pg, pF)
    out
end

signext(v, bits) = v >= 1 << (bits - 1) ? v - (1 << bits) : v

# Signed fields; the value −2^(bits−1) is not a valid encoding.
function decode_signed(bytes, n, bits)
    u = unpack_bits(bytes, n, bits); u === nothing && return nothing
    s = signext.(u, bits)
    wipe!(u)
    any(==(-(1 << (bits - 1))), s) ? (wipe!(s); nothing) : s
end

# G is not stored: G = (q + g·F)/f mod q, exact once centered since |G| < 128. The decoded
# f, g, F, G, the product g·F and every transform the helpers make are wiped on every exit; G is
# centered in place, so no uncentered copy (G mod q) is left.
function decode_sk(skb::AbstractVector{UInt8})
    isempty(skb) && return nothing
    lg = Int(skb[1]) - 0x50
    1 <= lg <= 10 || return nothing
    n = 1 << lg; b = FG_BITS[lg]
    length(skb) == sk_bytes(n) || return nothing
    lf = (n * b) ÷ 8
    f = g = F = nothing
    try                                                     # a read of `skb` can throw
        f = decode_signed(@view(skb[2:1+lf]), n, b)
        g = decode_signed(@view(skb[2+lf:1+2lf]), n, b)
        F = decode_signed(@view(skb[2+2lf:end]), n, FFGG_BITS)
    catch
        for x in (f, g, F); x === nothing || wipe!(x); end
        rethrow()
    end
    if f === nothing || g === nothing || F === nothing
        for x in (f, g, F); x === nothing || wipe!(x); end
        return nothing
    end
    sk = nothing; gF = Int[]; G = Int[]
    try
        if is_invertible(f)
            gF = poly_mul_modq(g, F); gF[1] = mod(gF[1] + q, q)
            G = poly_div_modq(gF, f)
            G .= centermod.(G)
            if fits(G, FFGG_BITS) && ntru_check_int(f, g, F, G)
                sk = secret_key(f, g, F, G)
            end
        end
    finally
        wipe!(f, g, F, G, gF)
    end
    sk
end

# Padded signature: header 0x30+logn ‖ salt(40) ‖ compress(s1), fixed length.
encode_sig(sig, n) = catbytes(UInt8(0x30 + logn(n)), sig.salt, sig.enc)

function decode_sig(sigb::AbstractVector{UInt8}, n)
    p = params(n)
    (length(sigb) == p.sigbytes && sigb[1] == 0x30 + logn(n)) || return nothing
    s1 = FE.decompress(@view(sigb[2+SALT_LEN:end]), (p.sigbytes - 1 - SALT_LEN) * 8, n)
    s1 === nothing ? nothing : (; salt=sigb[2:1+SALT_LEN], s1)
end

# ── Byte-level API ──────────────────────────────────────────────────────────
struct ExpandedKey{S<:NamedTuple}
    sk::S
    gs::FP.Setup                                   # integer-emulated signer (falcon_fpr.jl)
end
# Secret material: printing shows the degree only.
Base.show(io::IO, ek::ExpandedKey) = print(io, "Falcon.ExpandedKey(n = ", ek.sk.n, ", <secret>)")
Base.show(io::IO, ::FFLeaf) = print(io, "Falcon.FFLeaf(<secret>)")
Base.show(io::IO, T::FFNode) = print(io, "Falcon.FFNode(n = ", length(T.l10), ", <secret>)")
expand_sk(sk::NamedTuple) = ExpandedKey(sk, FP.setup(sk, params(sk.n).σ, params(sk.n).σmin))
# The decoded key belongs to the ExpandedKey; if expansion refuses the key (a leaf outside
# [σmin, σmax], or the fault check), it is wiped before the error propagates.
function expand_sk(skb::AbstractVector{UInt8})
    sk = decode_sk(skb); sk === nothing && throw(ArgumentError("invalid Falcon secret key"))
    try
        return expand_sk(sk)
    catch
        wipe_sk!(sk)
        rethrow()
    end
end
wipe_sk!(sk) = (wipe!(sk.f, sk.g, sk.F, sk.G); sk)
function wipe!(ek::ExpandedKey)
    wipe_sk!(ek.sk)
    wipe!(ek.gs)
    ek
end

function keypair(n::Int; randombytes::R=sysrandom, certified=false, fips206=false, fixedpoint=false) where {R}
    sk = keygen(n; randombytes, certified, fips206, fixedpoint)
    try
        return encode_pk(sk.h, n), encode_sk(sk)
    finally
        wipe_sk!(sk)
    end
end

# The signing data built for one call (FFT basis and ffLDL tree) is wiped on every exit; the
# caller's key is left as it is.
function sign(msg::AbstractVector{UInt8}, sk::NamedTuple; randombytes::R=sysrandom) where {R}
    ek = expand_sk(sk)
    try
        return sign(msg, ek; randombytes)
    finally
        wipe!(ek.gs)
    end
end
# The key decoded and expanded for one call is wiped on every exit, also when signing throws.
function sign(msg::AbstractVector{UInt8}, skb::AbstractVector{UInt8}; randombytes::R=sysrandom) where {R}
    ek = expand_sk(skb)
    try
        return sign(msg, ek; randombytes)
    finally
        wipe!(ek)
    end
end
function sign(msg::AbstractVector{UInt8}, ek::ExpandedKey; randombytes::R=sysrandom) where {R}
    encode_sig(sign_poly(ek.sk, ek.gs, msg; randombytes), ek.sk.n)
end

function verify(msg::AbstractVector{UInt8}, sigb::AbstractVector{UInt8}, pkb::AbstractVector{UInt8})
    isempty(pkb) && return false
    lg = Int(pkb[1]); 1 <= lg <= 10 || return false
    n = 1 << lg
    h = decode_pk(pkb, n); h === nothing && return false
    sig = decode_sig(sigb, n); sig === nothing && return false
    verify_poly(h, n, msg, sig)
end

end # module
