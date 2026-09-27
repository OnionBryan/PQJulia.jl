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

using LinearAlgebra, SHA, Random
import ..FalconFFT as FF
import ..FalconSampler as FS
import ..FalconChaCha as CC
import ..FalconNTRUGen as NG
import ..FalconEncoding as FE
import ..FalconCertify as FC
import ..FalconFxp as FX

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
function ntt_fwd(a)
    n = length(a); ψ = FF.psi_2n(n); ω = powermod(ψ, 2, q)
    â = [Int(mod(a[i+1] * powermod(ψ, i, q), q)) for i in 0:n-1]
    FF.ntt_ct!(â, ω); â
end
function ntt_inv(A)
    n = length(A); ψ = FF.psi_2n(n); ψi = invmod(ψ, q); ωi = invmod(powermod(ψ, 2, q), q); ni = invmod(n, q)
    C = Int.(collect(A)); FF.ntt_ct!(C, ωi)
    [Int(mod(C[i+1] * ni % q * powermod(ψi, i, q), q)) for i in 0:n-1]
end
is_invertible(f) = all(!=(0), ntt_fwd(mod.(f, q)))
poly_div_modq(a, b) = (A = ntt_fwd(mod.(a, q)); B = ntt_fwd(mod.(b, q));
                       ntt_inv([mod(A[i] * invmod(B[i], q), q) for i in 1:length(A)]))
poly_mul_modq(a, b) = (A = ntt_fwd(mod.(a, q)); B = ntt_fwd(mod.(b, q));
                       ntt_inv([mod(A[i] * B[i], q) for i in 1:length(A)]))

# ── Babai reduce: shrink (F,G) keeping f·G − g·F = q (FFT, bitsize-scaled) ───
maxbits(v) = maximum(x -> (x == 0 ? 0 : ndigits(abs(x), base=2)), v)
function reduce_FG!(f, g, F, G)
    ff = FF.fft(Float64.(f)); gf = FF.fft(Float64.(g))
    den = real.(ff .* conj.(ff) .+ gf .* conj.(gf))
    while true
        sz = max(maxbits(F), maxbits(G)); base = max(maxbits(f), maxbits(g))
        sz <= base + 8 && break
        sh = sz - 52
        Fs = Float64.(F .>> sh); Gs = Float64.(G .>> sh)
        Ff = FF.fft(Fs); Gf = FF.fft(Gs)
        kf = (Ff .* conj.(ff) .+ Gf .* conj.(gf)) ./ den
        k = round.(Int, real.(FF.ifft(kf)))
        all(==(0), k) && break
        kb = BigInt.(k) .<< sh
        F .= F .- NG.negamul(kb, f); G .= G .- NG.negamul(kb, g)
    end
    F, G
end

# ── Key generation (spec Alg. 5 NTRUGen; tprest ntrugen.py) ─────────────────
const SIGMA_FG = 1.43300980528773                  # 1.17·√(q/8192)

# f ~ D_{ℤⁿ, σ_fg}: sum of 4096/n SamplerZ(0, 1.17√(q/8192)) draws per coefficient.
function gen_poly(n, rng)
    f0 = [FS.samplerz(0.0, SIGMA_FG, SIGMA_FG - 0.001, rng) for _ in 1:4096]
    k = 4096 ÷ n
    [sum(@view f0[(i-1)*k+1:i*k]) for i in 1:n]
end

# Squared Gram–Schmidt norm of [[g, −f], [G, −F]]: max(‖(f,g)‖², q²‖(ḡ, f̄)/(ff̄+gḡ)‖²).
function gs_norm2(f, g)
    ff = FF.fft(Float64.(f)); gf = FF.fft(Float64.(g))
    ffgg = ff .* conj.(ff) .+ gf .* conj.(gf)
    Ft = real.(FF.ifft(conj.(gf) ./ ffgg)); Gt = real.(FF.ifft(conj.(ff) ./ ffgg))
    max(sum(abs2, f) + sum(abs2, g), q^2 * (sum(abs2, Ft) + sum(abs2, Gt)))
end

fits(v, bits) = (lim = 1 << (bits - 1); all(x -> -lim < x < lim, v))

# fips206: GS bound 0.9999·1.17√q decided exactly, plus the exact leaf check (NIST FIPS 206
# status update, Perlner, Oct 2025).
# fixedpoint: no floating point (FalconFxp, ePrint 2023/290): table sampling of (f, g), fixed-point
# GS check and Babai reduction; n = 4..1024.
const FXP_GS_LIMIT_206 = Int64(round((9999 // 10000 * 117 // 100)^2 * q * big(2)^32))
function keygen(n::Int; randombytes=sysrandom, maxtries=1000, certified=false, fips206=false, fixedpoint=false)
    params(n)
    gsf = fips206 ? 0.9999 * 1.17 : 1.17
    rng = CC.ChaCha20(randombytes(SEED_LEN))
    for _ in 1:maxtries
        f, g = fixedpoint ? (FX.gauss_sample_poly(n, rng), FX.gauss_sample_poly(n, rng)) :
                            (gen_poly(n, rng), gen_poly(n, rng))
        (fits(f, FG_BITS[logn(n)]) && fits(g, FG_BITS[logn(n)])) || continue
        if fixedpoint
            FX.gs_ok(f, g; limit = fips206 ? FXP_GS_LIMIT_206 : FX.GS_LIMIT) || continue
        else
            gs_norm2(f, g) > gsf^2 * q && continue
        end
        is_invertible(f) || continue
        local F, G
        try
            F, G = (fixedpoint ? FX.ntru_solve : NG.ntru_solve)(BigInt.(f), BigInt.(g))
        catch
            continue
        end
        (fixedpoint ? FX.reduce! : reduce_FG!)(BigInt.(f), BigInt.(g), F, G)
        NG.ntru_check(BigInt.(f), BigInt.(g), F, G) || continue
        (fits(F, FFGG_BITS) && fits(G, FFGG_BITS)) || continue
        (certified || fips206) && !certify(f, g, Int.(F), Int.(G); fips206).ok && continue
        return secret_key(f, g, F, G)
    end
    error("keygen failed in $maxtries tries")
end

secret_key(f, g, F, G) = (; f=Int.(f), g=Int.(g), F=Int.(F), G=Int.(G),
                            h=poly_div_modq(g, f), n=length(f))

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
function ntru_basis(sk)
    g, f, G, F = Float64.(sk.g), Float64.(sk.f), Float64.(sk.G), Float64.(sk.F)
    [permutedims(anticirc(g)) permutedims(anticirc(-f));
     permutedims(anticirc(G)) permutedims(anticirc(-F))]
end

# ── Klein / GPV sampler over the GS basis ───────────────────────────────────
function klein_sample(B, Bgs, gsnorm2, target, σ, σmin, r)
    m = size(B, 1); c = copy(target); z = zeros(Int, m)
    for i in m:-1:1
        ci = dot(c, @view Bgs[i, :]) / gsnorm2[i]
        σi = σ / sqrt(gsnorm2[i])
        zi = FS.samplerz(ci, σi, σmin, r)
        z[i] = zi
        c .-= zi .* @view B[i, :]
    end
    z
end

# ── Hash message+salt to a point c ∈ ℤ_q^n (SHAKE256, rejection) ─────────────
function hash_to_point(msg::AbstractVector{UInt8}, salt::AbstractVector{UInt8}, n::Int)
    k = (1 << 16) ÷ q
    len = 4n + 64
    while true
        stream = SHA.shake256(vcat(salt, msg), UInt64(len))
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
    end
    (; B, Bgs, gsnorm2)
end

function sign_poly_klein(sk, gs, msg::AbstractVector{UInt8}; randombytes=sysrandom)
    n = sk.n; p = params(n)
    salt = randombytes(SALT_LEN)
    c = hash_to_point(msg, salt, n)
    target = Float64.(vcat(c, zeros(Int, n)))           # (c, 0)
    while true
        r = CC.ChaCha20(randombytes(SEED_LEN))
        z = klein_sample(gs.B, gs.Bgs, gs.gsnorm2, target, p.σ, p.σmin, r)
        s = target .- gs.B' * z                         # (c, 0) − z·B
        s0 = round.(Int, s[1:n]); s1 = round.(Int, s[n+1:2n])
        norm2 = sum(abs2, s0) + sum(abs2, s1)
        norm2 <= p.β2 && return (; salt, s1, s0, norm2)
    end
end

# ── ffSampling (Falcon spec §3.9; tprest ffsampling.py, fft_ratio = 1) ─────
# ffLDL tree of the FFT Gram matrix of B₀ = [[g, −f], [G, −F]]; leaves hold σ/√d.
struct FFLeaf
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
        return FFNode(l10, ffldl(a0, a1, a0, σ), ffldl(b0, b1, b0, σ))
    end
    FFNode(l10, FFLeaf(σ / sqrt(real(g00[1]))), FFLeaf(σ / sqrt(real(d11[1]))))
end

function ffsample(t0, t1, T, σmin, r)
    if T isa FFLeaf
        return (ComplexF64[FS.samplerz(real(t0[1]), T.σ, σmin, r)],
                ComplexF64[FS.samplerz(real(t1[1]), T.σ, σmin, r)])
    end
    z1 = FF.merge_fft(ffsample(FF.split_fft(t1)..., T.t1, σmin, r)...)
    t0b = t0 .+ (t1 .- z1) .* T.l10
    z0 = FF.merge_fft(ffsample(FF.split_fft(t0b)..., T.t0, σmin, r)...)
    (z0, z1)
end

# Per-key signing data: FFT of B₀ and its normalized ffLDL tree (reused across signatures).
function sign_setup(sk)
    b00 = FF.fft(Float64.(sk.g)); b01 = FF.fft(-Float64.(sk.f))
    b10 = FF.fft(Float64.(sk.G)); b11 = FF.fft(-Float64.(sk.F))
    g00 = b00 .* conj.(b00) .+ b01 .* conj.(b01)
    g01 = b00 .* conj.(b10) .+ b01 .* conj.(b11)
    g11 = b10 .* conj.(b10) .+ b11 .* conj.(b11)
    p = params(sk.n)
    T = ffldl(g00, g01, g11, p.σ)
    # Every leaf must lie in [σmin, σmax] or the signature distribution is wrong: refuse to sign
    # (NIST FIPS 206 status update, Perlner, Oct 2025; the round-3 spec relies on the GS-norm check alone).
    all(σ -> p.σmin <= σ <= FS.MAX_SIGMA, leaves(T)) ||
        throw(ArgumentError("Falcon key has an ffLDL leaf outside [σmin, σmax]; refusing to sign"))
    (; b00, b01, b10, b11, T)
end
leaves(T, out=Float64[]) = T isa FFLeaf ? push!(out, T.σ) : (leaves(T.t0, out); leaves(T.t1, out))

# (c, 0)·B₀⁻¹ in the FFT domain.
function target(gs, c)
    cf = FF.fft(Float64.(c))
    cf .* gs.b11 ./ q, .-cf .* gs.b01 ./ q
end

# One ffSampling attempt: (s0, s1) = (c, 0) − z·B₀, so s0 + s1·h ≡ c (mod q).
function preimage(gs, c, t0, t1, σmin, r)
    z0, z1 = ffsample(t0, t1, gs.T, σmin, r)
    v0 = round.(Int, real.(FF.ifft(z0 .* gs.b00 .+ z1 .* gs.b10)))
    v1 = round.(Int, real.(FF.ifft(z0 .* gs.b01 .+ z1 .* gs.b11)))
    c .- v0, .-v1
end

# ── Sign (spec Alg. 10): one salt, then resample until short and encodable ──
# Each attempt seeds a fresh ChaCha20 from `randombytes`, as the reference does.
function sign_poly(sk, gs, msg::AbstractVector{UInt8}; randombytes=sysrandom)
    n = sk.n; p = params(n)
    salt = randombytes(SALT_LEN)
    c = hash_to_point(msg, salt, n)
    t0, t1 = target(gs, c)
    while true
        s0, s1 = preimage(gs, c, t0, t1, p.σmin, CC.ChaCha20(randombytes(SEED_LEN)))
        norm2 = sum(abs2, s0) + sum(abs2, s1)
        norm2 <= p.β2 || continue
        enc = FE.compress(s1, (p.sigbytes - 1 - SALT_LEN) * 8)
        enc === nothing && continue
        return (; salt, s1, s0, norm2, enc)
    end
end

# ── Verify (math level) ─────────────────────────────────────────────────────
function verify_poly(pk_h, n::Int, msg::AbstractVector{UInt8}, sig)
    p = params(n)
    length(sig.s1) == n || return false
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
    out = UInt8[]; acc = UInt64(0); nacc = 0; mask = (UInt64(1) << bits) - 1
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
    vals = Vector{Int}(undef, n); acc = UInt64(0); nacc = 0; k = 0; mask = (UInt64(1) << bits) - 1
    for b in bytes
        acc = (acc << 8) | b; nacc += 8
        while nacc >= bits && k < n
            nacc -= bits; k += 1; vals[k] = Int((acc >> nacc) & mask)
        end
        acc &= (UInt64(1) << nacc) - 1
    end
    (k == n && acc == 0) ? vals : nothing
end

encode_pk(h, n) = vcat(UInt8(logn(n)), pack_bits(h, 14))

function decode_pk(pk::AbstractVector{UInt8}, n)
    (length(pk) == pk_bytes(n) && pk[1] == logn(n)) || return nothing
    h = unpack_bits(@view(pk[2:end]), n, 14)
    (h === nothing || any(>=(q), h)) ? nothing : h
end

function encode_sk(sk)
    n = sk.n; b = FG_BITS[logn(n)]
    vcat(UInt8(0x50 + logn(n)), pack_bits(sk.f, b), pack_bits(sk.g, b), pack_bits(sk.F, FFGG_BITS))
end

signext(v, bits) = v >= 1 << (bits - 1) ? v - (1 << bits) : v

# Signed fields; the value −2^(bits−1) is not a valid encoding.
function decode_signed(bytes, n, bits)
    u = unpack_bits(bytes, n, bits); u === nothing && return nothing
    s = signext.(u, bits)
    any(==(-(1 << (bits - 1))), s) ? nothing : s
end

# G is not stored: G = (q + g·F)/f mod q, exact once centered since |G| < 128.
function decode_sk(skb::AbstractVector{UInt8})
    isempty(skb) && return nothing
    lg = Int(skb[1]) - 0x50
    1 <= lg <= 10 || return nothing
    n = 1 << lg; b = FG_BITS[lg]
    length(skb) == sk_bytes(n) || return nothing
    lf = (n * b) ÷ 8
    f = decode_signed(@view(skb[2:1+lf]), n, b)
    g = decode_signed(@view(skb[2+lf:1+2lf]), n, b)
    F = decode_signed(@view(skb[2+2lf:end]), n, FFGG_BITS)
    (f === nothing || g === nothing || F === nothing) && return nothing
    is_invertible(f) || return nothing
    gF = poly_mul_modq(g, F); gF[1] = mod(gF[1] + q, q)
    G = centermod.(poly_div_modq(gF, f))
    fits(G, FFGG_BITS) || return nothing
    NG.ntru_check(BigInt.(f), BigInt.(g), BigInt.(F), BigInt.(G)) || return nothing
    secret_key(f, g, F, G)
end

# Padded signature: header 0x30+logn ‖ salt(40) ‖ compress(s1), fixed length.
encode_sig(sig, n) = vcat(UInt8(0x30 + logn(n)), sig.salt, sig.enc)

function decode_sig(sigb::AbstractVector{UInt8}, n)
    p = params(n)
    (length(sigb) == p.sigbytes && sigb[1] == 0x30 + logn(n)) || return nothing
    s1 = FE.decompress(@view(sigb[2+SALT_LEN:end]), (p.sigbytes - 1 - SALT_LEN) * 8, n)
    s1 === nothing ? nothing : (; salt=sigb[2:1+SALT_LEN], s1)
end

# ── Byte-level API ──────────────────────────────────────────────────────────
struct ExpandedKey
    sk::NamedTuple
    gs::NamedTuple
end
expand_sk(sk::NamedTuple) = ExpandedKey(sk, sign_setup(sk))
function expand_sk(skb::AbstractVector{UInt8})
    sk = decode_sk(skb); sk === nothing && throw(ArgumentError("invalid Falcon secret key"))
    expand_sk(sk)
end

function keypair(n::Int; randombytes=sysrandom, certified=false, fips206=false, fixedpoint=false)
    sk = keygen(n; randombytes, certified, fips206, fixedpoint)
    encode_pk(sk.h, n), encode_sk(sk)
end

sign(msg::AbstractVector{UInt8}, sk; randombytes=sysrandom) = sign(msg, expand_sk(sk); randombytes)
function sign(msg::AbstractVector{UInt8}, ek::ExpandedKey; randombytes=sysrandom)
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
