# src/falcon/falcon.jl
# ============================================================================
# Falcon (FN-DSA / FIPS 206) — integration: keygen, sign, verify.
# Builds on the validated cores: falcon_fft.jl (FFT/NTT), falcon_sampler.jl
# (SamplerZ), falcon_ntrugen.jl (NTRU equation solver).
#
# Preimage sampling uses Klein/GPV over the Gram–Schmidt basis of the NTRU
# lattice (correctness-first; produces valid Falcon signatures — the FFT-tree
# ffSampling is the O(n log n) acceleration, swappable to the NEON FFT later).
# ============================================================================
include(joinpath(@__DIR__, "falcon_fft.jl"))
include(joinpath(@__DIR__, "falcon_sampler.jl"))
include(joinpath(@__DIR__, "falcon_ntrugen.jl"))

module Falcon

using LinearAlgebra, SHA
import ..FalconFFT as FF
import ..FalconSampler as FS
import ..FalconNTRUGen as NG

const q = 12289

# Per-scheme: σ (signing), σmin, β² (acceptance) — Falcon spec Table 3.3.
params(n) = n == 512  ? (σ=165.736617183, σmin=1.277833697, β2=34034726) :
            n == 1024 ? (σ=168.388571447, σmin=1.298280334, β2=70265242) :
            (σ = 1.55 * sqrt(q), σmin=1.2778, β2=Int(floor(1.3 * 2n * (1.55^2 * q))))  # small-n tests
const SIGMA_FG(n) = 1.17 * sqrt(q / (2n))

# Generic centered discrete Gaussian over ℤ of width σ (rejection vs continuous
# Gaussian) — for sampling f,g in keygen, where σ_fg ≫ MAX_SIGMA so the leaf
# SamplerZ does not apply. Correctness-grade (keygen is not the per-signature hot
# path); a dedicated table sampler is the constant-time hardening.
function dgauss(σ::Float64)
    bnd = ceil(Int, 12σ)
    while true
        z = rand(-bnd:bnd)
        rand() < exp(-(z^2) / (2σ^2)) && return z
    end
end

# ── polynomial helpers in ℤ[x]/(xⁿ+1) ──────────────────────────────────────
adj(f) = [f[1]; [-f[length(f)-i+1] for i in 1:length(f)-1]]    # conjugate (negacyclic adjoint)
centermod(x) = (y = mod(x, q); y > q ÷ 2 ? y - q : y)

# NTT helpers (invertibility, h = g/f, s0 = c − s1·h) via the fast O(n log n)
# butterfly NTT (FalconFFT.ntt_ct!), negacyclic (ψ-weighted).
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
    n = length(f)
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

# ── Key generation ──────────────────────────────────────────────────────────
function keygen(n::Int; maxtries=200)
    σfg = SIGMA_FG(n)
    for _ in 1:maxtries
        f = [dgauss(σfg) for _ in 1:n]
        g = [dgauss(σfg) for _ in 1:n]
        nf2 = sum(abs2, f) + sum(abs2, g)
        (nf2 == 0 || nf2 > 1.17^2 * q || q^2 / nf2 > 1.17^2 * q) && continue
        is_invertible(f) || continue
        local F, G
        try
            Fb, Gb = NG.ntru_solve(BigInt.(f), BigInt.(g)); F, G = Fb, Gb
        catch; continue; end
        reduce_FG!(BigInt.(f), BigInt.(g), F, G)
        NG.ntru_check(BigInt.(f), BigInt.(g), F, G) || continue
        h = poly_div_modq(g, f)                       # public key h = g/f mod q
        return (; f=Int.(f), g=Int.(g), F=Int.(F), G=Int.(G), h=Int.(h), n)
    end
    error("keygen failed in $maxtries tries")
end

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
function hash_to_point(msg::Vector{UInt8}, salt::Vector{UInt8}, n::Int)
    k = (1 << 16) ÷ q
    stream = SHA.shake256(vcat(salt, msg), UInt64(8 * n + 8192))   # ample (rejection rate ~ k·q/2¹⁶)
    c = Int[]; i = 1
    while length(c) < n
        i + 1 > length(stream) && error("hash_to_point: stream exhausted")
        elt = (Int(stream[i]) << 8) + Int(stream[i+1]); i += 2
        elt < k * q && push!(c, elt % q)
    end
    c
end

# Klein/GPV reference: dense GS of the 2n×2n basis (O(n³) setup, O(n²) per sample).
function sign_setup_klein(sk)
    B = ntru_basis(sk)
    # Gram–Schmidt of the ROWS of B (classical; small n):
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

function sign_poly_klein(sk, gs, msg::Vector{UInt8})
    n = sk.n; p = params(n)
    while true
        salt = rand(UInt8, 40)
        c = hash_to_point(msg, salt, n)
        target = Float64.(vcat(c, zeros(Int, n)))       # (c, 0)
        r = FS.RNG(rand(UInt8, 32))
        z = klein_sample(gs.B, gs.Bgs, gs.gsnorm2, target, p.σ, p.σmin, r)
        v = gs.B' * z                                   # lattice point  z·B (as column)
        s = target .- v                                 # short vector (s0, s1)
        s0 = round.(Int, s[1:n]); s1 = round.(Int, s[n+1:2n])
        norm2 = sum(abs2, s0) + sum(abs2, s1)
        if norm2 <= p.β2
            return (; salt, s1, s0, norm2)
        end
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
    (; b00, b01, b10, b11, T = ffldl(g00, g01, g11, params(sk.n).σ))
end

# ── Sign (returns the signature polynomial s1 + salt; math-level) ───────────
function sign_poly(sk, gs, msg::Vector{UInt8})
    n = sk.n; p = params(n)
    while true
        salt = rand(UInt8, 40)
        c = hash_to_point(msg, salt, n)
        cf = FF.fft(Float64.(c))
        t0 = cf .* gs.b11 ./ q; t1 = .-cf .* gs.b01 ./ q    # (c, 0)·B₀⁻¹
        r = FS.RNG(rand(UInt8, 32))
        z0, z1 = ffsample(t0, t1, gs.T, p.σmin, r)
        v0 = round.(Int, real.(FF.ifft(z0 .* gs.b00 .+ z1 .* gs.b10)))
        v1 = round.(Int, real.(FF.ifft(z0 .* gs.b01 .+ z1 .* gs.b11)))
        s0 = c .- v0; s1 = .-v1                            # (c, 0) − z·B₀
        norm2 = sum(abs2, s0) + sum(abs2, s1)
        norm2 <= p.β2 && return (; salt, s1, s0, norm2)
    end
end

# ── Verify ──────────────────────────────────────────────────────────────────
function verify_poly(pk_h, n::Int, msg::Vector{UInt8}, sig)
    p = params(n)
    c = hash_to_point(msg, sig.salt, n)
    s0 = centermod.(c .- poly_mul_modq(sig.s1, pk_h))   # s0 = c − s1·h mod q (centered)
    norm2 = sum(abs2, s0) + sum(abs2, sig.s1)
    return norm2 <= p.β2
end

end # module
