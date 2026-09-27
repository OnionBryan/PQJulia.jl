# src/falcon/falcon_certify.jl
# ============================================================================
# Falcon — exact key certification. The ffLDL tree of Falcon spec Alg. 8–9 is
# defined over ℚ[x]/(xⁿ+1) (tprest/falcon.py ffsampling.py: `ldl`, `ffldl` in
# coefficient representation), so its leaves d are exact rationals. Signing
# needs every normalized leaf σ/√d in [σmin, σmax] (NIST FIPS 206 status update,
# Perlner 2025); squaring gives the root-free test σ²/σmax² ≤ d ≤ σ²/σmin²
# (the same rewriting as Kaihara et al., ePrint 2026/2046). The keygen
# Gram–Schmidt bound (spec Alg. 5 line 9; ntrugen.py gs_norm) is decided exactly
# as well. No floating point is used.
#
# Field elements are (num, den): num ∈ ℤ[x]/(xⁿ+1), den ∈ ℤ>0. Inverses use the
# norm tower: a·ā(−x) = N(a)(x²) (ntrugen.py field_norm), recursively to n = 1.
# ============================================================================
module FalconCertify

import ..FalconNTRUGen as NG

const q = 12289

struct QX
    num::Vector{BigInt}
    den::BigInt
end
QX(v::AbstractVector) = QX(BigInt.(v), big(1))

function normalize(a::QX)
    g = gcd(gcd(a.num), a.den)
    g == 0 && return QX(a.num, big(1))
    s = a.den < 0 ? -g : g
    QX(a.num .÷ s, a.den ÷ s)
end
Base.:*(a::QX, b::QX) = normalize(QX(NG.negamul(a.num, b.num), a.den * b.den))
Base.:+(a::QX, b::QX) = normalize(QX(a.num .* b.den .+ b.num .* a.den, a.den * b.den))
Base.:-(a::QX, b::QX) = normalize(QX(a.num .* b.den .- b.num .* a.den, a.den * b.den))
Base.:(==)(a::QX, b::QX) = (x = normalize(a); y = normalize(b); x.num == y.num && x.den == y.den)

# adj(f)(x) = f(1/x) = f0 − f_{n−1}x − … − f1 x^{n−1}  (fft.py adj, coefficient form)
adj(a::QX) = QX([a.num[1]; -reverse(a.num[2:end])], a.den)

# Integer inverse: returns (u, d) with a·u = d ∈ ℤ, via a(x)·a(−x) = N(a)(x²).
function inv_int(a::Vector{BigInt})
    n = length(a)
    if n == 1
        a[1] == 0 && throw(DivideError())
        return [sign(a[1])], abs(a[1])
    end
    u, d = inv_int(NG.field_norm(a))
    NG.negamul(NG.galois_conjugate(a), NG.lift(u)), d
end
function Base.inv(a::QX)
    u, d = inv_int(a.num)
    normalize(QX(u .* a.den, d))
end
Base.:/(a::QX, b::QX) = a * inv(b)

# Coefficient split f = f0(x²) + x·f1(x²)  (common.py split).
split(a::QX) = (normalize(QX(a.num[1:2:end], a.den)), normalize(QX(a.num[2:2:end], a.den)))

# Leaves of ffLDL(G) for G = [[g00, g01], [adj(g01), g11]] (ffsampling.py ffldl), left to right.
function ffldl_leaves!(out, g00::QX, g01::QX, g11::QX)
    l10 = adj(g01) / g00
    d11 = g11 - l10 * adj(l10) * g00
    if length(g00.num) > 2
        a0, a1 = split(g00); b0, b1 = split(d11)
        ffldl_leaves!(out, a0, a1, a0)
        ffldl_leaves!(out, b0, b1, b0)
    else
        push!(out, g00.num[1] // g00.den, d11.num[1] // d11.den)
    end
    out
end

"""
Exact leaves by exact field arithmetic in ℚ[x]/(xⁿ+1). A reference for small n only: each
tower inverse multiplies coefficient size by the degree, so cost explodes with n.
"""
function exact_leaves_tower(f, g, F, G)
    b00, b01, b10, b11 = QX(g), QX(-f), QX(G), QX(-F)
    g00 = b00 * adj(b00) + b01 * adj(b01)
    g01 = b00 * adj(b10) + b01 * adj(b11)
    g11 = b10 * adj(b10) + b11 * adj(b11)
    ffldl_leaves!(Rational{BigInt}[], g00, g01, g11)
end

# ── Modular ffLDL + CRT (multi-modular exactness as in Pornin, ePrint 2023/290) ──
# Mod a prime p ≡ 1 (mod 2n), xⁿ+1 splits at ζ_k = ψ^{2k+1}; with ψ_{n/2} = ψ_n² the evaluation
# order, split/merge and adjoint (ζ_k⁻¹ = ζ_{n−1−k}) mirror the complex FFT in falcon_fft.jl,
# so ffLDL runs pointwise in 𝔽_p. By Ducas–Prest (ePrint 2015/1014, Cor. 1) the leaves are the
# LDL* diagonal of the re-indexed Gram matrix, d_i = Δ_i/Δ_{i−1} with Δ Gram minors of the integer
# basis, so |num|, den ≤ H = ∏‖rows‖² (Hadamard) = ‖(f,g)‖²ⁿ·‖(F,G)‖²ⁿ; primes are added until
# ∏p > 2H², where rational reconstruction is unique.

powmod_(b, e, p) = powermod(b, e, p)
function ntt_primes(n, bits)
    ps = Int[]; acc = 0.0; k = (2^31 - 1) ÷ (2n)
    while acc < bits
        p = k * 2n + 1; k -= 1
        isprime_(p) || continue
        push!(ps, p); acc += log2(p)
    end
    ps
end
function isprime_(p)
    p < 2 && return false
    for a in (2, 3, 5, 7, 11, 13, 17)                          # deterministic for p < 3.4·10¹⁴
        p == a && return true
        p % a == 0 && return false
    end
    d = p - 1; s = 0
    while iseven(d); d >>= 1; s += 1; end
    for a in (2, 3, 5, 7, 11, 13, 17)
        x = powermod(a, d, p)
        (x == 1 || x == p - 1) && continue
        composite = true
        for _ in 1:s-1
            x = mulmod(x, x, p)
            x == p - 1 && (composite = false; break)
        end
        composite && return false
    end
    true
end
mulmod(a, b, p) = Int((Int128(a) * b) % p)

function psi(n, p)                                               # primitive 2n-th root of unity
    for g in 2:p-1
        ψ = powermod(g, (p - 1) ÷ (2n), p)
        powermod(ψ, n, p) == p - 1 && return ψ
    end
    error("no 2n-th root mod $p")
end

struct ModRing
    p::Int
    roots::Dict{Int,Vector{Int}}                                 # m ↦ [ψ_m^{2k+1}], ψ_m = ψ_n^{n/m}
    inv2::Int
end
function ModRing(n, p)
    ψ = psi(n, p); r = Dict{Int,Vector{Int}}(); m = n; ψm = ψ
    while m >= 1
        r[m] = [powermod(ψm, 2k + 1, p) for k in 0:m-1]
        m ÷= 2; ψm = mulmod(ψm, ψm, p)
    end
    ModRing(p, r, invmod(2, p))
end
function nttm(R::ModRing, f::AbstractVector{Int})
    length(f) == 1 && return [mod(f[1], R.p)]
    F0 = nttm(R, f[1:2:end]); F1 = nttm(R, f[2:2:end]); mergem(R, F0, F1)
end
function mergem(R::ModRing, F0, F1)
    h = length(F0); ζ = R.roots[2h]; p = R.p; F = Vector{Int}(undef, 2h)
    for k in 1:h
        t = mulmod(ζ[k], F1[k], p)
        F[k] = mod(F0[k] + t, p); F[k+h] = mod(F0[k] - t, p)
    end
    F
end
function splitm(R::ModRing, F)
    n = length(F); h = n ÷ 2; ζ = R.roots[n]; p = R.p
    F0 = Vector{Int}(undef, h); F1 = Vector{Int}(undef, h)
    for k in 1:h
        F0[k] = mulmod(F[k] + F[k+h], R.inv2, p)
        F1[k] = mulmod(mulmod(mod(F[k] - F[k+h], p), invmod(ζ[k], p), p), R.inv2, p)
    end
    F0, F1
end
adjm(F) = reverse(F)                                             # f(ζ_k⁻¹) = F[n−1−k]

# Returns false if some pivot vanishes mod p (that prime is skipped).
function ffldl_mod!(out, R::ModRing, g00, g01, g11)
    p = R.p; n = length(g00)
    any(iszero, g00) && return false
    inv00 = [invmod(x, p) for x in g00]
    a01 = adjm(g01)
    l10 = [mulmod(a01[k], inv00[k], p) for k in 1:n]
    al10 = adjm(l10)
    d11 = [mod(g11[k] - mulmod(mulmod(l10[k], al10[k], p), g00[k], p), p) for k in 1:n]
    if n > 2
        a0, a1 = splitm(R, g00); b0, b1 = splitm(R, d11)
        ffldl_mod!(out, R, a0, a1, a0) || return false
        return ffldl_mod!(out, R, b0, b1, b0)
    end
    push!(out, mulmod(g00[1] + g00[2], R.inv2, p), mulmod(d11[1] + d11[2], R.inv2, p))
    true
end

function leaves_mod(f, g, F, G, p)
    n = length(f); R = ModRing(n, p)
    b00, b01, b10, b11 = nttm(R, Int.(g)), nttm(R, -Int.(f)), nttm(R, Int.(G)), nttm(R, -Int.(F))
    ad(v) = adjm(v)
    g00 = mod.(Int.(Int128.(b00) .* ad(b00) .% p .+ Int128.(b01) .* ad(b01) .% p), p)
    g01 = mod.(Int.(Int128.(b00) .* ad(b10) .% p .+ Int128.(b01) .* ad(b11) .% p), p)
    g11 = mod.(Int.(Int128.(b10) .* ad(b10) .% p .+ Int128.(b11) .* ad(b11) .% p), p)
    out = Int[]
    ffldl_mod!(out, R, g00, g01, g11) ? out : nothing
end

# Rational reconstruction: the unique a/b ≡ r (mod M) with |a|, b ≤ N, N = ⌊√(M/2)⌋.
function ratrecon(r::BigInt, M::BigInt, N::BigInt)
    r0, r1 = M, mod(r, M); t0, t1 = big(0), big(1)
    while r1 > N
        qq = r0 ÷ r1
        r0, r1 = r1, r0 - qq * r1
        t0, t1 = t1, t0 - qq * t1
    end
    (t1 == 0 || abs(t1) > N) && error("rational reconstruction failed")
    t1 < 0 ? (-r1 // -t1) : (r1 // t1)
end

"Exact leaves d of the ffLDL tree of B₀ = [[g, −f], [G, −F]], in signing order (multi-modular)."
function exact_leaves(f, g, F, G)
    n = length(f)
    logH = n * (log2(sum(abs2, big.(f)) + sum(abs2, big.(g))) + log2(sum(abs2, big.(F)) + sum(abs2, big.(G))))
    ps = ntt_primes(n, 2logH + 8)
    res = Vector{Vector{Int}}(); used = Int[]
    for p in ps
        v = leaves_mod(f, g, F, G, p)
        v === nothing && continue
        push!(res, v); push!(used, p)
    end
    while sum(log2, used) < 2logH + 2                             # replace skipped primes
        p = ntt_primes(n, sum(log2, ps) + 40)[end]; push!(ps, p)
        v = leaves_mod(f, g, F, G, p); v === nothing || (push!(res, v); push!(used, p))
    end
    M = prod(big.(used)); N = isqrt(M ÷ 2)
    w = [(M ÷ p) * invmod(M ÷ p, big(p)) for p in used]          # CRT idempotents
    [ratrecon(mod(sum(w[j] * res[j][i] for j in eachindex(used)), M), M, N) for i in eachindex(res[1])]
end

"Exact squared Gram–Schmidt norm max(‖(g,−f)‖², q²‖(ḡ, f̄)/(ff̄+gḡ)‖²) (ntrugen.py gs_norm)."
function exact_gs_norm2(f, g)
    F, G = QX(f), QX(g)
    ffgg = F * adj(F) + G * adj(G)
    Ft = adj(G) / ffgg; Gt = adj(F) / ffgg
    sq(a::QX) = sum(a.num .^ 2) // a.den^2
    max(big(sum(abs2, f) + sum(abs2, g)) // 1, big(q)^2 * (sq(Ft) + sq(Gt)))
end

"""
    certify(sk; σ, σmin, σmax, gs_factor=117//100) -> NamedTuple

Exact certificate for a Falcon key: every ffLDL leaf satisfies σmin ≤ σ/√d ≤ σmax, and the
Gram–Schmidt norm is ≤ gs_factor·√q. The float parameters are taken at their exact binary values.
"""
function certify(f, g, F, G; σ, σmin, σmax, gs_factor=117//100)
    d = exact_leaves(f, g, F, G)
    # Each leaf covers two Gram–Schmidt vectors and det(B₀) = qⁿ, so ∏d = qⁿ exactly (self-check).
    prod(d) == big(q)^length(f) || error("ffLDL leaf reconstruction failed the determinant check")
    s2 = Rational{BigInt}(σ)^2
    lo = s2 / Rational{BigInt}(σmax)^2; hi = s2 / Rational{BigInt}(σmin)^2
    gs2 = exact_gs_norm2(f, g)
    (; leaves_ok = all(x -> lo <= x <= hi, d), gs_ok = gs2 <= Rational{BigInt}(gs_factor)^2 * q,
       leaves = d, gs_norm2 = gs2, bounds = (lo, hi))
end

end # module
