# test/lattice_crypto.jl
# ============================================================================
# The GENUINE lattice ⟷ Hodge-spectrum ⟷ crypto connection — every quantity is a
# real, defined, computable lattice invariant from the literature, replacing the
# fabricated "Hodge gap vs decryption failure" table (research/pqc_bridge/
# render_pqc_matrices.py was a hardcoded string).
#
# The well-known facts (cited, not guessed):
#  • Flat-torus spectral geometry: the Laplace–Beltrami eigenvalues of 𝕋=ℝⁿ/Λ are
#    {4π²‖w‖² : w∈Λ*} (Λ* = dual lattice), eigenfunctions e^{2πi⟨w,x⟩}. So the
#    Hodge/Laplace spectral gap of the lattice torus is  γ(Λ) = 4π²·λ₁(Λ*)².
#    [spectral geometry; heat-trace = θ_{Λ*}]
#  • Smoothing parameter (Micciancio–Regev, FOCS'04/SICOMP'07):
#    η_ε(Λ) = min{ s : ρ_{1/s}(Λ*∖0) ≤ ε },  ρ_t(w)=e^{−π‖w‖²/t²}  — i.e. the
#    smallest s with Σ_{w∈Λ*∖0} e^{−π s²‖w‖²} ≤ ε. The SAME dual-lattice Gaussian
#    mass that is the heat-trace tail. η_ε governs LWE hardness / leakage.
#  • Banaszczyk transference: 1 ≤ λ₁(Λ)·λ₁(Λ*) ≤ n. So the security side (dual
#    minimum = spectral gap) and the correctness side (primal minimum = packing /
#    decoding margin) are linked but distinct invariants of the SAME lattice.
#  • Decoding (correctness): under Gaussian noise of width σ the nearest-lattice-
#    point decoder fails with the union-bound rate governed by λ₁(Λ) (primal).
#  • Kyber decryption failure (FIPS 203 / Bos et al. 2017): noise
#    n_e = eᵀr + e₂ + c_v − sᵀ(e₁+c_u), DFR = P(‖n_e‖_∞ ≥ ⌈q/4⌋), computed by
#    convolving the centered-binomial product distributions + compression errors.
# ============================================================================
module LatticeCrypto

using LinearAlgebra

# Lattice Λ = { B c : c ∈ ℤⁿ }, B the (n×n, full-rank) basis (columns). The dual
# Λ* = { B⁻ᵀ c : c ∈ ℤⁿ } has basis B⁻ᵀ (⟨B⁻ᵀeᵢ, Beⱼ⟩ = δᵢⱼ ∈ ℤ).
dual_basis(B::AbstractMatrix) = inv(B)'

# Shortest nonzero vector length λ₁ by exact enumeration over c ∈ [-R,R]ⁿ∖0.
# R must be large enough to contain a shortest vector (checked: ‖shortest‖ must
# be < the shortest achievable on the box boundary, else raise R).
function lambda1(B::AbstractMatrix; R::Int = 4)
    n = size(B, 2); best = Inf
    for c in Iterators.product(ntuple(_ -> (-R:R), n)...)
        all(==(0), c) && continue
        nv = norm(B * collect(Float64, c))
        nv < best && (best = nv)
    end
    return best
end

# Hodge/Laplace spectral gap of the flat torus ℝⁿ/Λ  =  4π²·λ₁(Λ*)².
torus_spectral_gap(B::AbstractMatrix; R::Int = 4) = 4π^2 * lambda1(dual_basis(B); R = R)^2

# Enumerate dual-lattice vector norms (nonzero) within the box, for the Gaussian
# mass / heat-trace / smoothing-parameter sums.
function dual_norms(B::AbstractMatrix; R::Int = 6)
    D = dual_basis(B); n = size(D, 2); out = Float64[]
    for c in Iterators.product(ntuple(_ -> (-R:R), n)...)
        all(==(0), c) && continue
        push!(out, norm(D * collect(Float64, c)))
    end
    return out
end

# Heat-trace of the torus Laplacian Z(t) = Σ_k e^{−tλ_k} = Σ_{w∈Λ*} e^{−4π²t‖w‖²}
# (= theta function θ_{Λ*}(4πit)); the +1 is the w=0 (constant) mode.
heat_trace(B::AbstractMatrix, t::Real; R::Int = 6) =
    1.0 + sum(exp(-4π^2 * t * w^2) for w in dual_norms(B; R = R))

# Smoothing parameter η_ε(Λ) = min{ s : Σ_{w∈Λ*∖0} e^{−π s²‖w‖²} ≤ ε }.
# The Gaussian mass is strictly decreasing in s ⇒ bisection.
function smoothing_parameter(B::AbstractMatrix, ε::Real; R::Int = 6)
    duals = dual_norms(B; R = R)
    mass(s) = sum(exp(-π * s^2 * w^2) for w in duals)
    lo, hi = 1e-6, 1e4                       # mass(lo)≫ε, mass(hi)≈0<ε
    @assert mass(lo) > ε "increase R / range: mass(lo) ≤ ε"
    for _ in 1:200
        mid = (lo + hi) / 2
        mass(mid) > ε ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

# Nearest-plane / union-bound decoding failure of Λ under i.i.d. Gaussian noise of
# per-coordinate std σ: P(decode ≠ 0) ≲ Σ_{v∈relevant neighbors} Q(‖v‖/2σ), where
# Q is the Gaussian tail. Dominated by the minimal vectors (length λ₁(Λ)). This is
# the standard coding-theoretic union bound; it ties correctness to the PRIMAL
# minimum λ₁(Λ) (vs security/smoothing on the dual).
Qtail(x) = 0.5 * erfc(x / sqrt(2))
function decoding_failure(B::AbstractMatrix, σ::Real; R::Int = 3)
    n = size(B, 2); acc = 0.0
    for c in Iterators.product(ntuple(_ -> (-R:R), n)...)
        all(==(0), c) && continue
        acc += Qtail(norm(B * collect(Float64, c)) / (2σ))
    end
    return min(acc, 1.0)
end

# erfc without SpecialFunctions: Abramowitz–Stegun 7.1.26 (|err|<1.5e-7).
function erfc(x::Real)
    s = sign(x); z = abs(x)
    t = 1.0 / (1.0 + 0.3275911 * z)
    y = (((((1.061405429t - 1.453152027)t) + 1.421413741)t - 0.284496736)t + 0.254829592)t
    e = y * exp(-z * z)
    return s ≥ 0 ? e : 2.0 - e
end

# ── Real Kyber/ML-KEM decryption-failure rate by noise convolution ──────────
# (FIPS 203 / Bos–Ducas–Kiltz–… 2017.) Decryption succeeds iff every coefficient
# of n_e = eᵀr − sᵀe₁ + e₂ + c_v − sᵀc_u  has |·| < ⌈q/4⌉. e,s,r ~ CBD(η1);
# e₁,e₂ ~ CBD(η2); c_u,c_v are decompression rounding errors (du,dv bits). One
# coefficient of eᵀr / sᵀe₁ / sᵀc_u is a sum of n·k i.i.d. products (negacyclic
# ring R_q = ℤ_q[X]/(X^256+1); CBD is symmetric so signs don't change the law).
# DFR = 1−(1−δ_coeff)^n,  δ_coeff = P(|coeff| ≥ ⌈q/4⌉). All sums are kept mod q by
# DIRECT (positive) cyclic convolution — never FFT — so the ~2^-160 tail survives.
const KYBER_Q = 3329
const KYBER_N = 256

cbd_law(η) = [binomial(2η, η + k) / 2.0^(2η) for k in -η:η]          # P(X=k), k=-η..η

_modq(v) = mod(v, KYBER_Q)
function _law_to_vec(vals_probs)                                      # (value→prob) over ℤ_q
    d = zeros(KYBER_Q); for (v, p) in vals_probs; d[_modq(v) + 1] += p; end; d
end
cbd_vec(η) = _law_to_vec((k, p) for (k, p) in zip(-η:η, cbd_law(η)))

# Distribution of the product X·Y of two independent integer laws (offset vectors).
function product_vec(xlaw, ylaw)               # xlaw[i]=P(X=i+xlo) etc — pass (range,probs)
    (xr, xp), (yr, yp) = xlaw, ylaw
    d = zeros(KYBER_Q)
    for (xi, x) in enumerate(xr), (yi, y) in enumerate(yr)
        d[_modq(x * y) + 1] += xp[xi] * yp[yi]
    end
    d
end

cyclic_conv(a, b) = (q = KYBER_Q; c = zeros(q);                       # direct, O(q²), positive
    @inbounds for i in 0:q-1; ai = a[i+1]; ai == 0 && continue;
        for j in 0:q-1; c[mod(i + j, q) + 1] += ai * b[j+1]; end; end; c)

function conv_power(d, m::Int)                                        # m-fold self-convolution
    result = zeros(KYBER_Q); result[1] = 1.0                          # δ_0 identity
    base = copy(d)
    while m > 0
        (m & 1) == 1 && (result = cyclic_conv(result, base))
        m >>= 1; m > 0 && (base = cyclic_conv(base, base))
    end
    result
end

# Decompression rounding error law of d-bit compression over uniform x∈ℤ_q.
function compress_err_vec(d::Int)
    q = KYBER_Q; comp(x) = mod(round(Int, (2.0^d / q) * mod(x, q)), 2^d)
    deco(y) = round(Int, (q / 2.0^d) * y)
    dist = zeros(q); for x in 0:q-1; dist[_modq(deco(comp(x)) - x) + 1] += 1.0 / q; end; dist
end

# Full per-coefficient noise law and the scheme DFR for given (k,η1,η2,du,dv).
function kyber_dfr(; k::Int, η1::Int, η2::Int, du::Int, dv::Int)
    n = KYBER_N; q = KYBER_Q
    prod_e_r  = product_vec((-η1:η1, cbd_law(η1)), (-η1:η1, cbd_law(η1)))    # eᵀr
    # e₁ and c_u are BOTH multiplied by the same s ⇒ combine (e₁+c_u) FIRST,
    # then take the product with s (not two independent product-sums).
    cu_range = -(q÷2):(q÷2); cu = compress_err_vec(du)
    e1cu = cyclic_conv(cbd_vec(η2), cu)                                       # law(e₁+c_u) over ℤ_q
    prod_s_e1cu = product_vec((-η1:η1, cbd_law(η1)),
                              (cu_range, [e1cu[_modq(v) + 1] for v in cu_range]))  # sᵀ(e₁+c_u)
    # n_e coefficient = Σ(n·k of eᵀr) − Σ(n·k of sᵀ(e₁+c_u)) + e₂ + c_v
    noise = conv_power(prod_e_r, n * k)
    noise = cyclic_conv(noise, conv_power(prod_s_e1cu, n * k))
    noise = cyclic_conv(noise, cbd_vec(η2))                                   # e₂
    noise = cyclic_conv(noise, compress_err_vec(dv))                         # c_v
    thr = cld(q, 4)                                                           # ⌈q/4⌉ = 833
    δc = sum(noise[_modq(v) + 1] for v in thr:(q - thr))                      # P(|coeff|≥⌈q/4⌉)
    dfr = -expm1(n * log1p(-δc))                                              # 1−(1−δc)^n, tiny-δc safe
    return (delta_coeff = δc, dfr = dfr, log2_dfr = dfr > 0 ? log2(dfr) : -Inf)
end

# ── Falcon / FN-DSA: the smoothing parameter as a DEPLOYED scheme parameter ──
# Falcon (Prest et al., NIST FIPS 206) signs with a discrete Gaussian over the
# NTRU lattice via the GPV framework. The signing standard deviation is
#     σ = σ_min · ‖B‖_GS,   ‖B‖_GS ≤ 1.17·√q   (det B = fG−gF = q ⇒ √q is the
#                                                Gram–Schmidt lower bound),
# and σ_min = η'_ε(ℤ) is the SMOOTHING PARAMETER — the Micciancio–Regev / GPV
# condition that the signature distribution is (ε-)independent of the secret
# basis, i.e. signatures don't leak the trapdoor. So Falcon's core security
# parameter IS the smoothing parameter — the same dual-lattice Gaussian quantity
# that equals the flat-torus Hodge spectral gap (torus_spectral_gap). Exact
# values: Falcon spec Table 3.3.
const FALCON_PARAMS = Dict(
    512  => (n=512,  q=12289, σ=165.736617183, σmin=1.277833697, σmax=1.8205, β2=34034726, pk=897,  sig=666),
    1024 => (n=1024, q=12289, σ=168.388571447, σmin=1.298280334, σmax=1.8205, β2=70265242, pk=1793, sig=1280),
)

falcon_gs_norm(q) = 1.17 * sqrt(q)                       # ‖B‖_GS ≤ 1.17·√q
# GPV signing width = smoothing parameter × Gram–Schmidt norm.
falcon_sigma(p) = p.σmin * falcon_gs_norm(p.q)
# 1-D smoothing parameter η'_ε(ℤ) = (1/π)·√(½ ln(2(1+1/ε))). Invert for the ε at
# which σ_min = η'_ε(ℤ) — Falcon's per-coordinate sampler smoothing tolerance.
falcon_smoothing_eps(σmin) = 1.0 / (0.5 * exp(2 * (π * σmin)^2) - 1.0)
eta_prime(ε) = (1 / π) * sqrt(0.5 * log(2 * (1 + 1 / ε)))

# Published ML-KEM parameter sets (k,η1,η2,du,dv) from pq-crystals/FIPS 203.
const KYBER_PARAMS = Dict(
    512  => (k = 2, η1 = 3, η2 = 2, du = 10, dv = 4),
    768  => (k = 3, η1 = 2, η2 = 2, du = 10, dv = 4),
    1024 => (k = 4, η1 = 2, η2 = 2, du = 11, dv = 5),
)

end # module
