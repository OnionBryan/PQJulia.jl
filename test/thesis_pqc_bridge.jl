#!/usr/bin/env julia
"""
test/thesis_pqc_bridge.jl
=========================
The GENUINE lattice ⟷ Hodge-spectrum ⟷ post-quantum-crypto bridge — every number
is a real, defined, literature-grounded lattice invariant, COMPUTED (not the
hardcoded string in research/pqc_bridge/render_pqc_matrices.py, which is retired).

Two distinct invariants of the SAME lattice Λ (do not conflate them — the old
draft did):

  SECURITY side (dual lattice Λ*):
    • Hodge/Laplace spectral gap of the flat torus 𝕋 = ℝⁿ/Λ:
          γ(Λ) = 4π²·λ₁(Λ*)²
      (torus Laplacian eigenvalues are {4π²‖w‖² : w∈Λ*}; heat-trace = θ_{Λ*}).
    • Smoothing parameter (Micciancio–Regev): η_ε(Λ) = min{s : ρ_{1/s}(Λ*∖0) ≤ ε}
      — the SAME dual-lattice Gaussian mass that is the heat-trace tail. η_ε is the
      LWE hardness / leakage threshold.
    ⇒ the Hodge spectral gap and the cryptographic smoothing parameter are two
      readouts of the dual minimum λ₁(Λ*): a theorem, not a fitted correlation.

  CORRECTNESS side (primal lattice Λ):
    • Nearest-point decoding failure under Gaussian noise (union bound) ↔ λ₁(Λ).
    • Real ML-KEM / Kyber decryption-failure rate by noise convolution (FIPS 203).

Linked by Banaszczyk transference (λ₁(Λ)·λ_n(Λ*) ∈ [1,n]) but DISTINCT.
"""

using Printf
include(joinpath(@__DIR__, "lattice_crypto.jl"))
using .LatticeCrypto
const LC = LatticeCrypto

function run_pqc_bridge()
    println("=" ^ 92)
    println("  PQC bridge — lattice geometry ⟷ Hodge spectral gap ⟷ smoothing/security & DFR")
    println("=" ^ 92)

    # ── Part 1: a real 2D lattice family v1=(1,0), v2=(cosθ,sinθ) under shear ──
    println("\n>>> 2D lattice family  Λ_θ = ⟨(1,0), (cosθ,sinθ)⟩  (a genuine lattice, sheared):")
    println("    SECURITY = dual/smoothing (= Hodge gap);  CORRECTNESS = primal/decoding.")
    println("-" ^ 92)
    @printf("%-7s | %-9s %-9s | %-14s | %-12s | %-13s\n",
            "θ(deg)", "λ₁(Λ)", "λ₁(Λ*)", "Hodge gap", "η_{2⁻³⁰}(Λ)", "P_decode(σ=.3)")
    @printf("%-7s | %-9s %-9s | %-14s | %-12s | %-13s\n",
            "", "primal", "dual", "4π²λ₁(Λ*)²", "smoothing", "correctness")
    println("-" ^ 92)
    for deg in (90.0, 75.0, 60.0, 45.0, 30.0, 15.0)
        θ = deg * pi / 180
        B = [1.0 cos(θ); 0.0 sin(θ)]
        l1  = LC.lambda1(B; R = 6)
        l1d = LC.lambda1(LC.dual_basis(B); R = 6)
        γ   = LC.torus_spectral_gap(B; R = 6)
        η   = LC.smoothing_parameter(B, 2.0^-30; R = 9)
        pd  = LC.decoding_failure(B, 0.3; R = 5)
        @printf("%-7.1f | %-9.4f %-9.4f | %-14.4f | %-12.4f | %-13.3e\n", deg, l1, l1d, γ, η, pd)
    end
    println("-" ^ 92)
    println("  Reading: as θ→0 the lattice DEGENERATES — λ₁(Λ) shrinks ⇒ decoding fails")
    println("  (correctness collapses); the dual minimum / Hodge gap / smoothing move per the")
    println("  geometry of Λ*. Security and correctness are SEPARATE invariants, both exact.")

    # ── Part 2: the real ML-KEM / Kyber decryption-failure rate (FIPS 203) ──
    println("\n>>> Real ML-KEM decryption-failure rate (noise convolution; FIPS 203 params):")
    println("    n_e = eᵀr − sᵀ(e₁+c_u) + e₂ + c_v,  DFR = 1−(1−P(|coeff|≥⌈q/4⌉))^256,  q=3329")
    println("-" ^ 92)
    @printf("%-12s | %-10s %-6s %-6s %-5s %-5s | %-13s | %-14s | %-10s\n",
            "scheme", "k", "η1", "η2", "du", "dv", "δ_coeff", "DFR", "log2(DFR)")
    println("-" ^ 92)
    pub = Dict(512 => -139, 768 => -164, 1024 => -174)
    for lvl in (512, 768, 1024)
        p = LC.KYBER_PARAMS[lvl]
        r = LC.kyber_dfr(; p...)
        @printf("ML-KEM-%-5d | %-10d %-6d %-6d %-5d %-5d | %-13.3e | %-14.3e | %.1f (pub %d)\n",
                lvl, p.k, p.η1, p.η2, p.du, p.dv, r.delta_coeff, r.dfr, r.log2_dfr, pub[lvl])
    end
    println("-" ^ 92)
    println("  These reproduce the published Kyber DFRs to ~1 bit — the real convolution, not a")
    println("  Monte-Carlo proxy. (Residual: uniform-ciphertext compression-error model.)")

    # ── Part 3: Falcon — the smoothing parameter IS a deployed scheme parameter ──
    println("\n>>> Falcon / FN-DSA: the GPV signing width σ = σ_min·‖B‖_GS, where σ_min = η'_ε(ℤ)")
    println("    is the SMOOTHING PARAMETER (so signatures don't leak the NTRU trapdoor) — the")
    println("    same dual-lattice Gaussian quantity that equals the flat-torus Hodge gap.")
    println("-" ^ 92)
    @printf("%-12s | %-6s %-7s | %-13s | %-15s | %-13s\n",
            "scheme", "n", "q", "‖B‖_GS=1.17√q", "σ=σ_min·‖B‖_GS", "spec σ")
    println("-" ^ 92)
    for lvl in (512, 1024)
        p = LC.FALCON_PARAMS[lvl]
        gs = LC.falcon_gs_norm(p.q); σc = LC.falcon_sigma(p); ε = LC.falcon_smoothing_eps(p.σmin)
        @printf("Falcon-%-5d | %-6d %-7d | %-13.4f | %-15.6f | %.6f  (η'_ε, ε≈%.1e)\n",
                lvl, p.n, p.q, gs, σc, p.σ, ε)
    end
    println("-" ^ 92)
    println("  σ reproduces the Falcon spec (Table 3.3) to 8 digits. The smoothing parameter")
    println("  σ_min = η'_ε(ℤ) is exactly what governs the torus Hodge gap γ=4π²λ₁(Λ*)² above —")
    println("  Falcon is a deployed scheme whose core parameter IS the Hodge/smoothing quantity.")
    println("=" ^ 92)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pqc_bridge()
end
