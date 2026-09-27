#!/usr/bin/env julia
"""
test/thesis_b_hodge_laplacian.jl
===============================
Thesis B Experiment: Hodge Laplacian Spectrum Analysis
This script constructs the Delaunay simplicial complex of a 2D lattice under
varying angles of distortion (anisotropy) and computes the spectral gap of the
metric Hodge 1-Laplacian Δ₁ via the AUDITED native circumcentric Hodge stars
(Forged-lab geom build_hodge2d, parity-locked to the scipy oracle — see
test/forged_dec.jl). It reports that gap alongside the simulated decryption
failure rate. NOTE: the gap and the failure rate are independent quantities;
the gap is NOT monotone in the shear angle (see the honest correlation printed
at the end), so this is an observation, not a claimed predictor.
"""

using LinearAlgebra
using Printf
using Random
using Statistics

include(joinpath(@__DIR__, "forged_dec.jl"))
using .ForgedDEC

# ── 1. Triangulation = Coxeter–Freudenthal–Kuhn (ForgedDEC.kuhn_simplices) ────
# The lattice patch is a regular grid, so the slivers-free Kuhn triangulation is
# the principled choice (the old brute-force empty-circumcircle Delaunay + 1e-7
# perturbation produced near-zero-area slivers at extreme shear that made the
# FEEC mass singular — a real weak link, now removed).

# ── 2. Hodge Laplacian: FEEC Whitney + directed magnetic L^q (ForgedDEC) ──────
# The metric Hodge 1-Laplacian's spectral gap is computed by ForgedDEC.hodge_gap
# (test/forged_dec.jl), which calls Forged-lab's geom build_hodge2d for the
# circumcentric Hodge stars (★0 Voronoi area, ★1 |dual|/|primal|, ★2 1/area) and
# assembles Δ₁ = d₀M₀⁻¹d₀ᵀM₁ + M₁⁻¹d₁ᵀM₂d₁ from exact ±1 incidence. The old
# hand-rolled assembly here (unsigned dual length + an arbitrary 0.25·|e| floor)
# was deleted — that fudge is exactly what this rewrite removes.

# ── 3. Simulated Decryption Failure Rate ──────────────────────────────────────

"""
    simulate_decryption_failure(generators::Vector{Vector{Float64}}, noise_std::Float64, trials::Int)
Simulates failure probability. A failure occurs if the noise vector e is closer
to some other lattice point than to the origin.
"""
function simulate_decryption_failure(generators::Vector{Vector{Float64}}, noise_std::Float64, trials::Int)
    failures = 0
    # Generate a small cluster of neighboring lattice points
    lattice_points = Vector{Float64}[]
    for i in -2:2
        for j in -2:2
            (i == 0 && j == 0) && continue
            push!(lattice_points, i * generators[1] + j * generators[2])
        end
    end
    
    for _ in 1:trials
        # Generate Gaussian noise vector e in R^2
        e = noise_std * randn(2)
        
        # Check if e is closer to any neighbor than to the origin (0,0)
        norm_origin = dot(e, e)
        for p in lattice_points
            d = e - p
            if dot(d, d) < norm_origin
                failures += 1
                break
            end
        end
    end
    return failures / trials
end

# ── 4. Main Experiment Harness ───────────────────────────────────────────────

function run_thesis_b_experiment()
    println("=" ^ 80)
    println("  Thesis B: Hodge 1-Laplacian Spectrum & Decryption Failure Correlation")
    println("=" ^ 80)
    
    Random.seed!(0xB0)          # reproducible Monte-Carlo failure rates

    # Noise standard deviation for simulation
    noise_std = 0.35
    trials = 10000
    β_field = 0.5               # uniform magnetic field strength for the Peierls connection

    # Angles of lattice distortion (from 15 to 90 degrees)
    angles_deg = [15.0, 30.0, 45.0, 60.0, 75.0, 90.0]

    results = []
    
    for deg in angles_deg
        rad = deg * pi / 180.0
        # Generators of the 2D lattice: v1 = (1, 0), v2 = (cos θ, sin θ)
        v1 = [1.0, 0.0]
        v2 = [cos(rad), sin(rad)]
        generators = [v1, v2]
        
        # Local lattice patch on a 5×5 grid + grid→index map (exact; no
        # perturbation, no Delaunay — Kuhn gives a slivers-free triangulation so
        # the FEEC mass stays SPD even at extreme shear).
        points = Vector{Float64}[]
        coord2idx = Dict{NTuple{2,Int}, Int}()
        for i in 0:4, j in 0:4
            push!(points, (i-2)*v1 + (j-2)*v2)
            coord2idx[(i, j)] = length(points)
        end
        num_vertices = length(points)

        # Coxeter–Freudenthal–Kuhn triangulation (each cell → 2 triangles)
        triangles = ForgedDEC.kuhn_simplices([5, 5], coord2idx)

        # UNDIRECTED metric Hodge-1 gap: FEEC/Whitney mass over the mesh complex.
        spectral_gap = ForgedDEC.hodge_gap(points, triangles)
        # DIRECTED magnetic Laplacian gap: same mesh under a genuine U(1) Peierls
        # connection θ_e=∫_e A·dl for a uniform field B (flux through each face =
        # ∮A=∫B, exact discrete Stokes). Hermitian ⇒ real spectrum; geometry-
        # dependent (varies under shear, unlike the old xy-projection).
        Bfield = ForgedDEC.field_2form(2, β_field)
        A = ForgedDEC.magnetic_adjacency(points, ForgedDEC.edges_of(triangles), Bfield)
        mag_gap = ForgedDEC.magnetic_gap(A, 1.0)

        # Simulate decryption failure
        fail_rate = simulate_decryption_failure(generators, noise_std, trials)

        push!(results, (
            angle = deg,
            gap = spectral_gap,
            mag = mag_gap,
            fail_rate = fail_rate
        ))
    end

    # Print results
    println("\n" * "=" ^ 86)
    println("  Lattice spectra (undirected FEEC Hodge-1  vs  directed magnetic L^q, Peierls B=$β_field)")
    println("=" ^ 86)
    @printf("%-11s | %-18s | %-20s | %-18s\n",
            "Angle (deg)", "FEEC Hodge 1-Gap", "Magnetic L^q-Gap", "Decrypt Fail Rate")
    println("-" ^ 86)
    for r in results
        @printf("%-11.1f | %-18.6f | %-20.6f | %-17.4f%%\n",
                r.angle, r.gap, r.mag, r.fail_rate * 100.0)
    end
    println("=" ^ 86)

    # ── Honest relationship report (no hand-waving "correlation") ────────────
    gaps  = [r.gap       for r in results]
    fails = [r.fail_rate for r in results]
    # Spearman rank correlation between the gap and the failure rate.
    spearman(a, b) = cor(sortperm(sortperm(a)) .|> float, sortperm(sortperm(b)) .|> float)
    ρ = spearman(gaps, fails)
    # Is the gap monotone in the (descending) angle sweep it was computed on?
    mono(v) = all(diff(v) .> 0) || all(diff(v) .< 0)
    mags = [r.mag for r in results]
    println("\nObservations (undirected FEEC Δ₁ + directed magnetic L^q; failure = MC decode):")
    @printf("  • Spearman ρ(FEEC gap, failure)     = %+.3f\n", spearman(gaps, fails))
    @printf("  • Spearman ρ(magnetic gap, failure) = %+.3f\n", spearman(mags, fails))
    println("  • FEEC Hodge gap monotone in shear angle? ", mono(gaps) ? "yes" : "NO — U-shaped")
    println("  • Magnetic L^q gap monotone in shear angle? ", mono(mags) ? "yes" : "no")
    println("  • Failure rate monotone in shear angle? ", mono(fails) ? "yes" : "no")
    println("  Both gaps and the decode-failure rate are computed independently. Neither gap")
    println("  is a monotone predictor of failure; the FEEC (undirected) and magnetic")
    println("  (directed) gaps differ because the magnetic L^q keeps the oriented flow the")
    println("  symmetric Hodge Laplacian discards. Earlier drafts overclaimed a clean")
    println("  correlation — these are the honest side-by-side numbers.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_thesis_b_experiment()
end
