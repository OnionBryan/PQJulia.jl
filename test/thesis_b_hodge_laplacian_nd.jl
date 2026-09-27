#!/usr/bin/env julia
"""
test/thesis_b_hodge_laplacian_nd.jl
==================================
Generalizes the Hodge-Laplacian experiment to 3D and 4D distorted lattices.
Triangulates with the Coxeter–Freudenthal–Kuhn subdivision of the (sheared) grid
(ForgedDEC.kuhn_simplices) — a slivers-free, valid complex in any dimension; the
old brute-force empty-circumsphere "Delaunay" was ill-defined on the cospherical
lattice.

Both 3D and 4D report, side by side:
  • undirected metric Δ₁ gap via the FEEC / Whitney Galerkin mass (M₀,M₁,M₂ SPD,
    mesh-robust in ANY dimension — the native galerkin was extended to the k=2
    2-form mass; the diagonal circumcentric ★ that broke on right/obtuse Kuhn
    simplices is no longer used), solved as an exact symmetric-definite pencil.
  • directed magnetic Laplacian L^q gap (MagNet) — the asymmetric/chiral operator.

The decryption failure rate is an independent Monte-Carlo decode quantity, not a
function of the gap; the three are reported side by side, not "correlated".
"""

using LinearAlgebra
using Printf
using Random

include(joinpath(@__DIR__, "forged_dec.jl"))
using .ForgedDEC
# Coxeter–Freudenthal–Kuhn grid triangulation lives in ForgedDEC.kuhn_simplices.

# ── 1. Simplex-face enumeration (triangulation is Kuhn — see kuhn_simplices) ──

# Helper to generate combinations of size m from a vector a
function combinations(a, m)
    if m == 0
        return [Int[]]
    elseif m > length(a)
        return []
    elseif m == length(a)
        return [a]
    else
        res = []
        for i in 1:(length(a) - m + 1)
            for tail in combinations(a[(i+1):end], m - 1)
                push!(res, vcat(a[i], tail))
            end
        end
        return res
    end
end

# ── 2. Collect Faces and Assemble Boundary Operators ─────────────────────────

function collect_faces(d_simplices, dim)
    faces = Set{Tuple}()
    for simp in d_simplices
        for subset in combinations(collect(simp), dim+1)
            push!(faces, Tuple(sort(subset)))
        end
    end
    return sort(collect(faces))
end

function build_boundary_operator(p_simplices, p1_simplices)
    num_p = length(p_simplices)
    num_p1 = length(p1_simplices)
    
    p_map = Dict(simp => idx for (idx, simp) in enumerate(p_simplices))
    B = zeros(Float64, num_p, num_p1)
    
    for (s_idx, S) in enumerate(p1_simplices)
        S_vec = [S...]
        for j in 1:length(S)
            face_vec = vcat(S_vec[1:j-1], S_vec[j+1:end])
            face = Tuple(face_vec)
            f_idx = p_map[face]
            sign = (-1)^(j-1)
            B[f_idx, s_idx] = sign
        end
    end
    return B
end

# (The former assemble_metric_laplacian used ad-hoc 1/length² · 1/area² weights
# that are NOT a Hodge star; it was removed. 3D now uses the audited native
# circumcentric star via ForgedDEC.hodge_gap_3d, 4D uses the combinatorial gap.)

# ── 3. High-Dimensional Decryption Failure Simulation ────────────────────────

# Helper to generate all combinations of coefficients of length D
function combinations_with_replacement_nd(vals, D)
    if D == 1
        return [[v] for v in vals]
    end
    sub = combinations_with_replacement_nd(vals, D-1)
    res = Vector{Int}[]
    for s in sub
        for v in vals
            push!(res, vcat(s, [v]))
        end
    end
    return res
end

function simulate_decryption_failure_nd(generators::Vector{Vector{Float64}}, noise_std::Float64, trials::Int)
    D = length(generators)
    failures = 0
    
    # Generate neighbor lattice points with coefficients in {-1, 0, 1}
    lattice_points = Vector{Float64}[]
    coeffs = combinations_with_replacement_nd(-1:1, D)
    for c in coeffs
        if all(c .== 0)
            continue
        end
        p = zeros(D)
        for j in 1:D
            p += c[j] * generators[j]
        end
        push!(lattice_points, p)
    end
    
    for _ in 1:trials
        # Generate D-dimensional Gaussian noise vector
        e = noise_std * randn(D)
        
        # Check if e is closer to any neighbor than to the origin
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

# ── 4. Main Experiment ────────────────────────────────────────────────────────

function run_thesis_b_nd_experiment()
    Random.seed!(0xB0)
    noise_std = 0.35
    trials = 5000
    angles_deg = [15.0, 45.0, 90.0]
    β_field = 0.5               # uniform magnetic field strength for the Peierls connection

    println("=" ^ 80)
    println("  Thesis B: N-Dimensional Hodge Laplacian Spectrum & Decryption Failure")
    println("=" ^ 80)
    
    # ── 3D Case ──
    println("\n>>> Running 3D Experiments (2x2x2 grid, 27 vertices)...")
    results_3d = []
    for deg in angles_deg
        rad = deg * pi / 180.0
        # 3D Generators: v1=(1,0,0), v2=(cos θ, sin θ, 0), v3=(0,0,1)
        v1 = [1.0, 0.0, 0.0]
        v2 = [cos(rad), sin(rad), 0.0]
        v3 = [0.0, 0.0, 1.0]
        generators = [v1, v2, v3]
        
        # 3D sheared grid + its grid→index map (exact, no perturbation needed)
        points = Vector{Float64}[]
        coord2idx = Dict{NTuple{3,Int}, Int}()
        for i in 0:2, j in 0:2, k in 0:2
            push!(points, i*v1 + j*v2 + k*v3)
            coord2idx[(i, j, k)] = length(points)
        end

        # Coxeter–Freudenthal–Kuhn tetrahedralization (valid complex, no slivers)
        tets = ForgedDEC.kuhn_simplices([3, 3, 3], coord2idx)

        # Undirected metric Δ₁ gap from the FEEC Whitney mass (mesh-robust on the
        # right-angled Kuhn tets — the circumcentric ★ used to collapse here).
        spectral_gap = ForgedDEC.hodge_gap(points, tets)
        # Directed magnetic L^q gap (oriented 2-faces, xy-projected chirality).
        mag_gap = ForgedDEC.magnetic_gap(
            ForgedDEC.magnetic_adjacency(points, ForgedDEC.edges_of(tets), ForgedDEC.field_2form(3, β_field)), 1.0)

        fail_rate = simulate_decryption_failure_nd(generators, noise_std, trials)
        push!(results_3d, (angle = deg, gap = spectral_gap, mag = mag_gap, fail_rate = fail_rate))
    end

    # Print 3D Results
    println("  (undirected 3D Δ₁: FEEC Whitney metric Hodge-1  |  directed: magnetic L^q, Peierls B=$β_field)")
    println("-" ^ 82)
    @printf("%-11s | %-20s | %-20s | %-17s\n", "Angle (deg)", "Metric Hodge 1-Gap", "Magnetic L^q-Gap", "Decrypt Fail Rate")
    println("-" ^ 82)
    for r in results_3d
        @printf("%-11.1f | %-20.6f | %-20.6f | %-16.4f%%\n", r.angle, r.gap, r.mag, r.fail_rate * 100.0)
    end
    println("-" ^ 82)

    # ── 4D Case (combinatorial — no validated metric star exists for n>3) ──
    println("\n>>> Running 4D Experiments (1x1x1x1 grid, 16 vertices)...")
    results_4d = []
    for deg in angles_deg
        rad = deg * pi / 180.0
        # 4D Generators: v1=(1,0,0,0), v2=(cos θ, sin θ, 0, 0), v3=(0,0,1,0), v4=(0,0,0,1)
        v1 = [1.0, 0.0, 0.0, 0.0]
        v2 = [cos(rad), sin(rad), 0.0, 0.0]
        v3 = [0.0, 0.0, 1.0, 0.0]
        v4 = [0.0, 0.0, 0.0, 1.0]
        generators = [v1, v2, v3, v4]
        
        # 4D sheared grid (single 4-cube, 16 vertices) + grid→index map
        points = Vector{Float64}[]
        coord2idx = Dict{NTuple{4,Int}, Int}()
        for i in 0:1, j in 0:1, k in 0:1, l in 0:1
            push!(points, i*v1 + j*v2 + k*v3 + l*v4)
            coord2idx[(i, j, k, l)] = length(points)
        end
        
        # Coxeter–Freudenthal–Kuhn 4-simplices + 2-skeleton (for the magnetic flow)
        simplices_4d = ForgedDEC.kuhn_simplices([2, 2, 2, 2], coord2idx)
        triangles = collect_faces(simplices_4d, 2)

        # Undirected FEEC Whitney metric Hodge-1 gap — the Whitney 2-form mass is
        # well-defined in ANY dimension, so 4D gets a genuine metric Laplacian
        # (no circumcentric-star / "open for n>3" caveat — that limitation was the
        # diagonal star; the Galerkin mass sidesteps it).
        spectral_gap = ForgedDEC.hodge_gap(points, simplices_4d)
        # Directed magnetic L^q gap (Peierls connection; dimension-agnostic).
        mag_gap = ForgedDEC.magnetic_gap(
            ForgedDEC.magnetic_adjacency(points, ForgedDEC.edges_of(simplices_4d), ForgedDEC.field_2form(4, β_field)), 1.0)

        fail_rate = simulate_decryption_failure_nd(generators, noise_std, trials)
        push!(results_4d, (angle = deg, gap = spectral_gap, mag = mag_gap, fail_rate = fail_rate))
    end

    # Print 4D Results
    println("  (undirected 4D Δ₁: FEEC Whitney metric Hodge-1 (robust in any dim)  |")
    println("   directed: magnetic L^q, Peierls B=$β_field.)")
    println("-" ^ 82)
    @printf("%-11s | %-20s | %-20s | %-17s\n", "Angle (deg)", "FEEC Hodge 1-Gap", "Magnetic L^q-Gap", "Decrypt Fail Rate")
    println("-" ^ 82)
    for r in results_4d
        @printf("%-11.1f | %-20.6f | %-20.6f | %-16.4f%%\n", r.angle, r.gap, r.mag, r.fail_rate * 100.0)
    end
    println("=" ^ 80)
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_thesis_b_nd_experiment()
end
