# test/forged_dec.jl
# ============================================================================
# Julia FFI bridge to libforgeddec.dylib — Forged-lab's audited native DEC.
#
# Metric Hodge-1 Laplacian, ANY dimension, the FEEC / WHITNEY route — the one
# mesh-robust discretization. The native shim builds a linarm::SimplicialComplex
# over the actual mesh (signed boundaries B₁,B₂, ∂∂=0 verified) and the Whitney
# Galerkin mass matrices M₀,M₁,M₂ (Arnold–Falk–Winther; all SPD, no
# well-centeredness assumption — the diagonal circumcentric ★ breaks on
# right/obtuse simplices, the Galerkin ★ does not). The Hodge-1 Laplacian
#   Δ₁ = δ₂d₁ + d₀δ₁,  δ_k = M_{k-1}⁻¹ B_k M_k
# is assembled in the M^{1/2} similarity (M₀,M₁,M₂ SPD ⇒ clean square roots, no
# guards). Works uniformly for triangles (2D), tets (3D), pentatopes (4D).
#
# Plus the directed magnetic Laplacian L^q (MagNet) for asymmetric graphs and the
# native SU(2) Wilson-loop holonomy. Build: native/build_shim.sh.
# ============================================================================
module ForgedDEC

using LinearAlgebra

const LIB = abspath(joinpath(@__DIR__, "..", "native", "libforgeddec.dylib"))
isavailable() = isfile(LIB)

# ── FEEC Whitney building blocks over an arbitrary-dimension mesh ────────────
# points: vector-of-vectors (dim ≥ 2). simplices: top simplices as index tuples,
# 1-based (sv = dim+1 vertices each: 3 triangles, 4 tets, 5 pentatopes). Returns
# a NamedTuple in the complex's DOF order: signed boundaries B1 (nV×nE), B2
# (nE×nF); SPD Whitney masses M0,M1,M2; dd2 (=0, ∂₁∂₂ correctness invariant).
function feec_blocks(points::Vector{Vector{Float64}}, simplices)
    nV = length(points); dim = length(points[1])
    nS = length(simplices); sv = length(first(simplices))
    coords = Vector{Float64}(undef, nV * dim)
    @inbounds for i in 1:nV, d in 1:dim; coords[(i-1)*dim + d] = points[i][d]; end
    sflat = Vector{Int32}(undef, nS * sv)
    @inbounds for (i, t) in enumerate(simplices)
        s = sort(collect(t)); for j in 1:sv; sflat[(i-1)*sv + j] = Int32(s[j] - 1); end
    end
    rE = Ref{Cint}(0); rF = Ref{Cint}(0); rdd = Ref{Clonglong}(0)
    ccall((:forged_feec_hodge1, LIB), Cint,
          (Ptr{Float64}, Cint, Cint, Ptr{Cint}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Clonglong},
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          coords, nV, dim, sflat, nS, sv, rE, rF, rdd, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL) == 0 ||
        error("forged_feec_hodge1 query failed")
    nE = Int(rE[]); nF = Int(rF[])
    B1 = Matrix{Float64}(undef, nE, nV); B2 = Matrix{Float64}(undef, nF, nE)
    M0 = Matrix{Float64}(undef, nV, nV); M1 = Matrix{Float64}(undef, nE, nE)
    M2 = Matrix{Float64}(undef, nF, nF)
    ccall((:forged_feec_hodge1, LIB), Cint,
          (Ptr{Float64}, Cint, Cint, Ptr{Cint}, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Clonglong},
           Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          coords, nV, dim, sflat, nS, sv, rE, rF, rdd, B1, B2, M0, M1, M2)
    # native writes row-major; Julia read column-major ⇒ each is the transpose of
    # the intended logical matrix, so permutedims recovers it.
    return (B1 = permutedims(B1), B2 = permutedims(B2),
            M0 = permutedims(M0), M1 = permutedims(M1), M2 = permutedims(M2),
            dd2 = Int(rdd[]), nV = nV, nE = nE, nF = nF)
end

# Combinatorial Betti number β_k of the complex over `simplices` (native, unit
# Hodge). The metric-independent ground truth for the FEEC harmonic count.
function betti(simplices, k::Integer)
    nS = length(simplices); sv = length(first(simplices))
    sflat = Vector{Int32}(undef, nS * sv)
    @inbounds for (i, t) in enumerate(simplices)
        s = sort(collect(t)); for j in 1:sv; sflat[(i-1)*sv + j] = Int32(s[j] - 1); end
    end
    Int(ccall((:forged_betti, LIB), Cint, (Ptr{Cint}, Cint, Cint, Cint), sflat, nS, sv, k))
end

# Eigenvalues of the FEEC metric Hodge-1 Laplacian via the EXACT symmetric-
# definite pencil (S + M₁ G M₁) x = λ M₁ x, reduced by Cholesky M₁ = LLᵀ to the
# standard symmetric problem C = L⁻¹(S + M₁GM₁)L⁻ᵀ. G = B₁ᵀ M₀⁻¹ B₁ uses an exact
# Cholesky solve. NO thresholded pseudo-inverse — that was zeroing small-but-
# nonzero modes on obtuse meshes, creating spurious harmonics / corrupting the
# gap. Valid because M₀,M₁ are genuinely SPD on a non-degenerate mesh (Cholesky
# succeeds); a truly isolated DOF would make Cholesky throw, surfacing it.
function hodge1_eigs(points, simplices)
    b = feec_blocks(points, simplices)
    b.dd2 == 0 || error("∂₁∂₂ ≠ 0 (dd_max=$(b.dd2)) — broken complex")
    M0 = Symmetric(Matrix(b.M0)); M1 = Symmetric(Matrix(b.M1)); M2 = Symmetric(Matrix(b.M2))
    G = b.B1' * (cholesky(M0) \ b.B1)            # B₁ᵀ M₀⁻¹ B₁  (nE×nE)
    S = b.B2 * M2 * b.B2'                         # B₂ M₂ B₂ᵀ    (nE×nE)
    A = S + Matrix(M1) * G * Matrix(M1)
    L = cholesky(M1).L
    C = L \ Symmetric((A + A') / 2) / L'          # L⁻¹ A L⁻ᵀ
    return eigvals(Symmetric((C + C') / 2))
end

# Dimension of the FEEC harmonic 1-form space = #{zero eigenvalues of Δ₁} = β₁.
harmonic1_dim(points, simplices; tol = 1e-7) = count(x -> x <= tol, hodge1_eigs(points, simplices))

# SPD-robust matrix square roots (Moore–Penrose on rank-deficient masses, e.g. an
# isolated vertex with zero P1 mass — handled cleanly rather than crashing).
function _sqrt_psd(M)
    E = eigen(Symmetric(Matrix(M)))
    E.vectors * Diagonal(sqrt.(max.(E.values, 0.0))) * E.vectors'
end
function _invsqrt_psd(M; rtol = 1e-9)
    E = eigen(Symmetric(Matrix(M))); λmax = maximum(E.values)
    d = [λ > rtol * λmax ? 1.0 / sqrt(λ) : 0.0 for λ in E.values]
    E.vectors * Diagonal(d) * E.vectors'
end

# Symmetric Hodge–Dirac blocks in M^{1/2} coordinates (all FEEC masses SPD):
#   d̃₀ = M₁^½ d₀ M₀^{-½},  d̃₁ = M₂^½ d₁ M₁^{-½},  d₀=B₁ᵀ, d₁=B₂ᵀ.
# Returns (Δ̃₁, d̃₀, d̃₁) with Δ̃₁ = d̃₀d̃₀ᵀ + d̃₁ᵀd̃₁ symmetric PSD, similar to Δ₁.
# Serves both the spectral gap and the Hodge–Dirac √-relation check.
function hodge1_dirac_blocks(points, simplices)
    b = feec_blocks(points, simplices)
    b.dd2 == 0 || error("∂₁∂₂ ≠ 0 (dd_max=$(b.dd2)) — broken complex")
    M0ih = _invsqrt_psd(b.M0); M1h = _sqrt_psd(b.M1); M1ih = _invsqrt_psd(b.M1); M2h = _sqrt_psd(b.M2)
    d0t = M1h * b.B1' * M0ih                 # B1' = d₀ (nE×nV)
    d1t = M2h * b.B2' * M1ih                 # B2' = d₁ (nF×nE)
    Δ1 = d0t * d0t' + d1t' * d1t
    return Symmetric((Δ1 + Δ1') / 2), d0t, d1t
end

# Spectral gap (smallest nonzero eigenvalue) of the FEEC metric Hodge-1 Laplacian,
# from the exact Cholesky pencil (no pseudo-inverse thresholding artifacts).
function hodge_gap(points, simplices)
    nz = filter(x -> x > 1e-6, hodge1_eigs(points, simplices))
    isempty(nz) ? 0.0 : minimum(nz)
end

# ── FEEC Hodge-0 (P1) gap on the flat torus ℝ²/Λ ────────────────────────────
# N×N periodic grid; all cells congruent ⇒ per-triangle P1 stiffness/mass from the
# reference cell (scaled 1/N), scattered by periodic connectivity. The smallest
# nonzero generalized eigenvalue → the continuum Laplace gap 4π²·λ₁(Λ*)² at O(h²)
# (textbook P1 FEM), tying the discrete FEEC operator to the dual-lattice (=
# smoothing/security) quantity. v1,v2 are the lattice basis vectors.
function torus_hodge0_gap(v1::AbstractVector, v2::AbstractVector, N::Int)
    a = v1 / N; b = v2 / N
    function localKM(E)
        Area = 0.5 * abs(det(E)); G = inv(E)'
        g = [-(G[:,1] + G[:,2]), G[:,1], G[:,2]]
        K = [Area * dot(g[i], g[j]) for i in 1:3, j in 1:3]
        M = [Area * ((i == j) ? 2 : 1) / 12 for i in 1:3, j in 1:3]
        return K, M
    end
    K1, M1 = localKM(hcat(a, b)); K2, M2 = localKM(hcat(b, b - a))
    nv = N * N; idx(i, j) = mod(i, N) + N * mod(j, N) + 1
    K = zeros(nv, nv); M = zeros(nv, nv)
    for i in 0:N-1, j in 0:N-1
        t1 = (idx(i,j), idx(i+1,j), idx(i,j+1)); t2 = (idx(i+1,j), idx(i+1,j+1), idx(i,j+1))
        for p in 1:3, q in 1:3
            K[t1[p],t1[q]] += K1[p,q]; M[t1[p],t1[q]] += M1[p,q]
            K[t2[p],t2[q]] += K2[p,q]; M[t2[p],t2[q]] += M2[p,q]
        end
    end
    L = cholesky(Symmetric(M)).L
    nz = filter(x -> x > 1e-8, eigvals(Symmetric(L \ Symmetric(K) / L')))
    return minimum(nz)
end

# ── Directed-graph magnetic Laplacian L^q (native MagNet) ───────────────────
# Principled operator for an ASYMMETRIC graph: complex Hermitian, direction in
# the U(1) phase Θ^q = 2πq(A_uv − A_vu); real spectrum, no symmetrize hack.
function magnetic_laplacian(A::AbstractMatrix{<:Real}, q::Real)
    n = size(A, 1); @assert size(A, 2) == n
    Af = Matrix{Float32}(undef, n, n)
    @inbounds for i in 1:n, j in 1:n; Af[i, j] = Float32(A[i, j]); end
    Af = collect(Af')                       # row-major for the C ABI
    Lre = Vector{Float32}(undef, n * n); Lim = Vector{Float32}(undef, n * n)
    ccall((:forged_magnetic_laplacian, LIB), Cint,
          (Ptr{Float32}, Cint, Float32, Ptr{Float32}, Ptr{Float32}),
          Af, n, Float32(q), Lre, Lim) == 0 || error("magnetic_laplacian failed")
    L = Matrix{ComplexF64}(undef, n, n)
    @inbounds for i in 1:n, j in 1:n; L[i, j] = ComplexF64(Lre[(i-1)*n + j], Lim[(i-1)*n + j]); end
    return Hermitian(L)
end

function magnetic_gap(A::AbstractMatrix{<:Real}, q::Real)
    nz = filter(x -> x > 1e-6, eigvals(magnetic_laplacian(A, q)))
    isempty(nz) ? 0.0 : minimum(nz)
end

# Directed adjacency from a triangle list, oriented by GEOMETRIC winding (CCW
# positive-signed-area, xy-projection). A−Aᵀ drives the magnetic Laplacian /
# Kitaev i(A−Aᵀ); chirality from the plane, not vertex numbering.
function directed_adjacency(points, tris, nV::Integer)
    A = zeros(nV, nV)
    for t in tris
        a, b, c = t
        s = (points[b][1]-points[a][1])*(points[c][2]-points[a][2]) -
            (points[c][1]-points[a][1])*(points[b][2]-points[a][2])
        p, q, r = s >= 0 ? (a, b, c) : (a, c, b)
        A[p, q] += 1.0; A[q, r] += 1.0; A[r, p] += 1.0
    end
    return A
end

# ── Genuine U(1) magnetic connection (Peierls substitution), any dimension ──
# A real magnetic field is a 2-form B (antisymmetric n×n matrix). In the symmetric
# gauge A(x)=½Bx the Peierls phase along edge p→q is
#     θ_e = ∫_e A·dl = ½ (q−p)ᵀ B p
# (since (q−p)ᵀB(q−p)=0 for antisymmetric B). The flux through a 2-face is the
# oriented sum Σθ_e = ∮A = ∫∫B (discrete Stokes), i.e. the real magnetic flux —
# the same plaquette holonomy native/gauge.h computes. This is the principled
# directed structure in ANY dimension, replacing the ad-hoc xy-projected winding.
# Unique ascending edges of a set of simplices (any dimension).
function edges_of(simplices)
    s = Set{NTuple{2,Int}}()
    for t in simplices, i in 1:length(t), j in (i+1):length(t)
        a, b = minmax(t[i], t[j]); push!(s, (a, b))
    end
    return sort(collect(s))
end

# Antisymmetric field 2-form B = β·(e_1∧e_2) padded to n dimensions (the shear
# plane), the canonical uniform magnetic field for these xy-sheared lattices.
function field_2form(n::Integer, β::Real)
    B = zeros(n, n); B[1, 2] = β; B[2, 1] = -β; return B
end

peierls_phase(p, q, B) = 0.5 * dot(q .- p, B * p)

# Build the adjacency whose magnetic Laplacian carries the Peierls phases: under
# the native convention Θ^1 = 2π(A_uv−A_vu) we set A_uv = w + θ_e/(4π),
# A_vu = w − θ_e/(4π) so that ½(A_uv+A_vu)=w (edge magnitude) and the phase is θ_e.
# Call magnetic_gap(A, 1.0). Geometry-dependent ⇒ varies under shear.
function magnetic_adjacency(points, edges, B; weight = 1.0)
    n = maximum(maximum.(edges)); A = zeros(n, n)
    for (u, v) in edges
        θ = peierls_phase(points[u], points[v], B)
        δ = θ / (4π)
        A[u, v] = weight + δ; A[v, u] = weight - δ
    end
    return A
end

# Net Peierls flux through a 2-simplex (a,b,c): θ_ab + θ_bc + θ_ca = ∮A = ∫B.
face_flux(a, b, c, points, B) =
    peierls_phase(points[a], points[b], B) + peierls_phase(points[b], points[c], B) +
    peierls_phase(points[c], points[a], B)

# ── SU(2) Wilson-loop holonomy via native geom_su2_wilson_f32 ───────────────
# half_angle_vecs: length-3 vectors (θ/2)·n̂ per loop edge. Returns
# (W::NTuple{4}, trace, Θ): Θ the net SU(2) rotation angle, tr=2cos(Θ/2).
function su2_wilson(half_angle_vecs::Vector{<:AbstractVector})
    k = length(half_angle_vecs)
    flat = Vector{Float32}(undef, 3k)
    @inbounds for (i, v) in enumerate(half_angle_vecs), d in 1:3; flat[(i-1)*3 + d] = Float32(v[d]); end
    W = Vector{Float32}(undef, 4); tr = Ref{Float32}(0)
    ccall((:forged_su2_wilson, LIB), Cvoid,
          (Ptr{Float32}, Cint, Ptr{Float32}, Ptr{Float32}), flat, k, W, tr)
    Θ = 2 * atan(norm(@view W[1:3]), W[4])
    return (Tuple(Float64.(W)), Float64(tr[]), Θ)
end

# ── Coxeter–Freudenthal–Kuhn grid triangulation ─────────────────────────────
# The slivers-free triangulation of a (sheared) regular grid: each unit cube →
# D! simplices, a guaranteed-valid simplicial complex in ANY dimension, exact (no
# perturbation). Use this instead of empty-circumsphere "Delaunay", which on a
# cospherical lattice yields overlapping simplices + near-zero-area slivers that
# break the FEEC mass (a skipped degenerate simplex orphans a DOF → singular ★).
all_perms(n::Integer) = n == 1 ? [[1]] :
    reduce(vcat, [[insert!(copy(p), i, n) for i in 1:n] for p in all_perms(n - 1)])

# dims[d] = #grid points along axis d; coord2idx maps a 0-based grid coordinate
# tuple to the 1-based vertex index. Returns sorted (D+1)-simplices.
function kuhn_simplices(dims::Vector{Int}, coord2idx::Dict)
    D = length(dims)
    simplices = NTuple{D + 1, Int}[]
    for corner in Iterators.product((0:dims[d]-2 for d in 1:D)...)
        for perm in all_perms(D)
            c = collect(corner)
            verts = Int[coord2idx[Tuple(c)]]
            for k in 1:D
                c = copy(c); c[perm[k]] += 1
                push!(verts, coord2idx[Tuple(c)])
            end
            push!(simplices, NTuple{D + 1, Int}(sort(verts)))
        end
    end
    return unique(simplices)
end

end # module
