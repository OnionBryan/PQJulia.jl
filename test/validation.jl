#!/usr/bin/env julia
#= test/validation.jl
   The discriminating validation battery for the whole DEC / Hodge / crypto /
   physics stack — every check compares a computed quantity to an INDEPENDENT
   ground truth (topology, an analytic closed form, a published value, or a
   separate numerical oracle), so a real bug fails it. Run: julia test/validation.jl =#
using LinearAlgebra, Random, Printf
include(joinpath(@__DIR__, "forged_dec.jl"));     using .ForgedDEC
include(joinpath(@__DIR__, "lattice_crypto.jl")); using .LatticeCrypto
const FD = ForgedDEC; const LC = LatticeCrypto

pass = 0; fail = 0
chk(c, m) = (global pass, fail; c ? (pass += 1; println("  ok   ", m)) : (fail += 1; println("  FAIL ", m)))

# Test meshes of known topology.
disk(N) = (pts = [[0.0, 0.0]]; for i in 0:N-1; push!(pts, [cos(2π*i/N), sin(2π*i/N)]); end;
           (pts, [NTuple{3,Int}(sort([1, 2+i, 2+(i+1)%N])) for i in 0:N-1]))
function annulus(N; R=2.0, r=1.0)
    p = Vector{Float64}[]
    for i in 0:N-1; θ=2π*i/N; push!(p, [R*cos(θ), R*sin(θ)]); end
    for i in 0:N-1; θ=2π*i/N; push!(p, [r*cos(θ), r*sin(θ)]); end
    t = NTuple{3,Int}[]
    for i in 0:N-1; o0=1+i; o1=1+(i+1)%N; i0=N+1+i; i1=N+1+(i+1)%N
        push!(t, NTuple{3,Int}(sort([o0,o1,i0]))); push!(t, NTuple{3,Int}(sort([o1,i1,i0]))); end
    (p, t)
end
function torus(N, M; R=3.0, r=1.0)
    p = Vector{Float64}[]; id = Dict{NTuple{2,Int},Int}()
    for i in 0:N-1, j in 0:M-1; θ=2π*i/N; φ=2π*j/M
        push!(p, [(R+r*cos(φ))*cos(θ), (R+r*cos(φ))*sin(θ), r*sin(φ)]); id[(i,j)]=length(p); end
    t = NTuple{3,Int}[]
    for i in 0:N-1, j in 0:M-1
        a=id[(i,j)]; b=id[((i+1)%N,j)]; c=id[(i,(j+1)%M)]; d=id[((i+1)%N,(j+1)%M)]
        push!(t, NTuple{3,Int}(sort([a,b,c]))); push!(t, NTuple{3,Int}(sort([b,d,c]))); end
    (p, t)
end

println("── FEEC Hodge theorem  dim ker Δ₁ = β₁  (validates the operator incl. k=2 mass) ──")
for (nm, (p, t), β) in [("disk", disk(12), 0), ("annulus", annulus(16), 1), ("torus", torus(8,6), 2)]
    chk(FD.betti(t,1) == β && FD.harmonic1_dim(p,t) == β, "$nm: native β₁ = FEEC harmonic dim = $β")
end
let P = [[0.0,0,0,0],[1.0,0,0,0],[0.0,1,0,0],[0.0,0,1,0],[0.0,0,0,1]], s = [(1,2,3,4,5)]
    b = FD.feec_blocks(P, s)
    chk(b.dd2 == 0 && isposdef(Symmetric(Matrix(b.M2))), "4D pentatope: ∂∂=0 & Whitney M₂ (k=2,n=4) SPD")
end

println("── Whitney k=2 mass vs independent MC quadrature of ∫⟨W_σ,W_τ⟩ ──")
let P = [[0.0,0,0],[1.0,0.1,0.0],[0.2,1.0,0.1],[0.1,0.2,1.0]]
    M = hcat(P[2]-P[1], P[3]-P[1], P[4]-P[1]); Vol = abs(det(M))/6; G = inv(M)'
    g = [(-(G[:,1]+G[:,2]+G[:,3])), G[:,1], G[:,2], G[:,3]]
    function Wv(σ, λ)
        a, b, c = σ
        2*(λ[a+1]*cross(g[b+1],g[c+1]) - λ[b+1]*cross(g[a+1],g[c+1]) + λ[c+1]*cross(g[a+1],g[b+1]))
    end
    faces = [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]; Mq = zeros(4,4); Random.seed!(7); Ns = 2_000_000
    for _ in 1:Ns
        e = -log.(rand(4)); λ = e./sum(e); Ws = [Wv(f,λ) for f in faces]
        for i in 1:4, j in 1:4; Mq[i,j] += dot(Ws[i],Ws[j]); end
    end
    Mq .*= Vol/Ns
    evn = sort(eigvals(Symmetric(Matrix(FD.feec_blocks(P,[(1,2,3,4)]).M2)))); evq = sort(eigvals(Symmetric(Mq)))
    chk(maximum(abs.(evn.-evq)./abs.(evn)) < 5e-3, "k=2 mass eigvals match quadrature (rel < 5e-3)")
end

println("── Peierls connection: discrete Stokes (face flux = ∮A = ∫∇×A) ──")
let β = 0.3, B = 0.3*[0.0 1; -1 0], a = [0.0,0.0], b = [1.0,0.0], c = [0.3,1.2]
    sa = 0.5*((b[1]-a[1])*(c[2]-a[2]) - (c[1]-a[1])*(b[2]-a[2]))
    chk(isapprox(FD.face_flux(1,2,3,[a,b,c],B), -β*sa, atol=1e-12), "face flux = −β·area exactly")
end

println("── Lattice spectral gap & smoothing parameter vs closed forms ──")
let Z = [1.0 0; 0 1]
    chk(isapprox(LC.torus_spectral_gap(Z), 4π^2, atol=1e-6), "ℤ²: torus Hodge gap = 4π²")
    ε = 2.0^-10; mr = sqrt(log(2*2*(1+1/ε))/π)
    chk(isapprox(LC.smoothing_parameter(Z, ε; R=10), mr, rtol=2e-3),
        "ℤ²: smoothing η_{2⁻¹⁰} ≈ Micciancio–Regev bound $(round(mr,digits=4))")
end

println("── FEEC torus gap → continuum 4π²λ₁(Λ*)² at O(h²) (ties operator to crypto) ──")
let v1 = [1.0,0.0], v2 = [cos(π/3), sin(π/3)]
    cont = LC.torus_spectral_gap(hcat(v1,v2); R=6)
    e16 = abs(FD.torus_hodge0_gap(v1,v2,16) - cont)/cont
    e32 = abs(FD.torus_hodge0_gap(v1,v2,32) - cont)/cont
    chk(e32 < 5e-3 && e16/e32 > 3.0, @sprintf("sheared torus: gap→%.3f, rel err %.1e→%.1e (O(h²))", cont, e16, e32))
end

println("── FEEC gap: Cholesky pencil ≡ Hodge–Dirac similarity (no path discrepancy) ──")
let p = [[0.0,0,0],[1.0,0.1,0],[0.1,1.0,0.05],[0.05,0.1,1.0]], t = [(1,2,3,4)]
    gp = FD.hodge_gap(p, t)
    Δ, _, _ = FD.hodge1_dirac_blocks(p, t)
    gd = minimum(filter(x -> x > 1e-6, eigvals(Δ)))
    chk(isapprox(gp, gd, rtol=1e-8), @sprintf("pencil gap %.6f = Dirac-block gap %.6f", gp, gd))
end

println("── ML-KEM decryption-failure rate vs published δ (noise convolution) ──")
for (lvl, pubn) in ((512,-139), (768,-164), (1024,-174))
    r = LC.kyber_dfr(; LC.KYBER_PARAMS[lvl]...)
    chk(abs(r.log2_dfr - pubn) < 3, @sprintf("Kyber%d: log2(DFR)=%.1f (published %d, within 3 bits)", lvl, r.log2_dfr, pubn))
end

println("── Falcon: σ = σ_min·‖B‖_GS reproduces spec; σ_min = smoothing η'_ε(ℤ) ──")
for lvl in (512, 1024)
    p = LC.FALCON_PARAMS[lvl]
    chk(isapprox(LC.falcon_sigma(p), p.σ, rtol=1e-6) &&
        isapprox(LC.eta_prime(LC.falcon_smoothing_eps(p.σmin)), p.σmin, rtol=1e-9),
        @sprintf("Falcon-%d: σ=%.4f (spec %.4f); σ_min=η'_ε(ℤ)", lvl, LC.falcon_sigma(p), p.σ))
end

println("── QWZ local Chern marker vs known C(m) ──")
include(joinpath(@__DIR__, "qwz_chern.jl"))
for (mass, C) in ((1.0,1), (-1.0,-1), (3.0,0))
    mk = local_chern_markers(16, mass); v = [mk(ix,iy) for ix in 4:11 for iy in 4:11]
    chk(abs(sum(v)/length(v) - C) < 0.15, "QWZ m=$mass: bulk marker ≈ Chern $C")
end

println("\n", fail == 0 ? "ALL $pass VALIDATION CHECKS PASS" : "$fail FAILED / $pass passed")
exit(fail == 0 ? 0 : 1)
