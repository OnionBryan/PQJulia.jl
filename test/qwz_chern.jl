#!/usr/bin/env julia
"""
test/qwz_chern.jl
=================
Validates the Bianco–Resta LOCAL Chern marker against the KNOWN Chern number of
the Qi–Wu–Zhang (QWZ) 2-band Chern insulator — the verification the Thesis-C QWZ
section was missing (it only printed markers, never checked quantization).

QWZ Bloch Hamiltonian (Qi–Wu–Zhang, PRB 74, 085308 / cond-mat/0603414):
    H(k) = sin(k_x) σ_x + sin(k_y) σ_y + (m + cos k_x + cos k_y) σ_z
with the standard Chern phase diagram
    C = +1  for  0 < m < 2
    C = −1  for −2 < m < 0
    C =  0  for |m| > 2.

Real-space (inverse Fourier): onsite m·σ_z; hop r→r+x is T_x = ½(σ_z − i σ_x),
r→r+y is T_y = ½(σ_z − i σ_y) (+ h.c.). Half-filling = occupy the lower band.

Bianco–Resta local marker (PRB 84, 241106 "Mapping topological order in
coordinate space"): 𝔠(r) = −2π i ⟨r|[P X P, P Y P]|r⟩ (trace over the orbitals at
site r; unit cell area = 1). Its BULK average reproduces the global Chern number.
"""

using LinearAlgebra
using Printf

const σx = ComplexF64[0 1; 1 0]
const σy = ComplexF64[0 -im; im 0]
const σz = ComplexF64[1 0; 0 -1]

# Build the L×L open-boundary QWZ Hamiltonian (2 orbitals/site). Returns H and the
# per-orbital (x,y) site coordinates.
function qwz_hamiltonian(L::Int, mass::Real)
    ns = L * L
    site(ix, iy) = iy * L + ix + 1                       # 1-based, ix,iy ∈ 0:L-1
    blk(s) = (2s - 1):(2s)
    H = zeros(ComplexF64, 2ns, 2ns)
    # Fourier convention H[r,r+δ] = coeff of e^{+ik·δ}; the +i sign places C=+1 in
    # 0<m<2 (matching the standard QWZ orientation cond-mat/0603414). The opposite
    # sign just sends k→−k ⇒ C→−C; the quantization/flip/triviality are identical.
    Tx = 0.5 * (σz + im * σx)
    Ty = 0.5 * (σz + im * σy)
    xs = zeros(2ns); ys = zeros(2ns)
    for ix in 0:L-1, iy in 0:L-1
        s = site(ix, iy)
        H[blk(s), blk(s)] .+= mass * σz
        xs[2s-1] = xs[2s] = ix; ys[2s-1] = ys[2s] = iy
        if ix < L - 1
            sx = site(ix + 1, iy)
            H[blk(s), blk(sx)] .+= Tx; H[blk(sx), blk(s)] .+= Tx'
        end
        if iy < L - 1
            sy = site(ix, iy + 1)
            H[blk(s), blk(sy)] .+= Ty; H[blk(sy), blk(s)] .+= Ty'
        end
    end
    return Hermitian(H), xs, ys
end

# Occupied-band projector (lower band, half filling) and the Bianco–Resta marker.
function local_chern_markers(L::Int, mass::Real)
    H, xs, ys = qwz_hamiltonian(L, mass)
    F = eigen(H)
    ns = L * L
    P = zeros(ComplexF64, 2ns, 2ns)
    for m in 1:(ns)                                       # lower half = ns occupied states
        v = @view F.vectors[:, m]; P .+= v * v'
    end
    X = Diagonal(xs); Y = Diagonal(ys)
    PXP = P * X * P; PYP = P * Y * P
    # Bianco–Resta marker 𝔠 = 2π i [P X P, P Y P]; the prefactor SIGN is the
    # orientation convention (x∧y handedness). We fix it to reproduce the standard
    # QWZ C=+1 for 0<m<2 (cond-mat/0603414). Quantization, the sign-flip across
    # m=0, and triviality for |m|>2 are convention-independent.
    M = 2π * im * (PXP * PYP - PYP * PXP)
    site(ix, iy) = iy * L + ix + 1
    marker(ix, iy) = real(M[2*site(ix,iy)-1, 2*site(ix,iy)-1] + M[2*site(ix,iy), 2*site(ix,iy)])
    return marker
end

function run_qwz_validation()
    println("=" ^ 76)
    println("  QWZ Bianco–Resta local Chern marker  vs  known Chern number C(m)")
    println("=" ^ 76)
    L = 16
    @printf("  L=%d open lattice; bulk = central %d×%d block (edges excluded)\n", L, L÷2, L÷2)
    println("-" ^ 76)
    @printf("%-8s | %-22s | %-14s | %-10s\n", "mass m", "bulk marker (mean±sd)", "known C(m)", "verdict")
    println("-" ^ 76)
    allok = true
    for mass in (-3.0, -1.0, 1.0, 3.0, 0.5, -0.5)
        marker = local_chern_markers(L, mass)
        lo, hi = L÷4, 3L÷4 - 1                            # central bulk window
        vals = [marker(ix, iy) for ix in lo:hi for iy in lo:hi]
        μ, sd = sum(vals)/length(vals), std_(vals)
        Cknown = (0 < mass < 2) ? 1 : (-2 < mass < 0) ? -1 : 0
        ok = abs(μ - Cknown) < 0.15                        # quantized to the integer
        allok &= ok
        @printf("%-8.1f | %+8.4f ± %-10.4f | %-14d | %s\n", mass, μ, sd, Cknown, ok ? "✓ quantized" : "✗")
    end
    println("-" ^ 76)
    println(allok ? "  ✓ Bulk marker quantizes to the KNOWN QWZ Chern number across all phases\n    (+1 / −1 topological, 0 trivial) — discriminating, not just 'computed and shown'."
                  : "  ✗ marker does NOT match the known Chern number — investigate.")
    println("=" ^ 76)
    return allok
end

std_(v) = (μ = sum(v)/length(v); sqrt(sum((x-μ)^2 for x in v)/length(v)))

if abspath(PROGRAM_FILE) == @__FILE__
    run_qwz_validation()
end
