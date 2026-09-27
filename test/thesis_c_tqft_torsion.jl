#!/usr/bin/env julia
"""
test/thesis_c_tqft_torsion.jl
=============================
Thesis C Experiment: Digraph Asymmetry, Reidemeister Torsion, and Chirality.
This script:
1. Constructs the bipartite evaluation digraph of Shamir's secret sharing scheme.
2. Computes the symmetric and asymmetric parts of the digraph adjacency matrix.
3. Computes the exact Reidemeister torsion of the cochain complex over GF(p).
4. Verifies that permuting shares (chirality swap) flips the sign of the torsion.
"""

using PQJulia
using Test
using LinearAlgebra
using Printf

const PRIME_P = big(2)^521 - 1

# ── 1. Finite Field Linear Algebra Helpers ────────────────────────────────────

# Computes the determinant of a square matrix M over GF(p) using Gaussian elimination
function gf_det(M::Matrix{BigInt}, p::BigInt)
    n, m = size(M)
    @assert n == m "Matrix must be square"
    A = copy(M)
    det_val = big(1)
    for i in 1:n
        pivot_row = 0
        for r_idx in i:n
            if A[r_idx, i] != 0
                pivot_row = r_idx
                break
            end
        end
        if pivot_row == 0
            return big(0)
        end
        if i != pivot_row
            # Swap rows changes sign of determinant
            tmp = A[i, :]
            A[i, :] = A[pivot_row, :]
            A[pivot_row, :] = tmp
            det_val = mod(-det_val, p)
        end
        # Multiply by pivot
        det_val = mod(det_val * A[i, i], p)
        inv_pivot = PQJulia.mod_inverse(A[i, i], p)
        # Eliminate other rows below pivot
        for r_idx in (i+1):n
            factor = A[r_idx, i]
            factor == 0 && continue
            for j in i:n
                A[r_idx, j] = mod(A[r_idx, j] - factor * inv_pivot * A[i, j], p)
            end
        end
    end
    return det_val
end

# Computes the matrix inverse over GF(p) using Gaussian elimination
function gf_inv(M::Matrix{BigInt}, p::BigInt)
    n, m = size(M)
    @assert n == m "Matrix must be square"
    A = [M Matrix{BigInt}(I, n, n)]
    rows, cols = size(A)
    r = 0
    for c in 1:n
        pivot_row = 0
        for r_idx in (r+1):n
            if A[r_idx, c] != 0
                pivot_row = r_idx
                break
            end
        end
        if pivot_row == 0
            error("Matrix is singular over GF(p)")
        end
        r += 1
        if r != pivot_row
            tmp = A[r, :]
            A[r, :] = A[pivot_row, :]
            A[pivot_row, :] = tmp
        end
        # Scale pivot to 1
        inv_pivot = PQJulia.mod_inverse(A[r, c], p)
        for j in c:cols
            A[r, j] = mod(A[r, j] * inv_pivot, p)
        end
        # Eliminate other rows
        for r_idx in 1:n
            r_idx == r && continue
            factor = A[r_idx, c]
            factor == 0 && continue
            for j in c:cols
                A[r_idx, j] = mod(A[r_idx, j] - factor * A[r, j], p)
            end
        end
    end
    return A[:, (n+1):end]
end

# Lagrange parity-check matrix H ((n-k) x n) over GF(p)
function construct_parity_check(k::Int, n::Int, p::BigInt)
    H = zeros(BigInt, n - k, n)
    for row in 1:(n - k)
        indices = vcat(1:k, k + row)
        for i in indices
            den = big(1)
            for j in indices
                i == j && continue
                den = mod(den * (i - j), p)
            end
            c_i = PQJulia.mod_inverse(den, p)
            H[row, i] = c_i
        end
    end
    return H
end

# ── 2. Digraph & Torsion Experiment Harness ───────────────────────────────────

function run_thesis_c_experiment()
    println("=" ^ 80)
    println("  Thesis C: Digraph Asymmetry & Reidemeister Torsion Chirality")
    println("=" ^ 80)
    
    p = PRIME_P
    k = 3
    n = 5
    
    # 1. Construct the evaluation matrix G (n x k)
    G = Matrix{BigInt}(undef, n, k)
    for i in 1:n
        for j in 1:k
            G[i, j] = mod(big(i)^(j-1), p)
        end
    end
    
    # 2. Build adjacency matrix A of the bipartite digraph
    # Left nodes (1..k): coefficients of f(x)
    # Right nodes (k+1..k+n): shares s_i
    dim_A = k + n
    A = zeros(BigInt, dim_A, dim_A)
    # Directed edges from coefficients to shares
    for i in 1:n
        for j in 1:k
            A[k + i, j] = G[i, j]
        end
    end
    
    # 3. Decompose A into symmetric and asymmetric parts
    # A_sym = 1/2 * (A + A^T), A_asym = 1/2 * (A - A^T)
    inv_2 = PQJulia.mod_inverse(big(2), p)
    A_sym = mod.(inv_2 * (A + A'), p)
    A_asym = mod.(inv_2 * (A - A'), p)
    
    # Verify decomposition
    @assert mod.(A_sym + A_asym, p) == A "Digraph decomposition check failed!"
    # Verify symmetry & skew-symmetry
    @assert A_sym == mod.(A_sym', p) "Symmetric part is not symmetric!"
    @assert A_asym == mod.(-A_asym', p) "Asymmetric part is not skew-symmetric!"
    println("✓ Bipartite digraph adjacency matrix successfully decomposed.")
    println("  - A_sym represents undirected coefficient-share correlations.")
    println("  - A_asym represents the oriented information flow (causality).")
    
    # 4. Compute Reidemeister Torsion of the Exact Chain Complex
    H = construct_parity_check(k, n, p)
    @assert all(mod.(H * G, p) .== 0) "Chain complex exactness check failed (H * G != 0)!"
    
    # Right inverse s of H: s = H^T * (H * H^T)^-1
    HH_T = mod.(H * H', p)
    inv_HH_T = gf_inv(HH_T, p)
    s = mod.(H' * inv_HH_T, p)
    @assert mod.(H * s, p) == Matrix{BigInt}(I, n-k, n-k) "Right inverse check failed!"
    
    # Assemble transition matrix M = [G | s]
    M = [G s]
    tau = gf_det(M, p)
    println("\n✓ Reidemeister Torsion (τ) computed exactly: $tau")
    
    # 5. Permute shares to test Chirality (Orientation swap)
    # Swap share 1 and share 2 (rows 1 and 2 in G, columns 1 and 2 in H)
    P = Matrix{BigInt}(I, n, n)
    P[1, 1] = 0; P[2, 2] = 0
    P[1, 2] = 1; P[2, 1] = 1 # Transposition (odd permutation, sign = -1)
    
    G_perm = mod.(P * G, p)
    H_perm = mod.(H * P', p) # permute columns of H
    
    HH_T_perm = mod.(H_perm * H_perm', p)
    inv_HH_T_perm = gf_inv(HH_T_perm, p)
    s_perm = mod.(H_perm' * inv_HH_T_perm, p)
    
    M_perm = [G_perm s_perm]
    tau_perm = gf_det(M_perm, p)
    
    println("✓ Torsion after share transposition (τ_perm): $tau_perm")
    
    # Check that τ_perm == -τ mod p
    @assert tau_perm == mod(-tau, p) "Chirality sign swap check failed!"
    println("✓ Chirality conservation verified (τ_perm == -τ mod p holds exactly).")
    println("================================================================================")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_thesis_c_experiment()
end
