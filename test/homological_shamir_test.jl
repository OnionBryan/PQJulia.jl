#!/usr/bin/env julia
"""
PQJulia.jl Homological Secret Sharing Tests
==========================================
A rigorous homological verification suite for Shamir (k,n)-threshold sharing.
This suite verifies that the threshold properties are isomorphic to:
1. The dimension of the homology group H_0 (kernel) of the evaluation map.
2. The surjectivity of the secret projection operator (Shannon perfect secrecy).
3. The exactness of the full chain complex H_1 = ker H / im G = 0.
"""

using PQJulia
using Test
using Random

const PRIME_P = big(2)^521 - 1

# ── Finite Field Arithmetic Helpers ──────────────────────────────────────────

# Gaussian elimination over GF(p) to compute RREF, rank, and pivots
function gf_rref(M::Matrix{BigInt}, p::BigInt)
    rows, cols = size(M)
    A = copy(M)
    r = 0
    pivots = Int[]
    for c in 1:cols
        pivot_row = 0
        for r_idx in (r+1):rows
            if A[r_idx, c] != 0
                pivot_row = r_idx
                break
            end
        end
        if pivot_row == 0
            continue
        end
        r += 1
        push!(pivots, c)
        if r != pivot_row
            # Swap rows
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
        for r_idx in 1:rows
            r_idx == r && continue
            factor = A[r_idx, c]
            factor == 0 && continue
            for j in c:cols
                A[r_idx, j] = mod(A[r_idx, j] - factor * A[r, j], p)
            end
        end
    end
    return r, pivots, A
end

# Compute basis of the kernel (homology H_0) from RREF
function gf_kernel(M::Matrix{BigInt}, p::BigInt)
    rows, cols = size(M)
    rank, pivots, rref_M = gf_rref(M, p)
    free_vars = setdiff(1:cols, pivots)
    basis = Vector{BigInt}[]
    for f in free_vars
        v = zeros(BigInt, cols)
        v[f] = 1
        for (r_idx, p_col) in enumerate(pivots)
            v[p_col] = mod(-rref_M[r_idx, f], p)
        end
        push!(basis, v)
    end
    return basis
end

# Compute the rank of a matrix over GF(p)
function gf_rank(M::Matrix{BigInt}, p::BigInt)
    rank, _, _ = gf_rref(M, p)
    return rank
end

# Helper to construct subset matrices
function get_subset_matrix(G::Matrix{BigInt}, subset::AbstractVector{Int})
    return G[subset, :]
end

# ── Lagrange Parity Check Construction ────────────────────────────────────────
# Constructs the dual parity-check matrix H ((n-k) x n) over GF(p)
# using Lagrange interpolation coefficients.
function construct_parity_check(k::Int, n::Int, p::BigInt)
    H = zeros(BigInt, n - k, n)
    for row in 1:(n - k)
        # Select indices {1..k} and the extra index k+row
        indices = vcat(1:k, k + row)
        for i in indices
            den = big(1)
            for j in indices
                i == j && continue
                # Note: x_i = i, x_j = j
                den = mod(den * (i - j), p)
            end
            c_i = PQJulia.mod_inverse(den, p)
            H[row, i] = c_i
        end
    end
    return H
end

# Simple helper for combinations since we don't import Combinatorics.jl
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

# ── Main Test Suite ──────────────────────────────────────────────────────────

@testset "Rigorous Homological Threshold Analysis" begin
    p = PRIME_P

    # Sweep parameters k and n to ensure global algebraic invariants hold
    # without relying on specific parameter configurations.
    for n in 2:6
        for k in 1:n
            @testset "Threshold (k=$k, n=$n) Homological Profiling" begin
                
                # 1. Construct the evaluation matrix G (n x k)
                # G[i, j] = i^(j-1) mod p
                G = Matrix{BigInt}(undef, n, k)
                for i in 1:n
                    for j in 1:k
                        G[i, j] = mod(big(i)^(j-1), p)
                    end
                end

                # 2. Sweep all possible subsets of shares I ⊆ {1..n}
                # Verifying that dim H_0 = max(0, k - |I|) globally.
                for m in 0:n
                    # Generate all combinations of size m
                    for subset in combinations(1:n, m)
                        G_sub = get_subset_matrix(G, subset)
                        basis = gf_kernel(G_sub, p)
                        expected_dim = max(0, k - m)
                        
                        # Verify exact kernel dimension (H_0 rank identity)
                        @test length(basis) == expected_dim

                        # 3. Secrecy Projection Verification (Shannon Perfect Secrecy)
                        # If m < k, the projection map π_1 : ker G_sub -> GF(p)
                        # onto the first coordinate (the secret) must be surjective (rank 1).
                        if m < k
                            # To prove surjectivity, we must show there exists a vector v
                            # in the kernel where the first coordinate is non-zero.
                            # Since the field is large, we check if the first coordinate
                            # is not locked to zero across the kernel space.
                            has_secret_freedom = false
                            for v in basis
                                if v[1] != 0
                                    has_secret_freedom = true
                                    break
                                end
                            end
                            if !has_secret_freedom && length(basis) > 0
                                # Fallback: check if any linear combination yields non-zero
                                proj_sum = sum(v[1] for v in basis)
                                @test proj_sum != 0
                            else
                                @test has_secret_freedom
                            end
                        end
                    end
                end

                # 4. Chain Complex Exactness (H_1 = ker H / im G = 0)
                # If n > k, construct the parity-check matrix H and check exactness.
                if n > k
                    H = construct_parity_check(k, n, p)
                    
                    # Verify that H * G = 0 (image of G is inside kernel of H)
                    HG = mod.(H * G, p)
                    @test all(HG .== 0)

                    # Verify exactness at F_p^n:
                    # dim(im G) + dim(ker H) = n (by rank-nullity)
                    # We check that rank(G) == dim(ker H)
                    rank_G = gf_rank(G, p)
                    basis_H_ker = gf_kernel(H, p)
                    
                    @test rank_G == k
                    @test length(basis_H_ker) == k
                    
                    # This implies H_1 = ker H / im G has dimension 0 (exactness)
                    # proving the cochain complex is acyclic.
                end
            end
        end
    end
end


println("=" ^ 70)
println("ALL RIGOROUS HOMOLOGICAL INVARIANT TESTS PASSED")
println("=" ^ 70)
