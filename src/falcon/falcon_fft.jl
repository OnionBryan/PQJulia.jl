# src/falcon/falcon_fft.jl
# ============================================================================
# Falcon (FN-DSA / FIPS 206) — the FFT and NTT cores over ℤ[x]/(xⁿ+1).
#
# Faithful to the reference (tprest/falcon.py): Falcon uses a complex
# floating-point FFT (NOT the integer NTT used by Kyber/Dilithium) for fast
# Fourier sampling, plus an integer NTT mod q=12289 for the public key h=g/f and
# for verification. This module is the verifiable math core; the sampler,
# ffLDL tree, NTRU keygen, sign/verify build on it.
#
# FFT convention (mathematically fixed, internally consistent): fft(f)[k] =
# f(ζ_k) with ζ_k = exp(iπ(2k+1)/n) the n roots of xⁿ+1, k=0..n-1. Then mul_fft
# is pointwise (evaluation homomorphism) and split/merge_fft recurse to n/2 in
# the SAME convention — validated below against schoolbook negacyclic mult.
# ============================================================================
module FalconFFT

using LinearAlgebra

const q = 12 * 1024 + 1   # 12289

# ── Coefficient-domain helpers (common.py) ──────────────────────────────────
split(f) = (f[1:2:end], f[2:2:end])                       # (even, odd) coeffs
function merge(f0, f1)
    n = 2 * length(f0); f = Vector{eltype(f0)}(undef, n)
    @inbounds for i in 1:length(f0); f[2i-1] = f0[i]; f[2i] = f1[i]; end
    f
end
sqnorm(vs...) = sum(sum(abs2, v) for v in vs)

# Schoolbook negacyclic multiply  a*b mod (xⁿ+1)  (ground truth for validation).
function negacyclic_mul(a::AbstractVector, b::AbstractVector)
    n = length(a); c = zeros(eltype(a), n)
    @inbounds for i in 1:n, j in 1:n
        k = i + j - 2                                       # 0-based degree
        s = (k ÷ n) % 2 == 0 ? 1 : -1                       # x^n = -1 wrap sign
        c[(k % n) + 1] += s * a[i] * b[j]
    end
    c
end

# ── Complex FFT over R[x]/(xⁿ+1) ────────────────────────────────────────────
# ζ_k, the n roots of xⁿ+1; (2k+1)/n is exact for n a power of 2, so cispi rounds each root once.
const ROOTS = Dict{Int,Vector{ComplexF64}}()
roots(n) = get!(() -> ComplexF64[cispi((2k + 1) / n) for k in 0:n-1], ROOTS, n)

# fft(f) = merge_fft(fft(f_even), fft(f_odd)): O(n log n), each root used once per level.
function fft(f::AbstractVector)
    length(f) == 1 && return ComplexF64[f[1]]
    f0, f1 = split(f)
    merge_fft(fft(f0), fft(f1))
end
function ifft(F::AbstractVector)
    length(F) == 1 && return ComplexF64[F[1]]
    F0, F1 = split_fft(F)
    merge(ifft(F0), ifft(F1))
end

# FFT-domain arithmetic (all pointwise).
add_fft(a, b) = a .+ b
sub_fft(a, b) = a .- b
mul_fft(a, b) = a .* b
div_fft(a, b) = a ./ b
adj_fft(a)    = conj.(a)                                     # adjoint: f(x)→f(1/x)=conj on |ζ|=1

# Split/merge in the FFT domain (the recursion ffSampling needs). For F=fft(f)
# of length n, returns (fft(f0), fft(f1)) of length n/2 in the same convention,
# where f0,f1 are the even/odd coefficient halves.  ζ_{k+n/2}² = ζ_k² pairs them.
function split_fft(F::AbstractVector)
    n = length(F); h = n ÷ 2; ζ = roots(n)
    F0 = Vector{ComplexF64}(undef, h); F1 = Vector{ComplexF64}(undef, h)
    @inbounds for k in 1:h
        F0[k] = (F[k] + F[k+h]) / 2
        F1[k] = (F[k] - F[k+h]) * conj(ζ[k]) / 2          # 1/ζ = conj(ζ) on |ζ| = 1
    end
    F0, F1
end
function merge_fft(F0::AbstractVector, F1::AbstractVector)
    h = length(F0); n = 2h; ζ = roots(n)
    F = Vector{ComplexF64}(undef, n)
    @inbounds for k in 1:h
        F[k]   = F0[k] + ζ[k] * F1[k]
        F[k+h] = F0[k] - ζ[k] * F1[k]
    end
    F
end

# ── Integer NTT mod q (negacyclic, for h=g/f and verification) ──────────────
# Primitive 2n-th root ψ of unity mod q (q-1 = 12288 = 2¹²·3 ⇒ 2n | q-1 for
# n ≤ 1024). Found from a generator; validated by the negacyclic-mul check.
function psi_2n(n)
    @assert (q - 1) % (2n) == 0 "2n must divide q-1"
    e = (q - 1) ÷ (2n)
    for g in 2:q-1                                           # smallest g whose g^e has order 2n
        ψ = powermod(g, e, q)
        powermod(ψ, n, q) == q - 1 && return ψ               # ψ^n ≡ -1 ⇒ primitive 2n-th
    end
    error("no primitive 2n-th root")
end

# Negacyclic NTT: weight by ψ^i then length-n NTT with ω=ψ². Returns the
# point-values; intt inverts. mul mod (xⁿ+1) mod q = intt(ntt(a).*ntt(b)).
function ntt_mul(a::Vector{<:Integer}, b::Vector{<:Integer})
    n = length(a); ψ = psi_2n(n); ψi = invmod(ψ, q); ni = invmod(n, q)
    â = [mod(a[i+1] * powermod(ψ, i, q), q) for i in 0:n-1]
    b̂ = [mod(b[i+1] * powermod(ψ, i, q), q) for i in 0:n-1]
    A = _ntt(â, powermod(ψ, 2, q)); B = _ntt(b̂, powermod(ψ, 2, q))
    C = [mod(A[i] * B[i], q) for i in 1:n]
    c = _ntt(C, powermod(ψi, 2, q))
    return [mod(c[i+1] * ni % q * powermod(ψi, i, q), q) for i in 0:n-1]
end
function _ntt(a, ω)                                          # O(n²) DFT mod q (reference)
    n = length(a)
    [mod(sum(a[j+1] * powermod(ω, i*j, q) for j in 0:n-1), q) for i in 0:n-1]
end

# ── O(n log n) iterative Cooley–Tukey NTT mod q (exact; fills the perf gap) ──
function _bitrev!(a)
    n = length(a); j = 0
    for i in 1:n-1
        bit = n >> 1
        while j & bit != 0; j ⊻= bit; bit >>= 1; end
        j ⊻= bit
        if i < j; a[i+1], a[j+1] = a[j+1], a[i+1]; end
    end
    a
end
function ntt_ct!(a::Vector{Int}, ω::Int)                    # in-place forward NTT, root ω (order n)
    n = length(a); _bitrev!(a); len = 2
    while len <= n
        wlen = powermod(ω, n ÷ len, q)
        for i in 0:len:n-1
            w = 1
            for k in 0:(len ÷ 2 - 1)
                u = a[i+k+1]; v = Int(mod(a[i+k+len÷2+1] * w, q))
                a[i+k+1] = mod(u + v, q); a[i+k+len÷2+1] = mod(u - v + q, q)
                w = Int(mod(w * wlen, q))
            end
        end
        len <<= 1
    end
    a
end
# Negacyclic multiply a·b mod (xⁿ+1) mod q via the fast NTT (weight by ψ^i).
function ntt_mul_fast(a::Vector{<:Integer}, b::Vector{<:Integer})
    n = length(a); ψ = psi_2n(n); ω = powermod(ψ, 2, q)
    ψi = invmod(ψ, q); ωi = invmod(ω, q); ni = invmod(n, q)
    â = [Int(mod(a[i+1] * powermod(ψ, i, q), q)) for i in 0:n-1]
    b̂ = [Int(mod(b[i+1] * powermod(ψ, i, q), q)) for i in 0:n-1]
    ntt_ct!(â, ω); ntt_ct!(b̂, ω)
    C = [Int(mod(â[i] * b̂[i], q)) for i in 1:n]
    ntt_ct!(C, ωi)
    return [Int(mod(C[i+1] * ni % q * powermod(ψi, i, q), q)) for i in 0:n-1]
end

# ── Self-validation (run on include) ────────────────────────────────────────
function validate(; verbose=true)
    ok = true
    say(c, m) = (verbose && println("  ", c ? "ok  " : "FAIL", "  ", m); ok &= c)
    for n in (2, 4, 8, 16)
        a = randn(n); b = randn(n)
        say(maximum(abs.(ifft(fft(a)) .- a)) < 1e-9, "n=$n: ifft∘fft = id")
        # evaluation homomorphism: mul_fft = negacyclic mult
        c_fft = real.(ifft(mul_fft(fft(a), fft(b))))
        say(maximum(abs.(c_fft .- negacyclic_mul(a, b))) < 1e-8, "n=$n: mul_fft = neg-cyclic mult")
        # split/merge roundtrip and consistency with coefficient split
        F = fft(a); F0, F1 = split_fft(F)
        say(maximum(abs.(merge_fft(F0, F1) .- F)) < 1e-9, "n=$n: merge∘split_fft = id")
        a0, a1 = split(a)
        say(maximum(abs.(F0 .- fft(a0))) < 1e-8 && maximum(abs.(F1 .- fft(a1))) < 1e-8,
            "n=$n: split_fft gives (fft(even), fft(odd))")
        # NTT negacyclic mult mod q
        ia = rand(0:q-1, n); ib = rand(0:q-1, n)
        say(ntt_mul(ia, ib) == mod.(negacyclic_mul(ia, ib), q), "n=$n: NTT mult = neg-cyclic mod q")
    end
    return ok
end

end # module
