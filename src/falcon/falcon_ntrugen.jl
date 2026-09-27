# src/falcon/falcon_ntrugen.jl
# ============================================================================
# Falcon — NTRU key generation: solve  f·G − g·F = q  in ℤ[x]/(xⁿ+1).
# Faithful to tprest/falcon.py ntrugen.py. The tower recursion (field_norm /
# lift / galois_conjugate + xgcd base) produces (F,G) satisfying the NTRU
# equation EXACTLY; Babai `reduce` then shrinks (F,G) for short keys (it
# subtracts k·f from F and k·g from G, preserving f·G−g·F = q). BigInt
# throughout — F,G grow large before reduction.
#
# Lift identity (why the recursion is correct): lift is a ring homomorphism and
# f·f̄ = lift(field_norm(f)); with F=lift(Fp)·ḡ, G=lift(Gp)·f̄,
#   fG − gF = (f f̄)·lift(Gp) − (g ḡ)·lift(Fp) = lift(fp·Gp − gp·Fp) = lift(q) = q.
# ============================================================================
module FalconNTRUGen

import ..FalconFFT as FF

const q = 12289

# Exact negacyclic multiply a·b mod (xⁿ+1) over ℤ: Karatsuba full product, then fold x^n = −1.
function negamul(a::AbstractVector, b::AbstractVector)
    n = length(a)
    n < 32 && return negamul_school(a, b)
    ab = karamul(BigInt.(a), BigInt.(b))
    [ab[i] - ab[i+n] for i in 1:n]
end

# Full (2n-length) product of two length-n (n a power of 2) coefficient vectors, as ntrugen.py karamul.
function karamul(a::Vector{BigInt}, b::Vector{BigInt})
    n = length(a)
    if n <= 16
        c = zeros(BigInt, 2n)
        @inbounds for i in 1:n, j in 1:n
            c[i+j-1] += a[i] * b[j]
        end
        return c
    end
    h = n ÷ 2
    a0, a1, b0, b1 = a[1:h], a[h+1:n], b[1:h], b[h+1:n]
    a0b0 = karamul(a0, b0); a1b1 = karamul(a1, b1)
    mid = karamul(a0 .+ a1, b0 .+ b1) .- a0b0 .- a1b1
    c = zeros(BigInt, 2n)
    c[1:n] .+= a0b0; c[n+1:2n] .+= a1b1; c[h+1:h+n] .+= mid
    c
end

# Schoolbook negacyclic multiply (small n, and the reference for negamul).
function negamul_school(a::AbstractVector, b::AbstractVector)
    n = length(a); c = zeros(BigInt, n)
    @inbounds for i in 1:n, j in 1:n
        k = i + j - 2
        s = (k ÷ n) % 2 == 0 ? 1 : -1
        c[(k % n) + 1] += s * BigInt(a[i]) * BigInt(b[j])
    end
    c
end

galois_conjugate(a) = [iseven(i-1) ? a[i] : -a[i] for i in 1:length(a)]   # a(x)→a(−x)

function lift(a)                              # a(x) → a(x²): coeffs at even positions
    n = length(a); r = zeros(BigInt, 2n)
    @inbounds for i in 1:n; r[2i-1] = a[i]; end
    r
end

# field_norm: ℤ[x]/(xⁿ+1) → ℤ[y]/(y^{n/2}+1),  N(f) = ae(y)² − y·ao(y)².
function field_norm(f)
    n = length(f); h = n ÷ 2
    ae = [f[2i-1] for i in 1:h]; ao = [f[2i] for i in 1:h]
    ae2 = negamul(ae, ae); ao2 = negamul(ao, ao)
    yao2 = zeros(BigInt, h)                    # y·ao²  (shift up by 1, negacyclic)
    @inbounds for i in 1:h
        if i == 1; yao2[1] = -ao2[h]           # x^{h} = −1 wrap
        else; yao2[i] = ao2[i-1]; end
    end
    return ae2 .- yao2
end

# Extended gcd over ℤ:  u·a + v·b = g.
function xgcd(a::Integer, b::Integer)
    old_r, r = BigInt(a), BigInt(b); old_s, s = BigInt(1), BigInt(0); old_t, t = BigInt(0), BigInt(1)
    while r != 0
        Q = old_r ÷ r
        old_r, r = r, old_r - Q*r
        old_s, s = s, old_s - Q*s
        old_t, t = t, old_t - Q*t
    end
    return old_r, old_s, old_t      # g, u, v
end

# Solve f·G − g·F = q recursively. Returns (F,G) or throws if the base gcd ≠ 1
# (resultants not coprime ⇒ caller resamples f,g).
function ntru_solve(f, g)
    n = length(f)
    if n == 1
        d, u, v = xgcd(f[1], g[1])
        d == 1 || throw(ErrorException("non-coprime resultants"))
        return BigInt[-q*v], BigInt[q*u]
    end
    fp = field_norm(f); gp = field_norm(g)
    Fp, Gp = ntru_solve(fp, gp)
    F = negamul(lift(Fp), galois_conjugate(g))
    G = negamul(lift(Gp), galois_conjugate(f))
    return reduce!(f, g, F, G)             # per level, as ntrugen.py: keeps (F,G) near (f,g) size
end

# Babai reduce: (F,G) −= k·(f,g), k = round((F f̄ + G ḡ)/(f f̄ + g ḡ)) in FFT on the top
# 53 bits, repeated until (F,G) is within 8 bits of (f,g). Preserves f·G − g·F = q.
maxbits(v) = maximum(x -> (x == 0 ? 0 : ndigits(abs(x), base=2)), v)
function reduce!(f, g, F, G)
    fb = max(maxbits(f), maxbits(g))
    shf = max(fb - 53, 0)
    ff = FF.fft(Float64.(f .>> shf)); gf = FF.fft(Float64.(g .>> shf))
    den = ff .* conj.(ff) .+ gf .* conj.(gf)
    while true
        sz = max(maxbits(F), maxbits(G))
        sz <= fb + 8 && break
        sh = max(sz - 53, 0)
        Ff = FF.fft(Float64.(F .>> sh)); Gf = FF.fft(Float64.(G .>> sh))
        k = round.(BigInt, real.(FF.ifft((Ff .* conj.(ff) .+ Gf .* conj.(gf)) ./ den)))
        all(iszero, k) && break
        kb = k .<< (sh - shf)
        F .-= negamul(kb, f); G .-= negamul(kb, g)
    end
    F, G
end

# Verify the NTRU equation exactly.
ntru_check(f, g, F, G) = negamul(f, G) .- negamul(g, F) == [BigInt(q); zeros(BigInt, length(f)-1)]

end # module
