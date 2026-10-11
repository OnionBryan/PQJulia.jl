# FALCON-MRM, ePrint 2026/420. The spec publishes no test vectors and leaves the points below
# open; the choices are PQJulia's own, so signatures will not interoperate with other
# implementations until a revision of the spec fixes them.
#   domain separation  TAG_H1/TAG_H2 prefixed to the HashToPoint input in place of the salt
#   H1 input framing   pack14(ρ ‖ M1) as one 14-bit stream, then M2
#   M2                 bytes (the spec allows any bit string)
#   header byte        0x70 + log₂n (cc = 11, unused in round 3); sig = header ‖ e1 ‖ e2
#   ρ sampler          16-bit big-endian draws below 5q, reduced mod q (uniform_zq)
#   compression        either half longer than slen bytes: new attempt, fresh ρ and sampler seed
#   M1 range           sign rejects M1ᵢ ∉ [0, q)
#   odd k (§1.4)       encode sets the one element of ℤ_q^k that carries no bits to 0; decode requires 0
module FalconMRM

import ..Falcon
import ..FalconChaCha as CC
import ..FalconEncoding as FE
import ...Wipe: wipe!
import ...Bytes: catbytes

const q = 12289

# §1.1 table: λ, γ smallest allowed by Lemma 1, Eq. (1).
const PARAMS = Dict(512 => (λ = 19, γ = 15), 1024 => (λ = 38, γ = 24))
params(n) = haskey(PARAMS, n) ? PARAMS[n] : throw(ArgumentError("FALCON-MRM is defined for n = 512, 1024"))

const TAG_H1 = Vector{UInt8}("PQJulia FALCON-MRM H1")
const TAG_H2 = Vector{UInt8}("PQJulia FALCON-MRM H2")

m1_len(n) = (p = params(n); n - p.λ - p.γ)
slen(n) = Falcon.sig_bytes(n) - 1 - Falcon.SALT_LEN            # 625, 1239 (§1.3)
sig_bytes(n) = 1 + 2 * slen(n)                                 # 1251, 2479 (Table 1)
header(n) = UInt8(0x70 + Falcon.logn(n))

# ℤ_q vectors as 14-bit big-endian fields, zero-padded in the last byte's low bits (§1.1).
pack14(v) = Falcon.pack_bits(v, 14)
# HashToPoint with the input prefixed by a tag (§1.1); ρM1 is (ρ, M1) as one vector.
H1(ρM1, M2, λ) = Falcon.hash_to_point(catbytes(pack14(ρM1), M2), TAG_H1, λ)
H2(c1, k) = Falcon.hash_to_point(pack14(c1), TAG_H2, k)

# ρ ← ℤ_q^γ uniformly: 16-bit draws below 5q, reduced mod q.
function uniform_zq(k, randombytes)
    out = Int[]
    while length(out) < k
        b = randombytes(2); x = (Int(b[1]) << 8) | b[2]
        x < 5q && push!(out, x % q)
    end
    out
end

"Throw `ArgumentError` unless `M1` ∈ ℤ_q^m1_len(n): `m1_len(n)` elements, each in [0, q)."
function check_m1(M1::AbstractVector{<:Integer}, n::Int)
    k = m1_len(n)
    length(M1) == k || throw(ArgumentError("M1 must have $k elements of ℤ_q"))
    all(x -> 0 <= x < q, M1) || throw(ArgumentError("M1 elements must lie in [0, q)"))
    M1
end

"""
    sign(M1, M2, ek::Falcon.ExpandedKey; randombytes) -> Vector{UInt8}

Algorithm 3. A fresh ρ per attempt; returns header ‖ Compress(s1) ‖ Compress(s2), where
s1 + s2·h = c = (H1(ρ, M), H2(c1) + (ρ, M1)). Each attempt's sampler seed, ChaCha20 state and
keystream, targets t0, t1 (which give f and F from the public c) and candidate (s1, s2) are wiped
when the attempt ends: rejected, accepted or by an exception.
"""
function sign(M1::AbstractVector{<:Integer}, M2::AbstractVector{UInt8}, ek::Falcon.ExpandedKey;
              randombytes=Falcon.sysrandom)
    n = ek.sk.n; (; λ, γ) = params(n); fp = Falcon.params(n); L = slen(n) * 8
    check_m1(M1, n)
    Falcon.check_leaves(ek.gs, fp.σmin)
    ρM1 = Vector{Int}(undef, n - λ)                    # (ρ, M1), indexed from 1 whatever M1's axes
    for (i, x) in enumerate(M1)
        ρM1[γ + i] = x
    end
    while true
        copyto!(ρM1, 1, uniform_zq(γ, randombytes), 1, γ)
        c1 = H1(ρM1, M2, λ)
        h2 = H2(c1, n - λ)
        c = Vector{Int}(undef, n)                       # (c1, H2(c1) + (ρ, M1) mod q)
        copyto!(c, 1, c1, 1, λ)
        for i in 1:n-λ
            c[λ + i] = mod(h2[i] + ρM1[i], q)
        end
        sig = attempt(ek.gs, c, n, fp, L, randombytes)
        sig === nothing || return sig
    end
end

# One ffSampling attempt for the target c: header ‖ e1 ‖ e2, or `nothing` when the candidate is too
# long or does not compress (sign then draws a new ρ). t0, t1, the candidate and its encodings are
# wiped on every exit, including the rejections and an exception.
function attempt(gs, c::Vector{Int}, n::Int, fp, L::Int, randombytes)
    t0, t1 = Falcon.target(gs, c)
    try
        s1, s2 = sample_preimage(gs, c, t0, t1, fp.σmin, randombytes)
        e1 = e2 = nothing
        try
            sum(abs2, s1) + sum(abs2, s2) <= fp.β2 || return nothing
            e1 = FE.compress(s1, L); e1 === nothing && return nothing
            e2 = FE.compress(s2, L); e2 === nothing && return nothing
            return catbytes(header(n), e1, e2)
        finally
            wipe!(s1, s2)
            e1 === nothing || wipe!(e1)                 # a rejected attempt's e1 encodes its s1
            e2 === nothing || wipe!(e2)
        end
    finally
        wipe!(t0, t1)
    end
end

# ffSampling driven by ChaCha20 on a fresh 56-byte seed, as in Falcon.sign_poly; the seed is wiped
# once the generator holds it, and the generator (state and keystream) once sampling ends or throws.
function sample_preimage(gs, c, t0, t1, σmin, randombytes)
    seed = randombytes(Falcon.SEED_LEN)
    r = try CC.ChaCha20(seed) finally wipe!(seed) end
    try
        Falcon.preimage(gs, c, t0, t1, σmin, r)
    finally
        wipe!(r)
    end
end

"""
    verify(M2, sig, h) -> Union{Vector{Int}, Nothing}

Algorithm 4. Returns the recovered M1, or `nothing` (⊥).
"""
function verify(M2::AbstractVector{UInt8}, sig::AbstractVector{UInt8}, h::AbstractVector{<:Integer})
    n = length(h); haskey(PARAMS, n) || return nothing
    (; λ, γ) = params(n); fp = Falcon.params(n); ℓ = slen(n); o = firstindex(sig)
    (length(sig) == sig_bytes(n) && sig[o] == header(n)) || return nothing
    s1 = FE.decompress(@view(sig[o+1:o+ℓ]), 8ℓ, n); s1 === nothing && return nothing
    s2 = FE.decompress(@view(sig[o+1+ℓ:o+2ℓ]), 8ℓ, n); s2 === nothing && return nothing
    sum(abs2, s1) + sum(abs2, s2) <= fp.β2 || return nothing
    hq = Vector{Int}(undef, n)                          # h mod q, indexed from 1 whatever h's axes
    for (i, x) in enumerate(h)
        hq[i] = mod(x, q)
    end
    c = mod.(s1 .+ Falcon.poly_mul_modq(s2, hq), q)
    c1 = c[1:λ]
    ρM1 = mod.(c[λ+1:end] .- H2(c1, n - λ), q)
    H1(ρM1, M2, λ) == c1 ? ρM1[γ+1:end] : nothing
end

# ── §1.4 bit strings ↔ ℤ_q^k: 163 bits per 12 elements, then 27 bits per 2 ──
# For odd k the last element carries no bits: encode sets it to 0 and decode requires 0.
"B(k): bits encodable into ℤ_q^k by `encode` (6492 for k = 478, 13067 for k = 962)."
max_bits(k::Integer) = k >= 0 ? 163 * (k ÷ 12) + 27 * ((k % 12) ÷ 2) :
                                throw(ArgumentError("ℤ_q^k needs k ≥ 0, got k = $k"))

bits_to_int(b) = foldl((x, bit) -> (x << 1) | big(bit), b; init = big(0))
int_to_bits(x, ℓ) = [isodd(x >> (ℓ - 1 - i)) for i in 0:ℓ-1]
int_to_base(x, k) = (z = Int[]; for _ in 1:k; push!(z, Int(mod(x, q))); x = div(x, q); end; z)
base_to_int(z) = foldr((zi, x) -> x * q + zi, z; init = big(0))

"""Algorithm 7: encode exactly `max_bits(k)` bits into ℤ_q^k, a vector of k elements (for odd k the
last is 0); `nothing` on a length mismatch."""
function encode(μ::AbstractVector{Bool}, k::Integer)
    B = length(μ); B == max_bits(k) || return nothing
    z = Int[]; i = 0; o = firstindex(μ) - 1
    while B - i >= 163
        append!(z, int_to_base(bits_to_int(view(μ, o+i+1:o+i+163)), 12)); i += 163
    end
    while B - i >= 27
        append!(z, int_to_base(bits_to_int(view(μ, o+i+1:o+i+27)), 2)); i += 27
    end
    isodd(k) && push!(z, 0)
    i == B ? z : nothing
end

"""Algorithm 8, the inverse of `encode`: z ∈ ℤ_q^k back to its B bits. `nothing` (⊥) unless
B = max_bits(k), every element and block is in range, and for odd k the last element is 0."""
function decode(z::AbstractVector{<:Integer}, B::Integer)
    k = length(z)
    B == max_bits(k) || return nothing
    all(x -> 0 <= x < q, z) || return nothing
    μ = Bool[]; j = firstindex(z) - 1; r = B
    while r >= 163
        x = base_to_int(view(z, j+1:j+12)); x >= big(2)^163 && return nothing
        append!(μ, int_to_bits(x, 163)); j += 12; r -= 163
    end
    while r >= 27
        x = base_to_int(view(z, j+1:j+2)); x >= big(2)^27 && return nothing
        append!(μ, int_to_bits(x, 27)); j += 2; r -= 27
    end
    isodd(k) && !iszero(z[end]) ? nothing : μ
end

end # module
