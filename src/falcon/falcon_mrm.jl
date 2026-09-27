# FALCON-MRM, ePrint 2026/420. Not fixed by the spec: TAG_H1/TAG_H2, header 0x70 + log₂n,
# the check 0 ≤ M1ᵢ < q.
module FalconMRM

import ..Falcon
import ..FalconChaCha as CC
import ..FalconEncoding as FE

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
# HashToPoint with the input prefixed by a tag (§1.1).
H1(ρ, M1, M2, λ) = Falcon.hash_to_point(vcat(pack14(vcat(ρ, M1)), M2), TAG_H1, λ)
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

"""
    sign(M1, M2, ek::Falcon.ExpandedKey; randombytes) -> Vector{UInt8}

Algorithm 3. A fresh ρ per attempt; returns header ‖ Compress(s1) ‖ Compress(s2), where
s1 + s2·h = c = (H1(ρ, M), H2(c1) + (ρ, M1)).
"""
function sign(M1::AbstractVector{<:Integer}, M2::AbstractVector{UInt8}, ek::Falcon.ExpandedKey;
              randombytes=Falcon.sysrandom)
    n = ek.sk.n; (; λ, γ) = params(n); fp = Falcon.params(n); L = slen(n) * 8
    length(M1) == n - λ - γ || throw(ArgumentError("M1 must have $(n - λ - γ) elements of ℤ_q"))
    all(x -> 0 <= x < q, M1) || throw(ArgumentError("M1 elements must lie in [0, q)"))
    while true
        ρ = uniform_zq(γ, randombytes)
        c1 = H1(ρ, M1, M2, λ)
        c = vcat(c1, mod.(H2(c1, n - λ) .+ vcat(ρ, M1), q))
        t0, t1 = Falcon.target(ek.gs, c)
        s1, s2 = Falcon.preimage(ek.gs, c, t0, t1, fp.σmin, CC.ChaCha20(randombytes(Falcon.SEED_LEN)))
        sum(abs2, s1) + sum(abs2, s2) <= fp.β2 || continue
        e1 = FE.compress(s1, L); e1 === nothing && continue
        e2 = FE.compress(s2, L); e2 === nothing && continue
        return vcat(header(n), e1, e2)
    end
end

"""
    verify(M2, sig, h) -> Union{Vector{Int}, Nothing}

Algorithm 4. Returns the recovered M1, or `nothing` (⊥).
"""
function verify(M2::AbstractVector{UInt8}, sig::AbstractVector{UInt8}, h::AbstractVector{<:Integer})
    n = length(h); haskey(PARAMS, n) || return nothing
    (; λ, γ) = params(n); fp = Falcon.params(n); ℓ = slen(n)
    (length(sig) == sig_bytes(n) && sig[1] == header(n)) || return nothing
    s1 = FE.decompress(@view(sig[2:1+ℓ]), 8ℓ, n); s1 === nothing && return nothing
    s2 = FE.decompress(@view(sig[2+ℓ:end]), 8ℓ, n); s2 === nothing && return nothing
    sum(abs2, s1) + sum(abs2, s2) <= fp.β2 || return nothing
    c = mod.(s1 .+ Falcon.poly_mul_modq(s2, h), q)
    c1 = c[1:λ]
    ρM1 = mod.(c[λ+1:end] .- H2(c1, n - λ), q)
    ρ, M1 = ρM1[1:γ], ρM1[γ+1:end]
    H1(ρ, M1, M2, λ) == c1 ? M1 : nothing
end

# ── §1.4 bit strings ↔ ℤ_q^k: 163 bits per 12 elements, then 27 bits per 2 ──
"B(k): bits encodable into ℤ_q^k by `encode` (6492 for k = 478, 13067 for k = 962)."
max_bits(k) = 163 * (k ÷ 12) + 27 * ((k % 12) ÷ 2)

bits_to_int(b) = foldl((x, bit) -> (x << 1) | big(bit), b; init = big(0))
int_to_bits(x, ℓ) = [isodd(x >> (ℓ - 1 - i)) for i in 0:ℓ-1]
int_to_base(x, k) = (z = Int[]; for _ in 1:k; push!(z, Int(mod(x, q))); x = div(x, q); end; z)
base_to_int(z) = foldr((zi, x) -> x * q + zi, z; init = big(0))

"Algorithm 7: encode exactly `max_bits(k)` bits into ℤ_q^k; `nothing` on a length mismatch."
function encode(μ::AbstractVector{Bool}, k::Int)
    B = length(μ); B == max_bits(k) || return nothing
    z = Int[]; i = 0
    while B - i >= 163
        append!(z, int_to_base(bits_to_int(μ[i+1:i+163]), 12)); i += 163
    end
    while B - i >= 27
        append!(z, int_to_base(bits_to_int(μ[i+1:i+27]), 2)); i += 27
    end
    i == B ? z : nothing
end

"Algorithm 8: decode z ∈ ℤ_q^k into B bits; `nothing` (⊥) on an out-of-range block."
function decode(z::AbstractVector{<:Integer}, B::Int)
    all(x -> 0 <= x < q, z) || return nothing
    μ = Bool[]; j = 0; r = B
    while r >= 163
        j + 12 <= length(z) || return nothing
        x = base_to_int(z[j+1:j+12]); x >= big(2)^163 && return nothing
        append!(μ, int_to_bits(x, 163)); j += 12; r -= 163
    end
    while r >= 27
        j + 2 <= length(z) || return nothing
        x = base_to_int(z[j+1:j+2]); x >= big(2)^27 && return nothing
        append!(μ, int_to_bits(x, 27)); j += 2; r -= 27
    end
    r == 0 ? μ : nothing
end

end # module
