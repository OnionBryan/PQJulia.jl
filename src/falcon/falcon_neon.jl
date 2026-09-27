# src/falcon/falcon_neon.jl
# ============================================================================
# Falcon — NEON acceleration: the native (C/NEON) integer NTT over q=12289,
# exposed via libforgeddec (forged_ntt_*). Exact uint32 arithmetic (no f32 loss),
# so it is a drop-in accelerator for Falcon's mod-q polynomial multiply / h=g/f.
# The native NTT is cyclic (xⁿ−1); negacyclic (xⁿ+1) is the ψ-weighting here,
# with ψ = g^((p-1)/2n) so ψ² equals the native ω (= g^((p-1)/n)).
# ============================================================================
module FalconNEON

const LIB = abspath(joinpath(@__DIR__, "..", "..", "native", "libforgeddec.dylib"))
const P = UInt32(12289)
const ROOT = UInt32(11)                  # a primitive root of 12289
isavailable() = isfile(LIB)

twiddles(n) = (tw = Vector{UInt32}(undef, n); twi = Vector{UInt32}(undef, n);
    ccall((:forged_ntt_twiddles, LIB), Cvoid, (UInt32, UInt32, Csize_t, Ptr{UInt32}, Ptr{UInt32}),
          ROOT, P, n, tw, twi); (tw, twi))
ntt_forward!(a, tw) = ccall((:forged_ntt_forward, LIB), Cvoid,
    (Ptr{UInt32}, Ptr{UInt32}, UInt32, Csize_t), a, tw, P, length(a))
ntt_inverse!(a, twi, ni) = ccall((:forged_ntt_inverse, LIB), Cvoid,
    (Ptr{UInt32}, Ptr{UInt32}, UInt32, UInt32, Csize_t), a, twi, P, UInt32(ni), length(a))

# Negacyclic multiply a·b mod (xⁿ+1) mod q, via the native NEON NTT.
function neon_negamul(a::AbstractVector{<:Integer}, b::AbstractVector{<:Integer})
    n = length(a); p = Int(P)
    ψ = powermod(11, (p - 1) ÷ (2n), p); ψi = invmod(ψ, p); ni = invmod(n, p)
    tw, twi = twiddles(n)
    â = UInt32[mod(Int(mod(a[i+1], p)) * powermod(ψ, i, p), p) for i in 0:n-1]
    b̂ = UInt32[mod(Int(mod(b[i+1], p)) * powermod(ψ, i, p), p) for i in 0:n-1]
    ntt_forward!(â, tw); ntt_forward!(b̂, tw)
    Ĉ = UInt32[UInt32(mod(UInt64(â[i]) * UInt64(b̂[i]), p)) for i in 1:n]
    ntt_inverse!(Ĉ, twi, ni)
    return [Int(mod(Int(Ĉ[i+1]) * powermod(ψi, i, p), p)) for i in 0:n-1]
end

end # module
