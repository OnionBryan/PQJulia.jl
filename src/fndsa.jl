"""
Falcon / FN-DSA — NTRU-lattice hash-and-sign signatures (Falcon round-3 spec v1.2,
the basis of FN-DSA, to be standardized as NIST FIPS 206). Generates FNDSA.Falcon512 and FNDSA.Falcon1024.
The math-level implementation for every degree n = 2..1024 lives in FNDSA.Falcon.
"""
module FNDSA

include("falcon/falcon_fft.jl")
include("falcon/falcon_chacha.jl")
include("falcon/falcon_sampler.jl")
include("falcon/falcon_fpr.jl")
include("falcon/falcon_ntrugen.jl")
include("falcon/falcon_encoding.jl")
include("falcon/falcon_certify.jl")
include("falcon/falcon_fxp.jl")
include("falcon/falcon.jl")
include("falcon/falcon_mrm.jl")

import ..Wipe: wipe!

# The Falcon decoders index byte strings from 1. `onebased(x)` is `x`, or a 1-based copy when `x`
# has other axes (an OffsetArray, a view through Base.IdentityUnitRange). `with_onebased(f, sk)`
# hands the decoders a secret key as a plain Vector: any other AbstractVector (other axes, a view,
# a key store whose reads can fail) is copied inside the `try`, so the copy is wiped when `f`
# returns or throws and also when reading the key fails part-way.
onebased(x::AbstractVector{UInt8}) = Base.has_offset_axes(x) ? copyto!(Vector{UInt8}(undef, length(x)), x) : x
function with_onebased(f, x::AbstractVector{UInt8})
    x isa Vector{UInt8} && return f(x)
    y = Vector{UInt8}(undef, length(x))
    try
        copyto!(y, x)
        f(y)
    finally
        wipe!(y)
    end
end

for (name, n) in [(:Falcon512, 512), (:Falcon1024, 1024)]
    @eval module $name

    export falcon_keygen, falcon_sign, falcon_verify, falcon_expand_sk, falcon_wipe!, falcon_certify,
           falcon_mrm_sign, falcon_mrm_verify

    import ..Falcon, ..FalconMRM, ..onebased, ..with_onebased

    const N = $n
    const PK_BYTES = Falcon.pk_bytes(N)
    const SK_BYTES = Falcon.sk_bytes(N)
    const SIG_BYTES = Falcon.sig_bytes(N)
    const IDENTIFIER = "Falcon-$N"

    """Generate a keypair; returns `(pk, sk)` in the spec byte encodings. With `certified=true`, only
    keys that pass the exact certificate (`falcon_certify`) are returned; `fips206=true` also tightens
    the GS bound to 0.9999·1.17√q (NIST FIPS 206 status update). `fixedpoint=true` generates
    the key with no floating point (ePrint 2023/290)."""
    falcon_keygen(; certified::Bool=false, fips206::Bool=false, fixedpoint::Bool=false) =
        Falcon.keypair(N; certified, fips206, fixedpoint)

    """Exact certificate for `sk`: every ffLDL leaf σ/√d in [σmin, σmax] and Gram–Schmidt norm ≤ 1.17√q,
    decided without floating point (FalconCertify). Returns a NamedTuple with `ok` and the exact values."""
    function falcon_certify(sk::AbstractVector{UInt8}; fips206::Bool=false)
        length(sk) == SK_BYTES || throw(ArgumentError("$IDENTIFIER secret key must be $SK_BYTES bytes"))
        with_onebased(sk) do s
            k = Falcon.decode_sk(s); k === nothing && throw(ArgumentError("invalid $IDENTIFIER secret key"))
            try
                Falcon.certify(k; fips206)
            finally
                Falcon.wipe_sk!(k)                       # also when certify throws
            end
        end
    end

    "Decode `sk` once and precompute its ffLDL tree, for repeated signing."
    function falcon_expand_sk(sk::AbstractVector{UInt8})
        length(sk) == SK_BYTES || throw(ArgumentError("$IDENTIFIER secret key must be $SK_BYTES bytes"))
        with_onebased(Falcon.expand_sk, sk)
    end

    "Zero an expanded key's secret arrays (decoded key, FFT basis, ffLDL tree); it is unusable afterwards."
    falcon_wipe!(ek::Falcon.ExpandedKey) = Falcon.wipe!(ek)

    "Sign `msg`; returns a padded $(Falcon.sig_bytes($n))-byte signature."
    function falcon_sign(msg::AbstractVector{UInt8}, sk::AbstractVector{UInt8})
        ek = falcon_expand_sk(sk)
        try
            Falcon.sign(onebased(msg), ek)
        finally
            Falcon.wipe!(ek)                             # also when signing throws
        end
    end
    falcon_sign(msg::AbstractVector{UInt8}, ek::Falcon.ExpandedKey) =
        ek.sk.n == N ? Falcon.sign(onebased(msg), ek) : throw(ArgumentError("key is not a $IDENTIFIER key"))

    "Return `true` iff `sig` is a valid signature of `msg` under `pk`."
    falcon_verify(msg::AbstractVector{UInt8}, sig::AbstractVector{UInt8}, pk::AbstractVector{UInt8}) =
        length(pk) == PK_BYTES && Falcon.verify(onebased(msg), onebased(sig), onebased(pk))

    # FALCON-MRM (ePrint 2026/420): M1 ∈ ℤ_q^MRM_M1_LEN is recovered from the signature, M2 is sent in clear.
    const MRM_M1_LEN = FalconMRM.m1_len(N)
    const MRM_SIG_BYTES = FalconMRM.sig_bytes(N)
    const MRM_MAX_BITS = FalconMRM.max_bits(MRM_M1_LEN)

    "Sign (M1, M2) with message recovery; returns an $(FalconMRM.sig_bytes($n))-byte signature."
    function falcon_mrm_sign(M1::AbstractVector{<:Integer}, M2::AbstractVector{UInt8}, sk::AbstractVector{UInt8})
        FalconMRM.check_m1(M1, N)                        # before the key is decoded
        ek = falcon_expand_sk(sk)
        try
            FalconMRM.sign(M1, M2, ek)
        finally
            Falcon.wipe!(ek)                             # also when signing throws
        end
    end
    falcon_mrm_sign(M1::AbstractVector{<:Integer}, M2::AbstractVector{UInt8}, ek::Falcon.ExpandedKey) =
        ek.sk.n == N ? FalconMRM.sign(M1, M2, ek) : throw(ArgumentError("key is not a $IDENTIFIER key"))

    "Verify with message recovery; returns M1, or `nothing` if the signature is invalid."
    function falcon_mrm_verify(M2::AbstractVector{UInt8}, sig::AbstractVector{UInt8}, pk::AbstractVector{UInt8})
        h = length(pk) == PK_BYTES ? Falcon.decode_pk(onebased(pk), N) : nothing
        h === nothing ? nothing : FalconMRM.verify(M2, sig, h)
    end

    end # module $name
end

end # module FNDSA
