"""
ML-DSA (FIPS 204) — Multi-level implementation via @eval module generation.
Generates MLDSA.Category2 (ML-DSA-44), MLDSA.Category3 (ML-DSA-65), MLDSA.Category5 (ML-DSA-87).
"""
module MLDSA

using SHA
using Random

# Shared primitives — same for all security levels
include("dilithium_core.jl")

# Level parameters: (cat_num, K, L, TAU, ETA, OMEGA, GAMMA2_DIV, LG_GAMMA1, CTILDEBYTES)
# GAMMA2_DIV: 32 means GAMMA2=(Q-1)/32, 88 means GAMMA2=(Q-1)/88
const CATEGORY_PARAMS = [
    :Category2 => (2, 4, 4, 39, 2, 80, 88, 17, 32),   # ML-DSA-44
    :Category3 => (3, 6, 5, 49, 4, 55, 32, 19, 48),   # ML-DSA-65
    :Category5 => (5, 8, 7, 60, 2, 75, 32, 19, 64),   # ML-DSA-87
]

function derived_sizes(K, L, ETA, OMEGA, GAMMA2_DIV, LG_GAMMA1, CTILDEBYTES)
    polyeta_packed = ETA == 2 ? 96 : 128  # 3-bit: 96, 4-bit: 128
    polyz_packed = LG_GAMMA1 == 17 ? 576 : 640  # 18-bit: 576, 20-bit: 640
    polyw1_packed = GAMMA2_DIV == 88 ? 192 : 128  # 6-bit: 192, 4-bit: 128
    polyt1_packed = 320  # always 10-bit
    polyt0_packed = 416  # always 13-bit
    pk_bytes = 32 + K * polyt1_packed
    sk_bytes = 2*32 + 64 + L*polyeta_packed + K*polyeta_packed + K*polyt0_packed
    sig_bytes = CTILDEBYTES + L*polyz_packed + OMEGA + K
    return (; polyeta_packed, polyz_packed, polyw1_packed, polyt1_packed, polyt0_packed,
              pk_bytes, sk_bytes, sig_bytes)
end

for (category, params) in CATEGORY_PARAMS
    @eval module $category

    export dilithium_keygen, dilithium_keygen_derand
    export dilithium_sign, dilithium_sign_derand
    export dilithium_verify

    # Import shared primitives from parent
    import ..MLDSA: Q, N, D, ZETAS
    import ..MLDSA: montgomery_reduce, reduce32, caddq, freeze
    import ..MLDSA: ntt!, invntt!
    import ..MLDSA: poly_pointwise!, poly_add!, poly_sub!
    import ..MLDSA: poly_reduce!, poly_caddq!, poly_shiftl!, poly_chknorm
    import ..MLDSA: poly_uniform!, power2round
    import ..MLDSA: derived_sizes

    import SHA
    using Random

    # Per-level constants
    const (CAT_NUM, K, L, TAU, ETA, OMEGA, _GAMMA2_DIV, _LG_GAMMA1, CTILDEBYTES) = $params
    const BETA = Int32(TAU * ETA)
    const GAMMA1 = Int32(1 << _LG_GAMMA1)
    const GAMMA2 = Int32(div(Q - 1, _GAMMA2_DIV))
    const SEEDBYTES = 32
    const CRHBYTES = 64
    const TRBYTES = 64

    const SIZES = derived_sizes(K, L, ETA, OMEGA, _GAMMA2_DIV, _LG_GAMMA1, CTILDEBYTES)
    const POLYETA_PACKED = SIZES.polyeta_packed
    const POLYZ_PACKED = SIZES.polyz_packed
    const POLYW1_PACKED = SIZES.polyw1_packed
    const POLYT1_PACKED = SIZES.polyt1_packed
    const POLYT0_PACKED = SIZES.polyt0_packed
    const PK_BYTES = SIZES.pk_bytes
    const SK_BYTES = SIZES.sk_bytes
    const SIG_BYTES = SIZES.sig_bytes

    const IDENTIFIER = "ML-DSA-$K$L"

    # Include parameterized level-specific functions
    include("dilithium_level.jl")

    end # module $category
end

end # module MLDSA
