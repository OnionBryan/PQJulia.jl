"""
ML-KEM (FIPS 203) — Multi-level implementation via @eval module generation.
Generates MLKEM.Category1 (ML-KEM-512), MLKEM.Category3 (ML-KEM-768), MLKEM.Category5 (ML-KEM-1024).

Architecture follows NistyPQC.jl pattern:
- Shared primitives (NTT, reduce, basemul, CBD, compress) included once
- Per-level constants + KEM functions compiled 3x with different constants
"""
module MLKEM

using SHA
using Random

# Shared primitives — same for all security levels
include("kyber_core.jl")

# Level-specific parameter tuples: (category_number, K, ETA1, ETA2, DU, DV)
# Derived from pq-crystals/kyber/ref/params.h
const CATEGORY_PARAMS = [
    :Category1 => (1, 2, 3, 2, 10, 4),   # ML-KEM-512
    :Category3 => (3, 3, 2, 2, 10, 4),   # ML-KEM-768
    :Category5 => (5, 4, 2, 2, 11, 5),   # ML-KEM-1024
]

function derived_sizes(K, DU, DV)
    polyvecbytes = K * KYBER_POLYBYTES
    polycompressedbytes = DV == 4 ? 128 : 160  # 128 for d=4, 160 for d=5
    polyveccompressedbytes = K * (DU == 10 ? 320 : 352)  # 320 for d=10, 352 for d=11
    indcpa_pk = polyvecbytes + KYBER_SYMBYTES
    indcpa_sk = polyvecbytes
    indcpa_ct = polyveccompressedbytes + polycompressedbytes
    kem_pk = indcpa_pk
    kem_sk = indcpa_sk + indcpa_pk + 2 * KYBER_SYMBYTES
    kem_ct = indcpa_ct
    return (; polyvecbytes, polycompressedbytes, polyveccompressedbytes,
              indcpa_pk, indcpa_sk, indcpa_ct, kem_pk, kem_sk, kem_ct)
end

for (category, params) in CATEGORY_PARAMS
    @eval module $category

    export kyber_kem_keypair, kyber_kem_keypair_derand
    export kyber_kem_enc, kyber_kem_enc_derand
    export kyber_kem_dec

    # Import shared primitives from parent MLKEM module (defined in kyber_core.jl)
    import ..MLKEM: KYBER_N, KYBER_Q, KYBER_SYMBYTES, KYBER_POLYBYTES, KYBER_ZETAS, KYBER_MONT, KYBER_QINV
    import ..MLKEM: montgomery_reduce, barrett_reduce, fqmul, caddq
    import ..MLKEM: kyber_ntt!, kyber_invntt!
    import ..MLKEM: kyber_poly_basemul_montgomery!, kyber_poly_add!, kyber_poly_sub!
    import ..MLKEM: kyber_poly_reduce!, kyber_poly_tomont!
    import ..MLKEM: kyber_poly_tobytes!, kyber_poly_frombytes!
    import ..MLKEM: kyber_poly_frommsg!, kyber_poly_tomsg!
    import ..MLKEM: kyber_poly_compress!, kyber_poly_decompress!
    import ..MLKEM: kyber_poly_getnoise_eta1!, kyber_poly_getnoise_eta2!
    import ..MLKEM: kyber_compress, kyber_decompress
    import ..MLKEM: kyber_hash_h, kyber_hash_g
    import ..MLKEM: kyber_prf, kyber_xof
    import ..MLKEM: kyber_sample_uniform!, kyber_rej_uniform!
    import ..MLKEM: kyber_cbd2!, kyber_cbd3!
    import ..MLKEM: kyber_verify, kyber_cmov!
    import SHA
    import ..MLKEM: derived_sizes

    using Random

    # Per-level constants from parameter tuple
    const (CATEGORY_NUMBER, KYBER_K, KYBER_ETA1, KYBER_ETA2, KYBER_DU, KYBER_DV) = $params

    # Derived sizes
    const SIZES = derived_sizes(KYBER_K, KYBER_DU, KYBER_DV)
    const KYBER_POLYVECBYTES = SIZES.polyvecbytes
    const KYBER_POLYCOMPRESSEDBYTES = SIZES.polycompressedbytes
    const KYBER_POLYVECCOMPRESSEDBYTES = SIZES.polyveccompressedbytes
    const KYBER_INDCPA_PUBLICKEYBYTES = SIZES.indcpa_pk
    const KYBER_INDCPA_SECRETKEYBYTES = SIZES.indcpa_sk
    const KYBER_INDCPA_BYTES = SIZES.indcpa_ct
    const KYBER_PUBLICKEYBYTES = SIZES.kem_pk
    const KYBER_SECRETKEYBYTES = SIZES.kem_sk
    const KYBER_CIPHERTEXTBYTES = SIZES.kem_ct
    const KYBER_SSBYTES = 32
    const KYBER_INDCPA_MSGBYTES = KYBER_SYMBYTES

    const IDENTIFIER = "ML-KEM-$(KYBER_K * KYBER_N)"

    include("kyber_kem.jl")

    end # module $category
end

end # module MLKEM
