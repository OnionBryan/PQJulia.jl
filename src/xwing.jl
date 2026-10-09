"""
X-Wing: the X25519 + ML-KEM-768 hybrid KEM of draft-connolly-cfrg-xwing-kem-11 (§5), an
individual Internet-Draft in the Independent Submission stream since -07 (2025-05-26); -11
refreshed the expired -10 with no normative change.
IND-CCA if either X25519 (gap-CDH) or ML-KEM-768 is secure, with SHA3 as a random oracle (§6;
proof in ePrint 2024/039). Against a quantum adversary, security rests on ML-KEM-768 (with
SHA3-256 as a PRF), not on X25519.
Keys and ciphertexts are the draft's fixed-length byte strings.
"""
module XWing

using SHA, Random
using ..Wipe: wipe!
import ..Keccak
import ..MLKEM.Category3 as MK
import ..X25519: x25519, x25519_base

export xwing_keypair, xwing_keypair_derand, xwing_encaps, xwing_encaps_derand, xwing_decaps

const PK_BYTES = 1216
const SK_BYTES = 32
const CT_BYTES = 1120
const SS_BYTES = 32
const LABEL = UInt8[0x5c, 0x2e, 0x2f, 0x2f, 0x5e, 0x5c]          # "\./" ‖ "/^\"

# §5.2 expandDecapsulationKey: SHAKE256(sk, 96 bytes) → ML-KEM-768 (d, z) and the X25519 scalar.
function expand(sk::AbstractVector{UInt8})
    length(sk) == SK_BYTES || throw(ArgumentError("X-Wing decapsulation key must be $SK_BYTES bytes"))
    s = collect(sk)
    e = Keccak.shake256(s, UInt64(96))
    dz = e[1:64]
    pkM, skM = MK.kyber_kem_keypair_derand(dz)
    skX = e[65:96]
    wipe!(s, e, dz)
    (; skM, skX, pkM, pkX = x25519_base(skX))
end

# §5.3
function combiner(ssM, ssX, ctX, pkX)
    inp = vcat(ssM, ssX, ctX, pkX, LABEL)
    ss = Keccak.sha3_256(inp)
    wipe!(inp)
    ss
end

"§5.2.1 GenerateKeyPairDerand: `(pk, sk)` from a 32-byte seed (sk is the seed)."
function xwing_keypair_derand(seed::AbstractVector{UInt8})
    k = expand(seed)
    wipe!(k.skM, k.skX)
    vcat(k.pkM, k.pkX), collect(seed)
end

"§5.2 GenerateKeyPair: `(pk, sk)` with a fresh seed from the OS CSPRNG."
function xwing_keypair()
    seed = rand(RandomDevice(), UInt8, SK_BYTES)
    kp = xwing_keypair_derand(seed)
    wipe!(seed)
    kp
end

"§5.4.1 EncapsulateDerand: `(ct, ss)`; eseed[1:32] is the ML-KEM message, eseed[33:64] the X25519 scalar."
function xwing_encaps_derand(pk::AbstractVector{UInt8}, eseed::AbstractVector{UInt8})
    length(pk) == PK_BYTES || throw(ArgumentError("X-Wing encapsulation key must be $PK_BYTES bytes"))
    length(eseed) == 64 || throw(ArgumentError("X-Wing eseed must be 64 bytes"))
    pkM = collect(pk[1:1184]); pkX = pk[1185:1216]
    m = eseed[1:32]; ekX = eseed[33:64]
    ctX = x25519_base(ekX)
    ssX = x25519(ekX, pkX)
    ctM, ssM = MK.kyber_kem_enc_derand(pkM, m)    # throws on a failed FIPS 203 §7.2 check
    ss = combiner(ssM, ssX, ctX, pkX)
    wipe!(m, ekX, ssX, ssM)
    vcat(ctM, ctX), ss
end

"§5.4 Encapsulate: `(ct, ss)` with fresh randomness."
function xwing_encaps(pk::AbstractVector{UInt8})
    eseed = rand(RandomDevice(), UInt8, 64)
    res = xwing_encaps_derand(pk, eseed)
    wipe!(eseed)
    res
end

"§5.5 Decapsulate: the 32-byte shared secret."
function xwing_decaps(ct::AbstractVector{UInt8}, sk::AbstractVector{UInt8})
    length(ct) == CT_BYTES || throw(ArgumentError("X-Wing ciphertext must be $CT_BYTES bytes"))
    k = expand(sk)
    ctM = collect(ct[1:1088]); ctX = ct[1089:1120]
    ssM = MK.kyber_kem_dec(ctM, k.skM)
    ssX = x25519(k.skX, ctX)
    ss = combiner(ssM, ssX, ctX, k.pkX)
    wipe!(k.skM, k.skX, ssM, ssX)
    ss
end

end # module
