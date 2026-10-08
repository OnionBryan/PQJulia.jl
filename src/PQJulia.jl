"""
    PQJulia.jl — Post-Quantum Cryptography for Julia

Implementations of NIST FIPS post-quantum cryptographic standards:
- **ML-KEM** (FIPS 203): Module-Lattice Key Encapsulation (Kyber) — 512/768/1024
- **ML-DSA** (FIPS 204): Module-Lattice Digital Signatures (Dilithium) — 44/65/87
- **FN-DSA** (Falcon): NTRU-lattice signatures — Falcon-512 / Falcon-1024; exact key certificate,
  float-free keygen, message recovery (FALCON-MRM)
- **X25519** (RFC 7748) and **X-Wing** (X25519 + ML-KEM-768 hybrid KEM)
- **Shamir**: (k,n)-threshold secret sharing over GF(2^521-1)

ML-KEM and ML-DSA pass 855 NIST ACVP vectors across every interface; Falcon signing is
byte-exact with the round-3 C reference implementation.

## Quick Start

```julia
using PQJulia

# ML-KEM-768 Key Encapsulation
pk, sk = MLKEM.Category3.kyber_kem_keypair()
ct, ss_enc = MLKEM.Category3.kyber_kem_enc(pk)
ss_dec = MLKEM.Category3.kyber_kem_dec(ct, sk)
@assert ss_enc == ss_dec

# ML-DSA-65 Digital Signatures
pk, sk = MLDSA.Category3.dilithium_keygen()
msg = Vector{UInt8}("hello")
sig = MLDSA.Category3.dilithium_sign(msg, sk)
MLDSA.Category3.dilithium_verify(msg, sig, pk)  # true

# Falcon-512 Signatures
pk, sk = FNDSA.Falcon512.falcon_keygen()
sig = FNDSA.Falcon512.falcon_sign(msg, sk)
FNDSA.Falcon512.falcon_verify(msg, sig, pk)  # true
```
"""
module PQJulia

include("wipe.jl")
using .Wipe
export wipe!

# Keccak-f[1600]: SHAKE128/256 and SHA3-256/512 for every scheme
include("keccak.jl")

# ML-KEM (FIPS 203) — Kyber Key Encapsulation
include("mlkem.jl")
using .MLKEM

# ML-DSA (FIPS 204) — Dilithium Digital Signatures
include("mldsa.jl")
using .MLDSA

# FN-DSA (Falcon) — NTRU-lattice hash-and-sign signatures
include("fndsa.jl")
using .FNDSA

# X25519 (RFC 7748) and the X-Wing hybrid KEM (X25519 + ML-KEM-768)
include("x25519.jl")
using .X25519
include("xwing.jl")
using .XWing

# Re-export modules
export MLKEM, MLDSA, FNDSA, X25519, XWing

# Shamir Secret Sharing
include("shamir.jl")
export shamir_share, shamir_reconstruct, shamir_share_bytes, shamir_reconstruct_bytes

end # module PQJulia
