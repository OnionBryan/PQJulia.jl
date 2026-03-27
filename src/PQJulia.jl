"""
    PQJulia.jl — Post-Quantum Cryptography for Julia

Complete implementations of NIST FIPS post-quantum cryptographic standards:
- **ML-KEM** (FIPS 203): Module-Lattice Key Encapsulation (Kyber)
- **ML-DSA** (FIPS 204): Module-Lattice Digital Signatures (Dilithium)
- **SLH-DSA** (FIPS 205): Stateless Hash-Based Digital Signatures (SPHINCS+)
- **Shamir**: (k,n)-threshold secret sharing over GF(2^127-1)

All ML-KEM and ML-DSA implementations pass NIST ACVP Known Answer Tests (750+ vectors).
ML-KEM is C cross-validated (byte-exact with pq-crystals/kyber reference).

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
sig = MLDSA.Category3.dilithium_sign(msg, sk)
MLDSA.Category3.dilithium_verify(msg, sig, pk)  # true

# SLH-DSA-SHAKE-128f Hash-Based Signatures
pk, sk = SLHDSA.slh_keygen(SLHDSA.SHAKE_128F)
sig = SLHDSA.slh_sign(msg, sk, SLHDSA.SHAKE_128F)
SLHDSA.slh_verify(msg, sig, pk, SLHDSA.SHAKE_128F)  # true
```
"""
module PQJulia

# ML-KEM (FIPS 203) — Kyber Key Encapsulation
include("mlkem.jl")
using .MLKEM

# ML-DSA (FIPS 204) — Dilithium Digital Signatures
include("mldsa.jl")
using .MLDSA

# Re-export modules
export MLKEM, MLDSA

# Shamir Secret Sharing
include("shamir.jl")
export shamir_share, shamir_reconstruct

end # module PQJulia
