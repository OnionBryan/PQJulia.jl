# PQJulia.jl

[![CI](https://github.com/OnionBryan/PQJulia.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/OnionBryan/PQJulia.jl/actions/workflows/ci.yml)

Post-quantum cryptography in pure Julia: the NIST lattice standards ML-KEM (FIPS 203) and
ML-DSA (FIPS 204), and Falcon, the NTRU-lattice signature NIST is standardizing as FN-DSA (FIPS 206).

| Standard | Algorithm | Parameter sets |
|----------|-----------|----------------|
| FIPS 203 | ML-KEM (Kyber) | 512 / 768 / 1024 |
| FIPS 204 | ML-DSA (Dilithium) | 44 / 65 / 87 — pure, HashML-DSA, internal and external-μ interfaces |
| FN-DSA (draft FIPS 206) | Falcon (round-3 spec v1.2) | Falcon-512 / Falcon-1024 |
| RFC 7748 | X25519 | Curve25519 Diffie–Hellman |
| draft-connolly-cfrg-xwing-kem-11 | X-Wing | X25519 + ML-KEM-768 hybrid KEM |

Also includes Shamir (k,n)-threshold secret sharing over GF(2^521 − 1).

## Validation

Everything below runs in CI (`Pkg.test()`):

- **ML-KEM:** 240 NIST ACVP vectors: KeyGen, Encaps, Decaps, and the FIPS 203 encapsulation-
  and decapsulation-key input checks, at all three levels.
- **ML-DSA:** 615 NIST ACVP vectors: KeyGen; SigGen deterministic and hedged, pure and
  pre-hash, internal and external-μ; SigVer for every interface, at all three levels.
- **Falcon:** byte-exact against the round-3 C implementation. From the same key and
  randomness stream, all 120 reference signatures (n = 2 … 1024) are reproduced bit for bit,
  and SamplerZ matches all 3,072 reference vectors. The vectors come from
  [tprest/falcon.py](https://github.com/tprest/falcon.py) (MIT).
- **Falcon key certificate:** the exact ffLDL leaves and Gram–Schmidt norm are computed with no
  floating point. The multi-modular result matches an independent ℚ(x) tower, the determinant
  identity ∏d = qⁿ holds, and all 120 reference keys certify. A deliberately oversized key is
  rejected by both the certificate and the signer.
- **Falcon extensions:**
  - **Float-free keygen.** The fixed-point FFT reproduces ntrugen's root table exactly and
    matches the exact negacyclic product. The fixed-point GS check agrees with the exact GS
    norm, and every generated key passes the NTRU equation and the certificate.
  - **Message recovery (FALCON-MRM).** Tested for round trips, tampering and the sizes of
    ePrint 2026/420.
- **Wycheproof** (C2SP): 2,863 ML-DSA and ML-KEM edge-case vectors, including malformed keys
  and hints, out-of-range secret keys, boundary norms, and modulus-overflow encapsulation keys.
- **X25519:** the RFC 7748 vectors (including 1,000 iterations; 1,000,000 with
  `PQJULIA_LONG_TESTS=1`) and all 518 Wycheproof vectors (twist, low-order and non-canonical
  inputs). It also matches 164 vectors from an independent engine, the complete twisted
  Edwards law in paper5-haskell (`Edwards.hs`), whose prime is proved by a Pratt certificate.
  The 51-bit-limb engine agrees with a BigInt transcription of RFC 7748 §5 on random inputs.
- **X-Wing:** the three Appendix C vectors of draft-11, byte-exact (pk, ct, ss).
- **CCTV** (C2SP): the ML-KEM strcmp, "unlucky" XOF and 2,595 invalid-modulus vectors, and the
  accumulated ML-DSA vectors (10,000 seeded keys per level in CI, byte-exact with the reference).
- Property tests: round-trips, tamper rejection, malformed inputs.

`test/interop/liboqs_interop.jl` (opt-in; needs liboqs) checks two-way interoperability with
liboqs 0.15 for every scheme. Keys, signatures and ciphertexts produced by either side are
accepted by the other, and so are secret keys used for signing.

### Timing test

`test/timing/dudect.jl` (opt-in) runs the fixed-vs-random Welch t-test of
[dudect](https://eprint.iacr.org/2016/1069) on the functions [SECURITY.md](SECURITY.md) lists as
constant-time and on secret-key paths:
- **ML-KEM:** decapsulation of valid and tampered ciphertexts, encapsulation, noise sampling with
  the NTT, `kyber_verify` and `kyber_cmov!`.
- **ML-DSA:** secret-key unpacking, the NTT, `decompose` and `make_hint`.
- **X25519:** the ladder and its output encoding.

|t| > 4.5 on any percentile crop flags a timing difference. Run it on an idle machine with
`julia --project=. test/timing/dudect.jl` (`DUDECT_N` sets the sample count, default 100,000).
Fixed inputs are fresh copies, and the `decompose` and `make_hint` classes are both random,
differing only in the branch a leak would follow: repeating identical data runs measurably faster
on its own, on both Apple and Intel CPUs.

`julia test/timing/dispatch.jl` (runs in CI; installs [JET](https://github.com/aviatesk/JET.jl)
into a temporary environment) checks key generation, encapsulation, decapsulation and signing for
ML-KEM and ML-DSA at every level, Falcon key expansion, X25519 and X-Wing for runtime dispatch. A
boxed or type-unstable value on secret data makes timing depend on the value, which is how the
X25519 encoding leaked; JET finds it without timing noise.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/OnionBryan/PQJulia.jl")
```

Requires Julia 1.12. The only runtime dependencies are the standard libraries `SHA`,
`Random` and `LinearAlgebra`.

## Usage

```julia
using PQJulia
msg = Vector{UInt8}("hello")

# Key encapsulation (ML-KEM-768)
pk, sk = MLKEM.Category3.kyber_kem_keypair()
ct, shared_secret = MLKEM.Category3.kyber_kem_enc(pk)
shared_secret == MLKEM.Category3.kyber_kem_dec(ct, sk)   # true

# Signatures (ML-DSA-65); hedged by default, pass hedged=false for deterministic
pk, sk = MLDSA.Category3.dilithium_keygen()
sig = MLDSA.Category3.dilithium_sign(msg, sk; context=Vector{UInt8}("app-v1"))
MLDSA.Category3.dilithium_verify(msg, sig, pk; context=Vector{UInt8}("app-v1"))   # true

# Pre-hash signatures (HashML-DSA)
sig = MLDSA.Category3.dilithium_sign_prehash(msg, sk, "SHA2-512")
MLDSA.Category3.dilithium_verify_prehash(msg, sig, pk, "SHA2-512")   # true

# Signatures (Falcon-512)
pk, sk = FNDSA.Falcon512.falcon_keygen()
sig = FNDSA.Falcon512.falcon_sign(msg, sk)
FNDSA.Falcon512.falcon_verify(msg, sig, pk)   # true
ek = FNDSA.Falcon512.falcon_expand_sk(sk)     # decode once, then sign many messages fast
sig = FNDSA.Falcon512.falcon_sign(msg, ek)
pk, sk = FNDSA.Falcon512.falcon_keygen(certified=true)   # only keys with an exact certificate
FNDSA.Falcon512.falcon_certify(sk).ok                    # true: leaves in [σmin, σmax], GS norm ≤ 1.17√q
pk, sk = FNDSA.Falcon512.falcon_keygen(fixedpoint=true)  # no floating point in keygen (ePrint 2023/290)
pk, sk = FNDSA.Falcon512.falcon_keygen(fips206=true)     # FIPS 206 GS bound 0.9999·1.17√q

# Message recovery (FALCON-MRM, ePrint 2026/420): M1 travels inside the signature
using PQJulia.FNDSA: FalconMRM
bits = rand(Bool, FNDSA.Falcon512.MRM_MAX_BITS)              # 6492 bits
M1 = FalconMRM.encode(bits, FNDSA.Falcon512.MRM_M1_LEN)
sig = FNDSA.Falcon512.falcon_mrm_sign(M1, UInt8[], sk)        # 1251 bytes
FalconMRM.decode(FNDSA.Falcon512.falcon_mrm_verify(UInt8[], sig, pk), length(bits)) == bits   # true

# Hybrid key encapsulation (X-Wing: X25519 + ML-KEM-768)
pk, sk = XWing.xwing_keypair()             # 1216-byte pk, 32-byte sk
ct, ss = XWing.xwing_encaps(pk)            # 1120-byte ct, 32-byte shared secret
ss == XWing.xwing_decaps(ct, sk)           # true

# Secret sharing
shares = shamir_share_bytes(rand(UInt8, 32), 3, 5)
secret = shamir_reconstruct_bytes(shares[[1, 3, 5]], 3, 32)
```

| Level | ML-KEM | ML-DSA | FN-DSA |
|-------|--------|--------|--------|
| 1 | `MLKEM.Category1` (512) | `MLDSA.Category2` (44) | `FNDSA.Falcon512` |
| 3 | `MLKEM.Category3` (768) | `MLDSA.Category3` (65) | |
| 5 | `MLKEM.Category5` (1024) | `MLDSA.Category5` (87) | `FNDSA.Falcon1024` |

Key and signature sizes are available as constants on each module, e.g.
`MLDSA.Category3.SIG_BYTES` and `FNDSA.Falcon512.PK_BYTES`. Deterministic `_derand`
variants take their randomness explicitly, for testing against known-answer vectors.

## Security

[SECURITY.md](SECURITY.md) covers randomness, timing behavior and fixed issues.

## Repository layout

- `src/`: the package.
- `test/`: the test suite and the vendored ACVP and Falcon vectors (`test/kat/`).
- `benchmarks/bench.jl`: BenchmarkTools suite with its own environment
  (`julia benchmarks/bench.jl`).

## References

- [FIPS 203 — ML-KEM](https://csrc.nist.gov/pubs/fips/203/final)
- [FIPS 204 — ML-DSA](https://csrc.nist.gov/pubs/fips/204/final)
- [Falcon specification v1.2](https://falcon-sign.info/falcon.pdf)
- [pq-crystals reference implementations](https://github.com/pq-crystals)
- [RFC 7748 — Elliptic Curves for Security (X25519)](https://www.rfc-editor.org/rfc/rfc7748)
- [X-Wing: general-purpose hybrid post-quantum KEM, draft-connolly-cfrg-xwing-kem-11](https://datatracker.ietf.org/doc/draft-connolly-cfrg-xwing-kem/)
- [tprest/falcon.py](https://github.com/tprest/falcon.py)
- Ducas & Prest, [Fast Fourier Orthogonalization](https://eprint.iacr.org/2015/1014) (ffLDL leaves as Gram–Schmidt norms)
- Pornin, [Improved Key Pair Generation for Falcon, BAT and Hawk](https://eprint.iacr.org/2023/290) and [ntrugen](https://github.com/pornin/ntrugen) (multi-modular arithmetic, fixed-point keygen)
- Günther, Lyubashevsky, Schmidt, [FALCON with message recovery](https://eprint.iacr.org/2026/420)
- Perlner, [FIPS 206 Status Update](https://csrc.nist.gov/presentations/2025/fips-206-fn-dsa-falcon) (NIST, 2025)

## License

MIT
