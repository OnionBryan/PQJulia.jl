# Security Policy

## Scope

A pure-Julia implementation of NIST FIPS 203 (ML-KEM), FIPS 204 (ML-DSA), Falcon (round-3
specification v1.2, the basis of FN-DSA / FIPS 206), X25519 (RFC 7748) and X-Wing. ML-KEM and
ML-DSA pass 855 NIST ACVP vectors covering every interface. Falcon signing is byte-exact with
the round-3 C implementation on 120 reference vectors.

## Randomness

All secret randomness comes from the operating system CSPRNG (`Random.RandomDevice`):
ML-KEM keygen and encapsulation, ML-DSA keygen and hedged `rnd`, Falcon keygen, salts and
per-signature sampler seeds (expanded with ChaCha20, as in the reference), X25519 and X-Wing
keys, and Shamir coefficients. Julia's default task-local RNG (Xoshiro) is never used for key
material. ML-DSA signing is hedged by default, as FIPS 204 recommends; `hedged=false` selects
the deterministic variant.

## Timing

Julia compiles through LLVM, which may introduce data-dependent branches; source-level
constant-time code is not guaranteed to stay constant-time after compilation.

### Constant-time at the source level

| Function | Mechanism | File |
|----------|-----------|------|
| `kyber_verify` | OR-accumulate pairwise-XOR differences, `(-UInt64(r)) >>> 63` | kyber_core.jl |
| `kyber_cmov!` | Byte-wise mask with `UInt8(0) - b` | kyber_core.jl |
| `kyber_poly_tomsg!` | Barrett multiply-shift, no division | kyber_core.jl |
| `kyber_poly_compress!` | Barrett multiply-shift (80635>>28, 40318>>27) | kyber_core.jl |
| ML-KEM Decaps rejection | Always computes both paths, `cmov` selects | kyber_kem.jl |
| X25519 ladder | Mask swaps on 51-bit limbs | x25519.jl |

### Variable-time

| Operation | Note |
|-----------|------|
| Julia GC pauses, JIT compilation | First call compiles; warm up before timing-sensitive use |
| `decompose` (ML-DSA) | Multiply-shift on secret data (no UDIV); GAMMA2 branches are public constants. RUSTSEC-2025-0144 does not apply |
| SHAKE (SHA.jl) | Timing depends on input length, which is public |
| Falcon signing (ffSampling, SamplerZ) | Floating-point FFT and Gaussian sampling on secret data |
| X25519 scalar decode, final inversion and encoding | BigInt arithmetic |
| Falcon keygen (NTRU solver) | BigInt arithmetic; `fixedpoint=true` removes floating point (ePrint 2023/290) |

## Known Issues Addressed

| Issue | Source | Status |
|-------|--------|--------|
| XOF buffer exhaustion (poly_uniform, poly_uniform_eta, poly_challenge, kyber_sample_uniform) | C ref comparison + audit | Fixed — re-squeeze loops added |
| KyberSlash (division by q) | kyberslash.cr.yp.to | Fixed — multiply-shift used; dead `kyber_compress` removed |
| uint16 nonce overflow in signing | pq-crystals/dilithium#110 | Counter uses Int; `% UInt16` at the SHAKE call site wraps at iteration 9362 for ML-DSA-87, as in the C reference (probability ~2^-23400) |
| Hint index monotonicity | GHSA-5x2r-hc65-25f9 | Correct (strict <) |
| sign_internal_msg domain separator | FIPS 204 §5.2 | Fixed |
| Context string length | FIPS 204 §5.2 | Validated (≤255 bytes) |
| Shamir prime too small for 256-bit keys | Functional testing | Fixed — p = 2^521-1 |
| ML-DSA keygen seed and hedged `rnd` drawn from Julia's default (non-cryptographic) RNG | Audit, v0.2.0 | Fixed — `RandomDevice` |
| ML-DSA signer emitted an unverifiable signature (~1 in 20,000 for ML-DSA-44): MakeHint was given w instead of HighBits(w) | Audit, v0.2.0 (20,000-signature sweep) | Fixed — single FIPS 204 Sign_internal core; regression test |
| `dilithium_sign_internal_msg` prepended the pure-mode prefix, so it was not Sign_internal | ACVP internal-interface vectors | Fixed — hashes M′ as given (FIPS 204 Alg. 7) |
| ML-KEM missing FIPS 203 input checks (ek modulus check, dk hash check, lengths) | pq-crystals/kyber#107, ACVP keyCheck vectors | Fixed — enforced in Encaps/Decaps; `kyber_ek_check`/`kyber_dk_check` exported |
| ML-DSA Verify threw on context > 255 bytes and on malformed keys | FIPS 204 Alg. 3 | Fixed — returns false |
| ML-DSA signing accepted secret keys with s1/s2 outside [−η, η] | Wycheproof `InvalidPrivateKey` | Fixed — rejected |
| Falcon FFT lost ~8 bits of the 53-bit mantissa (naive O(n²) power accumulation) | Comparison against a 256-bit reference ffLDL tree | Fixed — split/merge FFT with exactly-rounded roots |
| Falcon signed with any key whose ffLDL leaves passed the GS-norm test only on paper | NIST FIPS 206 status update (Oct 2025) | Fixed — signing refuses keys with a leaf outside [σmin, σmax]; `falcon_keygen(certified=true)` and `falcon_certify` decide the leaf and GS-norm bounds exactly |

## Falcon Key Certificate

`falcon_certify(sk)` recomputes every ffLDL leaf d as an exact rational: the pointwise LDL runs
in ℤ/p for NTT primes p ≡ 1 (mod 2n), then CRT and rational reconstruction under a Hadamard
bound give d (multi-modular method of Pornin, ePrint 2023/290). Each leaf is tested as
σ²/σmax² ≤ d ≤ σ²/σmin² with σ, σmin and σmax at their exact binary values, and the GS norm
against 1.17²·q in ℚ. Every call checks ∏d = qⁿ. The tests compare the leaves against a ℚ(x)
field tower for n ≤ 16.

A certificate takes about 5 s at n = 512 and 25 s at n = 1024; `certified=true` keygen
certifies each candidate key.

## FIPS 206

`fips206=true` applies the NIST FIPS 206 status update (Perlner, 2025):

- the keygen GS bound 0.9999·1.17√q, decided exactly by the certificate;
- the signer's refusal of keys with a leaf outside [σmin, σmax];
- randomized signing;
- G recomputed from (f, g, F);
- the public key in NTT form, ĥ = NTT(g)/NTT(f) (`Falcon.pk_ntt`, `Falcon.verify_poly_ntt`).

## FALCON-MRM

FALCON-MRM follows ePrint 2026/420. Implementation choices:

- **Domain-separation tags** for H1 and H2 (`FalconMRM.TAG_H1`, `FalconMRM.TAG_H2`).
- **Header byte:** 0x70 + log₂n, the unused cc = 11 value of the round-3 header.

## X25519 and X-Wing

`X25519.x25519` is the RFC 7748 function: it masks the top bit of u, accepts non-canonical u,
and returns the all-zero value for small-order inputs (Wycheproof `ZeroSharedSecret`). For
standalone Diffie–Hellman, reject an all-zero result (RFC 7748 §6.1). X-Wing needs no such
check: its combiner hashes the X25519 ciphertext and public key (draft-11 §6).

## Reporting Vulnerabilities

Report security issues to the repository owner via GitHub private vulnerability reporting.
Do not open public issues for security vulnerabilities.

## References

- [FIPS 203 — ML-KEM](https://csrc.nist.gov/pubs/fips/203/final)
- [FIPS 204 — ML-DSA](https://csrc.nist.gov/pubs/fips/204/final)
- [SP 800-227 — Recommendations for KEMs](https://csrc.nist.gov/News/2025/nist-publishes-sp-800-227)
- [KyberSlash FAQ](https://kyberslash.cr.yp.to/faq.html)
- [RUSTSEC-2025-0144 — ML-DSA Decompose timing](https://rustsec.org/advisories/RUSTSEC-2025-0144.html)
- [CVE-2024-37880 — Compiler-introduced branch in poly_frommsg](https://nvd.nist.gov/vuln/detail/CVE-2024-37880)
