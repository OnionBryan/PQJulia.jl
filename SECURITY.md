# Security Policy

## Status

PQJulia.jl is research software. It has not received independent security audit.
Do not use it to protect sensitive data in production without additional review.

## What This Library Is

A pure-Julia implementation of NIST FIPS 203 (ML-KEM) and FIPS 204 (ML-DSA),
validated against 677 NIST ACVP Known Answer Tests and cross-validated against
the pq-crystals C reference implementation.

## What This Library Is NOT

- **Not constant-time.** Julia compiles via LLVM JIT. The compiler may introduce
  data-dependent branches, speculative execution, or memory access patterns that
  leak secret information through timing side channels. Source-level constant-time
  idioms (e.g., `kyber_verify`, `kyber_cmov!`) are best-effort but not guaranteed
  to survive compilation. See Trail of Bits' `__builtin_ct_select` RFC (Dec 2025)
  for the current state of LLVM CT support.

- **Not hardened against physical side channels.** Power analysis, electromagnetic
  emanation, and fault injection attacks are out of scope. Recent work shows
  ML-DSA key recovery from rejected signatures in <100 power traces on ARM
  Cortex-M4 (ePrint 2026/056, 2025/582). Hedged (non-deterministic) mode does
  not protect against these attacks (ePrint 2026/056).

- **Not FIPS 140-3 validated.** FIPS 140-3 module validation requires side-channel
  resistance at the module level. This library does not meet that requirement.

## Threat Model

PQJulia.jl is appropriate for:
- Research and prototyping
- Algorithm correctness verification
- Educational purposes
- Applications where timing side channels are not in the threat model
  (e.g., offline batch processing, local-only tools)

PQJulia.jl is NOT appropriate for:
- TLS servers or other network-facing cryptographic services
- Hardware security modules or embedded devices
- Any context where an attacker can measure operation timing
- FIPS 140-3 validated environments

## Constant-Time Analysis

### What IS constant-time at the source level

| Function | Mechanism | File |
|----------|-----------|------|
| `kyber_verify` | OR-accumulate pairwise-XOR differences, `(-UInt64(r)) >>> 63` | kyber_core.jl |
| `kyber_cmov!` | Byte-wise mask with `UInt8(0) - b` | kyber_core.jl |
| `kyber_poly_tomsg!` | Barrett multiply-shift, no division | kyber_core.jl |
| `kyber_poly_compress!` | Barrett multiply-shift (80635>>28, 40318>>27) | kyber_core.jl |
| ML-KEM Decaps rejection | Always computes both paths, `cmov` selects | kyber_kem.jl |

### What is NOT constant-time

| Operation | Risk | Mitigation |
|-----------|------|------------|
| Julia GC pauses | Hundreds of μs during crypto ops | None in pure Julia |
| JIT compilation | First call slower than subsequent | Pre-warm before use |
| LLVM optimizations | May re-introduce branches (CVE-2024-37880) | Cannot control in Julia |
| `decompose` (ML-DSA) | Uses multiply-shift on secret data (no UDIV); branches on GAMMA2 are public constants folded at JIT time. RUSTSEC-2025-0144 does not apply to this implementation. | No action needed |
| SHA.jl SHAKE implementation | Pure Julia, timing depends on input length | Acceptable (length is public) |

### For production use

Use `ccall` into a validated C implementation:
- [liboqs](https://github.com/open-quantum-safe/liboqs) (Open Quantum Safe)
- [mlkem-native](https://github.com/pq-code-package/mlkem-native) (formally verified)
- [HACL*](https://hacl-star.github.io/) (verified in F*)

## Known Issues Addressed

| Issue | Source | Status |
|-------|--------|--------|
| XOF buffer exhaustion (poly_uniform, poly_uniform_eta, poly_challenge, kyber_sample_uniform) | C ref comparison + audit | Fixed — re-squeeze loops added |
| KyberSlash (division by q) | kyberslash.cr.yp.to | Mitigated — multiply-shift used; dead `kyber_compress` removed |
| uint16 nonce overflow in signing | pq-crystals/dilithium#110 | Partially mitigated — counter uses Int but `% UInt16` at SHAKE call site reintroduces wrap at iteration 9362 for ML-DSA-87. Same limitation as the C reference. Probability of reaching 9362 iterations: ~2^-23400 |
| Hint index monotonicity | GHSA-5x2r-hc65-25f9 | Already correct (strict <) |
| sign_internal_msg domain separator | FIPS 204 §5.2 | Fixed |
| Context string length | FIPS 204 §5.2 | Validated (≤255 bytes) |
| Shamir prime too small for 256-bit keys | Functional testing | Fixed — p = 2^521-1 |

## Reporting Vulnerabilities

Report security issues to the repository owner via GitHub private vulnerability reporting.
Do not open public issues for security vulnerabilities.

## References

- [FIPS 203 — ML-KEM](https://csrc.nist.gov/pubs/fips/203/final)
- [FIPS 204 — ML-DSA](https://csrc.nist.gov/pubs/fips/204/final)
- [SP 800-227 — Recommendations for KEMs](https://csrc.nist.gov/News/2025/nist-publishes-sp-800-227)
- [KyberSlash FAQ](https://kyberslash.cr.yp.to/faq.html)
- [Trail of Bits — LLVM CT support](https://blog.trailofbits.com/2025/12/02/introducing-constant-time-support-for-llvm-to-protect-cryptographic-code/)
- [ePrint 2026/056 — Non-profiling SCA on ML-DSA](https://eprint.iacr.org/2026/056)
- [ePrint 2025/582 — Rejected signatures SCA](https://eprint.iacr.org/2025/582)
- [ePrint 2024/843 — Formally verified ML-KEM](https://eprint.iacr.org/2024/843)
- [RUSTSEC-2025-0144 — ML-DSA Decompose timing](https://rustsec.org/advisories/RUSTSEC-2025-0144.html)
- [CVE-2024-37880 — Compiler-introduced branch in poly_frommsg](https://nvd.nist.gov/vuln/detail/CVE-2024-37880)
- [liboqs Security Policy](https://github.com/open-quantum-safe/liboqs/blob/main/SECURITY.md)
