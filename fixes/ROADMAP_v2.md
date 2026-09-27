# PQJulia.jl — Roadmap v2 (Post-Audit)

All critical issues from v1 are fixed and pushed. 677/677 KATs pass.

## Current State

| Component | Lines | Status |
|-----------|-------|--------|
| ML-KEM (FIPS 203) | 770 | Complete, all 3 levels, KAT-validated |
| ML-DSA (FIPS 204) | 1,213 | Complete, all 3 levels, pure + prehash, KAT-validated |
| Shamir SSS | 38 | Functional, minimal |
| Total | ~2,657 | 677 ACVP KATs passing |

## What Was Fixed (pushed)

1. `poly_uniform!` — SHAKE128 re-squeeze loop (C ref: poly.c:359-367)
2. `poly_uniform_eta!` — SHAKE256 re-squeeze loop (C ref: poly.c:449-452)
3. `poly_challenge!` — SHAKE256 re-squeeze for SampleInBall
4. `kyber_sample_uniform!` — extend XOF output instead of re-hashing same input
5. `polyt0_pack/unpack!` — extracted from 3 inline copies to shared functions
6. Signing nonce `UInt16` → `Int` (pq-crystals/dilithium#110)
7. Context string length validation (FIPS 204 §5.2, ≤255 bytes)
8. Removed dead `kyber_compress` function (KyberSlash-vulnerable if ever called)

## What's Next — Paper-Backed Directions

### 1. Hybrid KEM Demo (X25519 + ML-KEM)

**Why now:** Two 2025-2026 papers validate the approach:
- Chen et al. (Entropy 2025): CPA-secure KEMs suffice for TLS 1.3 hybrid — 44.8% speedup
  by skipping FO re-encryption. CPA ML-KEM is simpler to implement correctly.
- Sharma et al. (COMSNETS 2026): User-space ML-KEM-768 as WireGuard PSK injection,
  108ms mean latency, minimal CPU overhead.

**Concrete task:** Build `hybrid_kem(x25519_pk, mlkem_pk)` that combines shared secrets via
HKDF. Needs a Curve25519 dependency (or FFI to libsodium).

**Effort:** 3 (needs external dep) | **Impact:** 5

### 2. SLH-DSA (FIPS 205) — If Pursuing

**State of the art:** ePrint 2025/2236 extends SPHINCS+ framework with variable tree heights
and chain lengths, finding 8-26% smaller signatures than standard SPHINCS+-128s.

**What it requires:**
- WOTS+ one-time signatures (hash chains)
- FORS few-time signatures (forest of random subsets)
- Hypertree: Merkle trees at multiple layers
- Address compression scheme (ADRS)
- All SHAKE256, no NTT — completely different from Kyber/Dilithium

**Honest assessment:** The math is simpler (trees + hashes, no lattices) but the implementation
is fiddly (address bytes must be exact, hypertree traversal order matters). If you got Dilithium's
rejection sampling right, you can do this. Estimated 1,500-2,000 lines.

**Effort:** 5 | **Impact:** 5 (completes NIST trifecta)

### 3. Benchmark Suite

**Why:** Zheng et al. (arXiv 2024, 10 citations) shows AVX-512 ML-KEM achieves 1.64x speedup
and batch keygen 3.5-4.9x speedup. Julia's SIMD support could achieve similar gains if the
NTT is the bottleneck — but we need benchmarks to know.

**Concrete task:**
```julia
using BenchmarkTools
@benchmark MLKEM.Category3.kyber_kem_keypair()
@benchmark MLKEM.Category3.kyber_kem_enc(pk)
@benchmark MLDSA.Category3.dilithium_sign(msg, sk)
# 33 benchmarks: 6 KEM ops + 5 DSA ops × 3 levels
```

**Effort:** 2 | **Impact:** 4

### 4. Side-Channel Documentation

**Paper evidence (do NOT try to fix in Julia):**
- ePrint 2026/056: Non-profiling attack on ML-DSA via public templates, 96 traces for
  challenge recovery, 300 for key recovery on ARM Cortex-M4
- ePrint 2025/582: Rejected signatures reduce traces by 50%+; rejected-only key recovery
  in <30 traces
- ePrint 2025/1629: Cauchy regression breaks MASKED Dilithium at all NIST levels in <2 min
- ePrint 2024/843: Formally verified ML-KEM in EasyCrypt with provably CT Jasmin
  implementation — this is the gold standard PQJulia cannot match

**What to document:**
- Julia's JIT makes CT guarantees impossible (CVE-2024-37880 class)
- `kyber_verify`/`kyber_cmov!` are CT at source level but LLVM may optimize away
- Hedged mode does NOT protect against SCA (ePrint 2026/056 confirms)
- For production: use `ccall` to libOQS, HACL*, or mlkem-native

**Effort:** 2 | **Impact:** 4

### 5. Shamir Hardening

**Current gaps:**
- No duplicate x-index guard (throws opaque `ArgumentError` from `mod_inverse`)
- No byte encoding layer (raw BigInt in/out)
- Uses `powermod` evaluation instead of Horner (correct but slower)
- No Feldman VSS commitment scheme

**Feldman VSS literature:** No recent IACR results for PQ Feldman specifically. Classical
Feldman works over groups where DLP is hard — needs rethinking for PQ setting. Pedersen
commitments are an alternative but same issue. Skip for now unless there's a use case.

**Effort:** 2 | **Impact:** 3

### 6. Property-Based Tests

**What to test:**
- Pack/unpack round-trip for all polynomial serialization formats
- NTT/INTT round-trip: `invntt!(ntt!(copy(a))) ≈ a` for random polynomials
- Compress/decompress approximate round-trip within tolerance
- Sign/verify: `verify(msg, sign(msg, sk), pk) == true` for random msg

**Effort:** 2 | **Impact:** 4

### 7. CI Matrix

**Task:** GitHub Actions with Julia 1.12 stable. Run KAT tests on push.
Add warning if `test/kat/` directory is missing (currently silently skips all 677 KATs).

**Effort:** 1 | **Impact:** 4

## Revised Priority (Impact²/Effort)

| Direction | Impact | Effort | Score |
|-----------|--------|--------|-------|
| CI + kat-missing warning | 4 | 1 | 16.0 |
| Benchmark suite | 4 | 2 | 8.0 |
| Property-based tests | 4 | 2 | 8.0 |
| CT documentation | 4 | 2 | 8.0 |
| Hybrid KEM demo | 5 | 3 | 8.3 |
| Shamir hardening | 3 | 2 | 4.5 |
| SLH-DSA (FIPS 205) | 5 | 5 | 5.0 |

## References

- Chen et al. "On the Security and Efficiency of TLS 1.3 Handshake with Hybrid Key Exchange
  from CPA-Secure KEMs" (Entropy 2025) — doi:10.3390/e27121242
- Sharma et al. "User-Space PQ Key Exchange in WireGuard Using ML-KEM-768"
  (COMSNETS 2026) — doi:10.1109/COMSNETS67989.2026.11418293
- Tu & Xie. "EMINEM: Efficient FPGA Implementation of Mixed-RadIx NTT Hardware
  AccElerators for NIST PQC" (ACM 2025) — doi:10.1145/3771287
- Mandal & Roy. "Winograd for NTT" (TCSI 2024) — doi:10.1109/TCSI.2024.3470335
- ePrint 2025/2236: "Extending the SPHINCS+ Framework" (SLH-DSA improvements)
- ePrint 2026/056: "Rejection Matters" (non-profiling ML-DSA SCA, 96 traces)
- ePrint 2025/582: "Release the Power of Rejected Signatures" (50% trace reduction)
- ePrint 2025/1629: "Solving Concealed ILWE" (breaks masked Dilithium via Cauchy regression)
- ePrint 2024/843: "Formally verifying ML-KEM in EasyCrypt" (gold standard CT proof)
- pq-crystals/dilithium#110: uint16_t nonce overflow (2026-03-31)
- pq-crystals/kyber#107: missing FIPS 203 input validations
- RUSTSEC-2025-0144: ML-DSA Decompose timing side-channel
- GHSA-5x2r-hc65-25f9: ML-DSA hint index monotonicity
- C2SP/CCTV: ML-KEM "unlucky" XOF test vectors
- SP 800-227 ipd: NIST PQC implementation guidance