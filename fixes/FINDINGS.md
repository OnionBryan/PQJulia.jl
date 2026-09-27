# PQJulia.jl — Security Research Findings

Compiled 2026-04-02 from: NIST FIPS 203/204 specs, pq-crystals GitHub issues,
RUSTSEC/GHSA advisories, KyberSlash disclosure, C2SP CCTV test vectors,
IETF security considerations drafts, and ePrint papers (2023-2026).

Baseline: 677/677 KATs passing.

## Critical Issues (must fix before public)

### C1. `poly_uniform!` missing re-squeeze loop (Dilithium)
- **File:** `dilithium_core.jl`
- **Bug:** Fixed 1280-byte SHAKE128 buffer, `error()` on exhaustion
- **C ref:** `poly.c:359-367` — `while(ctr < N)` re-squeeze loop, one block at a time
- **Fix ready:** `fixes/poly_uniform_fix.jl`

### C2. `poly_uniform_eta!` missing re-squeeze loop (Dilithium)
- **File:** `dilithium_level.jl`
- **Bug:** Fixed 512-byte SHAKE256 buffer, `error()` on exhaustion
- **C ref:** `poly.c:449-452` — `while(ctr < N)` re-squeeze loop
- **Fix ready:** `fixes/poly_uniform_eta_fix.jl`

### C3. `polyt0_pack` inlined twice (Dilithium)
- **File:** `dilithium_level.jl`
- **Bug:** Copy-pasted in keygen and sign, divergence risk
- **Fix ready:** `fixes/polyt0_pack_fix.jl`

### C4. `uint16_t nonce` overflow in signing loop (Dilithium)
- **File:** `dilithium_level.jl`
- **Bug:** `nonce::UInt16` wraps at 65535. For ML-DSA-87 (L=7), `nonce += L` wraps
  from 65534→5 after 9362 iterations, causing nonce reuse in mask generation.
- **Source:** pq-crystals/dilithium#110 (2026-03-31), confirmed in Microsoft SymCrypt
- **Probability:** ~2^-23400 for random keys, but exists and is a spec violation
- **Fix:** Change `nonce::UInt16` to `nonce::UInt32` (or `nonce::Int`)

### C5. KyberSlash check — division by q in compress/decompress (Kyber)
- **File:** `kyber_core.jl`
- **Bug:** If `kyber_compress` or `kyber_poly_tomsg!` uses integer division `÷ q`
  on secret-derived values, timing leaks the secret key
- **Source:** kyberslash.cr.yp.to, TCHES 2024
- **Status:** Agent analysis shows PQJulia uses multiply-shift (`80635 >> 28`).
  **VERIFY** — if already fixed, document as "mitigated"

## High Priority (should fix before public)

### H1. Hint index monotonicity — strict `<` not `<=`
- **File:** `dilithium_level.jl` (verify path)
- **Bug:** If hint indices use `<=` instead of `<`, duplicate indices pass,
  breaking strong unforgeability
- **Source:** GHSA-5x2r-hc65-25f9 (RustCrypto, Jan 2026)
- **Fix:** One-character change, but verify against Wycheproof ML-DSA vectors

### H2. Context string length validation
- **File:** `dilithium_level.jl`
- **Bug:** FIPS 204 limits context string to ≤255 bytes. No bounds check exists.
- **Fix:** `length(context) > 255 && error("Context string must be ≤ 255 bytes")`

### H3. Shamir x=0 guard
- **File:** `shamir.jl`
- **Bug:** Share at x=0 IS the secret (constant term of the polynomial)
- **Fix:** `any(x -> x == 0, 1:n) && error("Share index 0 leaks the secret")`
  (shares are 1-indexed so this can't happen currently, but the API allows it)

## Moderate Priority (documentation / hardening)

### M1. Document constant-time boundary
- Julia JIT makes CT guarantees impossible (CVE-2024-37880 class)
- `kyber_verify`/`kyber_cmov!` are CT at source level but LLVM may optimize away
- Document: "This is research software. CT not guaranteed. For production, use ccall to liboqs."

### M2. `Decompose` uses integer division on secret-derived values
- **File:** `dilithium_level.jl`
- **Bug:** `decompose` function uses division that may be data-dependent on ARM/x86
- **Source:** RUSTSEC-2025-0144 (RustCrypto)
- **Fix for production:** Barrett reduction. For research software: document the risk.

### M3. Misleading masks in `polyt0_unpack` / `polyz_unpack`
- **Source:** pq-crystals/dilithium#55
- **Status:** Cosmetic only, no correctness impact. Low priority.

## Not Fixing (by design)

### N1. Signing iteration cap
- FIPS 204 specifies NO cap. Adding one leaks information.
- The loop is unbounded by design (Fiat-Shamir with Unbounded Aborts).
- **Decision:** Do not add a cap. This reverses my earlier recommendation.

### N2. GC/JIT timing variability
- Julia-specific. Cannot be fixed without leaving Julia.
- **Decision:** Document in security model, don't pretend to fix.

### N3. Rejected signature side-channel leakage
- ePrint 2025/582, 2025/214, 2026/056 — but hardware/embedded only.
- **Decision:** Out of scope for pure software.

## Sources
- FIPS 203/204 final standards (NIST CSRC)
- C2SP/CCTV ML-KEM test vectors (XOF buffer edge cases)
- C2SP/Wycheproof ML-DSA vectors (hint monotonicity)
- pq-crystals/kyber#107 (missing FIPS 203 input validations)
- pq-crystals/dilithium#110 (uint16 nonce overflow)
- KyberSlash FAQ (kyberslash.cr.yp.to)
- RUSTSEC-2025-0144 (Decompose timing)
- GHSA-5x2r-hc65-25f9 (hint index <=)
- IETF draft-connolly-cfrg-ml-dsa-security-considerations-01
- IETF draft-sfluhrer-cfrg-ml-kem-security-considerations-04
- SP 800-227 ipd (NIST implementation guidance)
- ePrint 2025/450, 2025/582, 2025/214, 2026/056, 2026/472