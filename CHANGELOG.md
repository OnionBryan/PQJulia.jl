# Changelog

## Unreleased

### Security

- Falcon key expansion computes each ffLDL leaf σ/√d twice, the second time on an operand passed
  through an empty asm statement, and refuses the key if the results differ: one glitch in this
  square root, followed by typically 1–2×10⁶ signatures, recovers the key (Kaihara et al., ePrint
  2026/2046). Key expansion is 2–2.5% slower.
- Falcon signing (`falcon_sign`, `falcon_mrm_sign`) re-checks, branch-free, that every leaf of the
  expanded key lies in [σmin, σmax] before each signature (2.6 µs at n = 512, 6.6 µs at n = 1024,
  about 0.06% of signing). A wiped `ExpandedKey` now raises `ArgumentError` instead of resampling
  forever. Signatures are unchanged bit for bit (`test/falcon_hardening.jl`).

### Tests

- `test/attack_regressions.jl`:
  - ML-KEM decapsulation of a one-bit flip in every ciphertext byte, at all three levels and
    through X-Wing, returns exactly SHAKE256(z ‖ c′, 32); the last c₂ coefficient is swept through
    all values (ePrint 2026/2239). A comparison that skips one middle byte passes the rest of the
    suite and fails here.
  - ML-DSA: every `invntt!` input during key generation, signing and verification, recorded by
    probe modules compiled from the shipped source, stays below q; with the key-generation or
    verification reduction removed it reaches 2–3q and the test fails (a missing key-generation
    reduction previously passed the whole suite). The reduce/inverse-NTT sequence is also checked
    against a schoolbook product on sign-aligned inputs that overflow Int32 without the reduction
    (ePrint 2026/1032).
  - The public `dilithium_sign` and `dilithium_sign_prehash` with `hedged=false` reproduce the 90
    deterministic ACVP SigGen vectors, which previously reached only the `_derand` and internal
    entry points (Bernstein, ML-DSA bug study, 2026).

### Documentation

- SECURITY.md separates timing from physical side channels: no code is masked, and the published
  profiled attacks on Falcon signing are cited (ePrint 2026/2124, 2025/2159, 2026/1366, 2026/2170),
  together with the fault and corruption checks above.
- SECURITY.md: X25519 offers no protection against a quantum adversary (2026 Shor-ECDLP resource
  estimates, arXiv 2603.28846, 2607.13816, 2609.05625); X-Wing's post-quantum security rests on
  ML-KEM-768. The skipped all-zero X25519 check is justified by the X-Wing proof (ePrint 2024/039).
- README: binary64 precision of the Falcon signer (TWFalcon, TCHES 2026); expand the key once for
  repeated signing, with the fault trade-off; the points ePrint 2026/420 leaves open for
  FALCON-MRM and PQJulia's choices.
- FIPS 206 is described as forthcoming (no NIST draft or final text as of 2026-10-09); X-Wing as an
  individual Internet-Draft in the Independent Submission stream.

### Performance

Medians on one x86-64 core, each scheme in a fresh process; outputs are unchanged bit for bit.

- Keccak-f[1600] in the package (`src/keccak.jl`): the 25-lane state is an immutable tuple and
  the round is straight-line code with constant rotations. SHAKE128 (672 bytes) takes 2.1 µs
  against 5.9 µs in SHA.jl, and the state leaves no heap copy. Every SHAKE and SHA3 call goes
  through it, including the SHA3 pre-hashes of HashML-DSA; SHA.jl remains for the SHA-2 ones.
- ML-KEM: `montgomery_reduce` drops the trailing `rem` by q, a no-op for every operand the
  scheme produces (inside the C reference's range −q·2¹⁵ ≤ a < q·2¹⁵, checked on all ACVP,
  Wycheproof and CCTV vectors); basemul
  works on scalars instead of SubArrays; compression no longer allocates per 8 coefficients;
  NTT, sampling and packing loops check bounds once. 3.0–3.5× at every level.
- ML-DSA: pointwise multiply-accumulate in one pass, bounds checked once in the NTT, sampling
  and packing loops. 2.0–3.2×.
- Falcon verification: per-degree NTT tables replace the per-call root search and `powermod`
  per coefficient; butterflies run in UInt32 with branch-free reduction; signature
  (de)compression reads and writes bits in place. 5.6× (Falcon-512) and 5.7× (Falcon-1024).
- Falcon signing: ChaCha20 runs on registers and refills its buffer in place, SamplerZ reads its
  bytes from that buffer, and the FFT roots come from a table built at load. 35% fewer
  allocations; 3–11% faster.

## 0.3.1

- Falcon secret-key decoding is constant-time: Fermat inversion with a fixed exponent replaces the
  extended Euclidean algorithm, and the NTRU check runs in Int64 instead of BigInt. dudect had
  flagged both on an AMD EPYC 9B14.

## 0.3.0

### Security

- Falcon signs without hardware floating point. Key expansion, the FFT, ffLDL, ffSampling and
  SamplerZ run on an integer emulation of binary64 (Pornin's `fpr`), branch-free; signatures are
  unchanged bit for bit.
- X25519 decodes the scalar and encodes the shared secret in 64-bit words, without BigInt, whose
  timing depended on the value.
- ML-DSA: `make_hint` is branch-free on x86-64; key generation no longer boxes the η-sampling
  counter.
- Secret buffers are zeroed after use. `wipe!` is exported for buffers the caller owns, and
  `falcon_wipe!` clears an expanded Falcon key.

### Testing

- `test/timing/dudect.jl`: dudect timing test for ML-KEM, ML-DSA, X25519 and Falcon.
- `test/timing/dispatch.jl`: JET runtime-dispatch check on every secret-key path, run in CI.
- CI on Linux x86-64 and ARM64, macOS and Windows, plus two-way interop with liboqs 0.16.0.

## 0.2.0

- X25519 (RFC 7748) and the X-Wing hybrid KEM (draft-connolly-cfrg-xwing-kem-11).
- Falcon: exact key certificate, float-free key generation, FALCON-MRM, FIPS 206 options.
