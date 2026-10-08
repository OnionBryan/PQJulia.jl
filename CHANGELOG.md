# Changelog

## Unreleased

### Performance

Medians on one x86-64 core, each scheme in a fresh process; outputs are unchanged bit for bit.

- Keccak-f[1600] in the package (`src/keccak.jl`): the 25-lane state is an immutable tuple and
  the round is straight-line code with constant rotations. SHAKE128 (672 bytes) takes 2.1 µs
  against 5.9 µs in SHA.jl, and the state leaves no heap copy. SHA.jl remains for the SHA-2
  pre-hashes of HashML-DSA.
- ML-KEM: `montgomery_reduce` drops the trailing `rem` by q, a no-op for every operand the
  scheme produces (|a| ≤ q·2¹⁵, checked on all ACVP, Wycheproof and CCTV vectors); basemul
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
