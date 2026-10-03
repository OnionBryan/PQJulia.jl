# Changelog

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
