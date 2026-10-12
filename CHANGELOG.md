# Changelog

## Unreleased

### Fixed

Defects found by an internal sweep of every subsystem; each has a regression test
(`test/regressions_mldsa.jl`, `test/regressions_mlkem_xwing.jl`, `test/regressions_falcon.jl`,
`test/regressions_falcon_mrm.jl`, `test/wipe.jl`, `test/runtests.jl`). The wipe checks run in a
child process and cover the success paths, every rejection path and injected exceptions.

- Secrets left in memory. Every secret buffer is now allocated before a `try` and zeroed in its
  `finally`, so it is wiped on every exit, an exception included.
  - ML-KEM: the basemul scratch held Â∘ŝ_k, t̂_k∘r̂_k or ŝ_k∘NTT(u)_k after every keygen,
    encapsulation and decapsulation, and ŝ_k is recoverable from the decapsulation one and the
    public ciphertext; decapsulation never wiped the re-encryption of m′, which is a function of
    the secret m′ when the ciphertext is rejected.
  - ML-DSA: each rejected attempt left its challenge input μ ‖ w1Encode(w1), c̃, c, NTT(c) and the
    SampleInBall stream; a secret key failing the s1/s2 range check left the decoded K, s1, s2
    and t0; a signing error left the hedged `rnd`.
  - X-Wing: an encapsulation key failing the ML-KEM check left the X25519 scalar, m and the X25519
    shared secret; the key is now checked before anything is derived from the seed.
  - FALCON-MRM: signing never wiped the targets t0, t1 (with the public c they give f and F), the
    sampler seed, the ChaCha20 state and keystream, or rejected candidates.
  - Falcon keygen left its 56-byte ChaCha20 seed, which alone regenerates the key, the generator,
    every candidate f, g and the sampler and floating-point scratch; decoding a byte key left
    f mod q, G mod q, NTT(f) and f·G; key expansion that refused a key left the decoded f, g, F, G;
    `Falcon.sign(msg, sk::NamedTuple)` never wiped the FFT basis and ffLDL tree it built; the byte-key
    API wiped the expanded key only when signing succeeded; `falcon_certify` left its temporaries;
    a rejected signing attempt left its candidate and partial encoding.
- `wipe!`:
  - wrote zeros into an immutable String through `codeunits`, corrupting string literals and
    Dict keys; read-only storage (code units, ranges, unknown array types) is now an
    `ArgumentError`;
  - on strided views, transposes and BitArrays the barrier received a null pointer, so LLVM
    dropped the element-wise zeroing as dead stores; it now receives the storage pointer;
  - cleared `length * sizeof(T)` bytes of a view or `Memory`, too few for padded element types.
- Shamir:
  - `shamir_reconstruct` with k ≤ 0 returned the secret 0;
  - `shamir_share` accepted n ≥ p (the share at x = p is the secret) and a composite p;
  - x-coordinates whose difference overflows `Int` reconstructed a wrong secret;
  - `shamir_share_bytes` left a hex `String` copy of the secret (bytes now go through
    `mpz_import`/`mpz_export`).
- Inputs with other axes: `kyber_ek_check` accepted an invalid encapsulation key and rejected
  every valid one when given a 0-based vector; `kyber_dk_check`, X25519, X-Wing, FALCON-MRM and
  the Falcon byte-key API threw or misread such inputs. `xwing_encaps_derand` threw
  `MethodError` for a `Memory` or range seed.
- Falcon keygen's bare `catch` around the NTRU solver swallowed `InterruptException` (Ctrl-C) and
  every other error; only the solver's `NTRUSolveFailure` is caught now.
- `FalconMRM.encode(μ, k)` returned k − 1 elements for odd k; it now returns k (the last is 0),
  and `decode` is its exact inverse.
- The runtime-dispatch check (`test/timing/dispatch.jl`) failed on Julia 1.12.1 on every ML-KEM,
  ML-DSA, Falcon-signing and X-Wing path: `Base.vcat` of byte vectors builds its assertion message
  through dynamic dispatch. Byte strings are joined with `catbytes` (src/bytes.jl); every row reads 0.
- Found by an adversarial review of the fixes above:
  - `catbytes` and the multi-buffer `wipe!` were written with a plain `xs...`, so Julia compiled
    them for a widened Vararg signature and every call on the ML-KEM keygen, ML-DSA signing and
    X-Wing paths splatted through runtime dispatch, which JET does not report. Both specialize on
    the argument count now; `test/specialization.jl` walks the compiled call graph of every
    secret path and fails on any widened Vararg specialization.
  - ML-DSA keygen wiped its seed expansion, s1, s2, NTT(s1), t and t0 only on success; it now
    wipes them, and the partial secret key, on every exit.
  - A Falcon secret key given as a non-`Vector` whose reads fail part-way (a key store) left a
    partial copy, or the decoded f and g, unwiped; Shamir reconstruction left a plain-data copy
    of `Int` shares.
  - `FalconFFT.fft` and `FalconCertify.nttm` threw on a range; `poly_mul_modq`, `poly_div_modq`
    and `is_invertible` read outside a vector with offset axes; `verify_poly` threw for a public
    key of another degree instead of returning `false`; `NTRUSolveFailure()` did not exist.
- `wipe!` on a non-`Array` argument recursed until a stack overflow; Falcon secret types printed
  their key material in the REPL (now `<secret>`); `shamir_reconstruct_bytes` ignored `nbytes`.

### Changed

- ML-DSA input errors (context over 255 bytes, unknown pre-hash algorithm, wrong lengths) raise
  `ArgumentError` instead of `ErrorException`.
- Shamir raises `ArgumentError` for k < 1, an empty byte secret, a non-prime modulus and n ≥ p;
  `shamir_reconstruct_bytes` raises it when the secret does not fit in `nbytes`.
- `wipe!` raises `ArgumentError` for read-only storage and for arrays it cannot write through.
- `FalconMRM.decode(z, B)` returns `nothing` unless B = max_bits(length(z)); `max_bits` and
  `encode` reject negative k.
- Falcon keygen zeroes the buffer its `randombytes` callback returns, as signing already did; a
  callback must return a fresh buffer on each call.
- `FNDSA.FalconNTRUGen.NTRUSolveFailure` signals an (f, g) without an NTRU solution; the
  fixed-point solver raises it where it raised `InexactError`. Keys are unchanged.

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
