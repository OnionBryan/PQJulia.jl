# Security Policy

## Scope

A pure-Julia implementation of NIST FIPS 203 (ML-KEM), FIPS 204 (ML-DSA), Falcon (round-3
specification v1.2, the basis of FN-DSA, to be standardized as FIPS 206), X25519 (RFC 7748) and
X-Wing (an Internet-Draft). ML-KEM and ML-DSA pass 855 NIST ACVP vectors covering every
interface. Falcon signing is byte-exact with the round-3 C implementation on 120 reference
vectors.

## Randomness

All secret randomness comes from the operating system CSPRNG (`Random.RandomDevice`):
ML-KEM keygen and encapsulation, ML-DSA keygen and hedged `rnd`, Falcon keygen, salts and
per-signature sampler seeds (expanded with ChaCha20, as in the reference), X25519 and X-Wing
keys, and Shamir coefficients. Julia's default task-local RNG (Xoshiro) is never used for key
material. ML-DSA signing is hedged by default, as FIPS 204 recommends; `hedged=false` selects
the deterministic variant.

## Timing

This section concerns execution time only; power, electromagnetic and fault attacks are covered
under [Physical side channels](#physical-side-channels). Julia compiles through LLVM, which may
introduce data-dependent branches; source-level constant-time code is not guaranteed to stay
constant-time after compilation.

### Constant-time at the source level

| Function | Mechanism | File |
|----------|-----------|------|
| `kyber_verify` | OR-accumulate pairwise-XOR differences, `(-UInt64(r)) >>> 63` | kyber_core.jl |
| `kyber_cmov!` | Byte-wise mask with `UInt8(0) - b` | kyber_core.jl |
| `kyber_poly_tomsg!` | Barrett multiply-shift, no division | kyber_core.jl |
| `kyber_poly_compress!` | Barrett multiply-shift (80635>>28, 40318>>27) | kyber_core.jl |
| ML-KEM Decaps rejection | Always computes both paths, `cmov` selects | kyber_kem.jl |
| `decompose` (ML-DSA) | Multiply-shift, no division (RUSTSEC-2025-0144 does not apply); GAMMA2 branches are on public constants | dilithium_level.jl |
| `make_hint` (ML-DSA) | Bitwise `\|`/`&` comparisons; compiles to `setcc`/`csel`, no branches, on x86-64 and AArch64 | dilithium_level.jl |
| X25519 | Scalar bits from 64-bit words, mask swaps on 51-bit limbs, branch-free canonical encoding | x25519.jl |
| Falcon signing (key expansion, FFT, ffLDL, ffSampling, SamplerZ) | Integer emulation of binary64 (Pornin's `fpr`): branch-free add, mul, div, sqrt and rounding; no hardware floating point; unmasked | falcon/falcon_fpr.jl |

### Variable-time

| Operation | Note |
|-----------|------|
| Julia GC pauses, JIT compilation | First call compiles; warm up before timing-sensitive use |
| `poly_uniform_eta!` (ML-DSA keygen) | Rejection sampling, as in the reference: timing shows which bytes were rejected, and those are independent of the kept coefficients |
| Matrix expansion (ML-KEM, ML-DSA) | Rejection sampling on the public seed ρ |
| SHAKE, SHA3 (keccak.jl) | Timing depends on input length, which is public |
| SamplerZ rejection and BerExp, signature retries (Falcon) | As in the reference: iteration counts and the bytes compared depend on fresh randomness |
| Falcon keygen (NTRU solver) | BigInt arithmetic; `fixedpoint=true` removes floating point (ePrint 2023/290) |

## Physical side channels

No code is masked, and resistance to an attacker who can measure or perturb the device (power,
electromagnetic or fault attacks) is out of scope. ML-KEM, ML-DSA, X25519 and X-Wing claim
nothing beyond the timing behavior above.

Falcon signing is the most exposed. `fpr_of`, `floor` and `mul` in `falcon_fpr.jl`, and
ffSampling and SamplerZ above them, are branch-free but unmasked: they compute the signs,
exponent classes and normalization masks of secret values as data, and they port the reference
C line for line. Published profiled attacks on these operations report key recovery from 20–100
signatures through the conversion of ⌊μ⌋ at the ffSampling leaves (ePrint 2026/2124), from one
trace of key expansion (ePrint 2025/2159), and from about 100 signatures through the signs of
SamplerZ outputs (ePrint 2026/1366); the first and last derive their Falcon figures from
simulated leakage labels. ePrint 2026/2170 adds a leak through the sign of SamplerZ's input
center. Signing from an `ExpandedKey` (`falcon_expand_sk`) runs key expansion once per key
instead of once per signature, which limits but does not remove the exposure to 2025/2159; a
fault injected during expansion, in contrast, persists in every signature made from that key.

Two checks target faults and corruption of the expanded key. Key expansion computes each ffLDL
leaf σ/√d twice, the second time on an operand passed through an empty asm statement so that the
compiler keeps both evaluations, and throws if the two differ: one glitch in this square root
followed by typically 1–2×10⁶ signatures recovers the key (Kaihara et al., ePrint 2026/2046), and
the leaf range check alone still leaves 34–65% of keys recoverable. Before every signature the
signer re-checks, branch-free, that each leaf lies in [σmin, σmax], so a corrupted or wiped
`ExpandedKey` is refused. Neither check detects a fault that alters both evaluations identically,
a fault earlier in the FFT or LDL, or a change that leaves a leaf inside [σmin, σmax].

## Secrets in memory

Every secret buffer the package allocates is zeroed before the call returns, on every exit: a
result, a rejection, or an exception. Secret buffers are allocated before a `try` and wiped in its
`finally`. `wipe!` writes zeros and then passes the storage pointer to an empty inline-assembly
statement, so the compiler cannot drop the zeroing as a dead store (`test/wipe.jl` checks the
compiled code, including strided views). It accepts arrays whose storage is an `Array`, `Memory`
or `BitArray` (views, reshapes, reinterpretations, transposes of them) and raises `ArgumentError`
for read-only storage such as a String's code units or a range. Secret keys are built into storage
reserved at full size, so no partial copy is left behind by a reallocation. The regression tests
(`test/regressions_*.jl`) check each path in a child process by recording the buffers the code
allocates and, for Falcon on Linux, by scanning the process memory for key material.

| Scheme | Wiped |
|--------|-------|
| ML-KEM | keygen seed expansion, σ, s and e, Â∘ŝ, PRF and noise buffers; encapsulation (K, r), r, e1, e2, Âᵀ∘r̂ and t̂ᵀ∘r̂; decapsulation ŝ, ŝᵀ∘û, m′, (K′, r′), the re-encryption of m′ and the rejection-key input; the basemul scratch on every path. Decapsulation reads dk_PKE and z in place, without a copy |
| ML-DSA | keygen ξ, seed expansion, ρ′, K, s1, s2, t and t0; when signing, the decoded K, s1, s2 and t0, ρ″, y, z, w0, w1 and the hints, and for every attempt, rejected or accepted, the challenge input μ ‖ w1Encode(w1), c̃, c, NTT(c) and the SampleInBall stream; the hedged `rnd` |
| X-Wing | the expanded ML-KEM key and X25519 scalar, the encapsulation seed and both partial shared secrets; the ML-KEM key check runs before anything is derived from the seed |
| X25519 | holds its scalar and field elements in immutable tuples, which leave no heap copy |
| SHAKE, SHA3 | the Keccak state is an immutable tuple (keccak.jl), which leaves no heap copy |
| Falcon | whenever the API decodes a byte key for one call: the decoded f, g, F, G, every mod-q, NTT and product copy that decoding computes, the FFT basis and the ffLDL tree, also when expansion refuses the key or signing throws; in keygen, the ChaCha20 seed (which alone regenerates the key) and generator, every candidate f, g, accepted or rejected, and the sampler, floating-point and NTT scratch; when signing (Falcon and FALCON-MRM), the targets t0, t1, the FFT and ffSampling arrays, the ChaCha20 seed, state and keystream, which SamplerZ reads in place, and every rejected candidate and its partial encoding |

What the caller owns stays the caller's to erase: secret keys, returned shared secrets, and an
expanded Falcon key (`falcon_wipe!(ek)`). `wipe!` is exported for this. Falcon keygen and signing
zero the buffers their `randombytes` callback returns, so a callback must return a fresh buffer on
each call. Not wiped: BigInt (GMP) values in Falcon keygen, in `falcon_certify` and in Shamir
(the byte API converts with `mpz_import`/`mpz_export` and leaves no String or byte copy), and a
helper's own temporary when an asynchronous exception (`InterruptException`, `OutOfMemoryError`)
interrupts it. Julia's garbage collector does not move objects, but it does not lock pages
either, so memory can reach swap or a core dump.

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
| X25519 decoded the scalar and encoded the shared secret through BigInt (GMP), whose timing depends on the value (dudect \|t\| ≈ 20 on Apple M5 and Intel Broadwell) | `test/timing/dudect.jl` | Fixed — 64-bit word arithmetic, no BigInt on secret data; type-stable, so no value-dependent boxing |
| ML-DSA keygen boxed the rejection counter of `poly_uniform_eta!` (a nested function reassigned an enclosing variable), so it ran through runtime dispatch on secret-derived data | `test/timing/dispatch.jl` (JET) | Fixed — top-level `rej_eta!`; type-stable |
| Falcon signed with hardware floating point on secret data, which is not constant-time | Round-3 reference (FALCON_FPEMU) | Fixed — integer-emulated binary64; signatures unchanged bit for bit |
| Falcon secret-key decoding inverted NTT coefficients of f with the extended Euclidean algorithm and checked the NTRU equation in BigInt, both variable-time on secret data (dudect \|t\| up to 27.6, AMD EPYC 9B14) | `test/timing/dudect.jl` | Fixed — Fermat inversion with a fixed exponent; NTRU check in Int64 |
| `make_hint` compiled to conditional branches on x86-64 | Disassembly after a dudect flag (\|t\| 7.4, Intel i7-8086K) | Fixed — bitwise form, branch-free |
| Falcon signed with any key whose ffLDL leaves passed the GS-norm test only on paper | NIST FIPS 206 status update (Oct 2025) | Fixed — signing refuses keys with a leaf outside [σmin, σmax]; `falcon_keygen(certified=true)` and `falcon_certify` decide the leaf and GS-norm bounds exactly |
| One fault in the ffLDL leaf square root at key expansion biases every later signature and recovers the key; the leaf range check alone leaves 34–65% of keys recoverable | Kaihara et al., ePrint 2026/2046 (CCS 2026) | Mitigated — each leaf computed twice, expansion refuses the key on a mismatch; leaf range re-checked, branch-free, before every signature ([Physical side channels](#physical-side-channels)) |
| A decapsulation check that leaves one ciphertext coordinate unverified permits key recovery (no-op `cmov` in pqc_kyber, truncated comparison in wolfSSL) | Das, ePrint 2026/2239 | Not affected — every byte compared, rejection key selected fail-closed; `test/attack_regressions.jl` tampers every ciphertext byte |
| Secret intermediates left in memory: the ML-KEM basemul scratch (ŝ_k recoverable from one decapsulation and its ciphertext) and the decapsulation re-encryption of m′; ML-DSA per-attempt challenge buffers; FALCON-MRM targets t0, t1 (f and F from the public c) and sampler seed; the Falcon keygen seed (regenerates the key) and decoding's mod-q and NTT copies; every buffer on error paths | Internal bug sweep (2026-10) | Fixed — wiped in `finally` on every exit; `test/regressions_*.jl` |
| `wipe!` wrote into immutable Strings through `codeunits`, lost the zeroing of strided views and transposes to dead-store elimination, and cleared too few bytes for padded element types | Internal bug sweep (2026-10) | Fixed — read-only storage refused; storage pointer passed to the barrier; `elsize` bytes |
| `kyber_ek_check` accepted an invalid encapsulation key (and rejected valid ones) passed as a 0-based vector; X25519, X-Wing and the Falcon byte API misread such inputs | Internal bug sweep (2026-10) | Fixed — indexed from `firstindex` or copied to 1-based storage |
| Shamir: k ≤ 0 reconstructed the secret 0; n ≥ p put the secret in a share; a composite p was accepted; `Int` overflow in x-coordinate differences gave a wrong secret | Internal bug sweep (2026-10) | Fixed — `ArgumentError`; BigInt interpolation |
| Falcon keygen's bare `catch` around the NTRU solver swallowed `InterruptException` and every other error | Internal bug sweep (2026-10) | Fixed — only `NTRUSolveFailure` is caught |
| `test/timing/dispatch.jl` failed on Julia 1.12.1: `Base.vcat` of byte vectors builds its assertion message through dynamic dispatch | Internal bug sweep (2026-10) | Fixed — `catbytes` (src/bytes.jl); every row 0 |
| ML-DSA builds that drop a load-bearing reduction before the inverse NTT overflow Int32 and still pass KATs (wolfSSL small build) | Lee, Lim, Yoon, ePrint 2026/1032 | Not affected — all reductions present; `test/attack_regressions.jl` checks every `invntt!` input during keygen, signing and verification and compares worst-case pipelines with a schoolbook product |

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

`fips206=true` applies the NIST FIPS 206 status update (Perlner, 2025). NIST's publication
listings show no FIPS 206 text, draft or final, as of 2026-10-09, so these options may change
when the standard appears:

- the keygen GS bound 0.9999·1.17√q, decided exactly by the certificate;
- the signer's refusal of keys with a leaf outside [σmin, σmax];
- randomized signing;
- G recomputed from (f, g, F);
- the public key in NTT form, ĥ = NTT(g)/NTT(f) (`Falcon.pk_ntt`, `Falcon.verify_poly_ntt`).

## FALCON-MRM

FALCON-MRM follows ePrint 2026/420, which publishes no test vectors and leaves open the
domain-separation tags of H1 and H2 (`FalconMRM.TAG_H1`, `FalconMRM.TAG_H2`), the framing of
H1's input, the byte encoding of M2, the header byte, the sampling of ρ and the handling of an overlong
compressed signature. PQJulia's choices are listed in the
[README](README.md#message-recovery-falcon-mrm); until a revision of the specification fixes
them, its FALCON-MRM signatures do not interoperate with other implementations.

## X25519 and X-Wing

`X25519.x25519` is the RFC 7748 function: it masks the top bit of u, accepts non-canonical u,
and returns the all-zero value for small-order inputs (Wycheproof `ZeroSharedSecret`). For
standalone Diffie–Hellman, reject an all-zero result (RFC 7748 §6.1). X-Wing performs no such
check, and the draft specifies none: its security proof (Barbosa et al., CiC 2024, §7.1) models
X25519 as the RFC 7748 function on arbitrary 32-byte strings, and the combiner hashes the X25519
ciphertext and public key.

X25519 gives no protection against a quantum adversary. Shor's algorithm computes elliptic-curve
discrete logarithms, so X25519 exchanges recorded today can be decrypted once a fault-tolerant
quantum computer of sufficient size exists (harvest now, decrypt later). 2026 estimates for one
256-bit prime-field instance are about 1,200–1,450 logical qubits and 40–90 million Toffoli gates
(Babbush et al., arXiv:2603.28846; Häner et al., arXiv:2609.05625), or 835 logical qubits at
2^30.9 Toffoli gates (Luo et al., arXiv:2607.13816); Babbush et al. expect many curves with a
256-bit modulus and group order to cost the same order of magnitude. These are resource estimates for hardware that does
not exist; no such computation has been carried out, and no classical result weakens X25519.
X-Wing's post-quantum security rests on ML-KEM-768 (with SHA3-256 as a PRF); X25519 is a
classical hedge against a failure of ML-KEM.

X-Wing is an individual Internet-Draft, not a CFRG document or an RFC. It has been in the
Independent Submission stream since revision -07 (2025-05-26); -11 (2026-09-23), which PQJulia
implements, refreshed the expired -10 with no normative change.

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
- Zhou, Wang, Sun, Yu, [Every Signing Leaks: Breaking Falcon via Floating-Point Conversion Leakage](https://eprint.iacr.org/2026/2124) (ePrint 2026/2124)
- Li, Ma, Dou, Guo, [One Fell Swoop: A Single-Trace Key-Recovery Attack on the Falcon Signing Algorithm](https://eprint.iacr.org/2025/2159) (ePrint 2025/2159; TCHES 2027)
- Brinkmann, Kraus, May, [Halfspace Learning for Lattice Signature Key Recovery from Signs](https://eprint.iacr.org/2026/1366) (ePrint 2026/1366; CRYPTO 2026)
- Kaihara, Abou Haidar, Tibouchi, Abe, [Square Root of All Evil: The Dangers of Falcon's Superfluous Square Roots](https://eprint.iacr.org/2026/2046) (ePrint 2026/2046; ACM CCS 2026)
- Lin, Tibouchi, Yu, [Swing the Lure: How to Cheaply Mitigate Sign Leakage in Falcon](https://eprint.iacr.org/2026/2170) (ePrint 2026/2170)
- Barbosa et al., [X-Wing: The Hybrid KEM You've Been Looking For](https://doi.org/10.62056/a3qj89n4e) (IACR Communications in Cryptology 1(1), 2024; ePrint 2024/039)
- [X-Wing, draft-connolly-cfrg-xwing-kem-11](https://www.ietf.org/archive/id/draft-connolly-cfrg-xwing-kem-11.txt) and its [history](https://datatracker.ietf.org/doc/draft-connolly-cfrg-xwing-kem/history/)
- Babbush et al., [arXiv:2603.28846](https://arxiv.org/abs/2603.28846); Häner et al. (IonQ), [arXiv:2609.05625](https://arxiv.org/abs/2609.05625); Luo et al., [arXiv:2607.13816](https://arxiv.org/abs/2607.13816) (Shor-ECDLP resource estimates)
