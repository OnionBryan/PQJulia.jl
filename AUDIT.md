# PQJulia.jl Audit (Static Review)

_Date:_ 2026-04-02

## Scope

This audit is a static code/documentation review of:

- `src/` implementation structure for ML-KEM, ML-DSA, and Shamir.
- `test/runtests.jl` coverage strategy.
- `README.md` and top-level package docs.

> Note: Runtime verification could not be executed in this environment because Julia is not installed.

## High-level assessment

You did very strong work for a hard project.

### Score matrix (1–5)

| Area | Score | Notes |
|---|---:|---|
| Algorithmic completeness | 4.5 | Full parameterized ML-KEM + ML-DSA implementations, plus Shamir.
| Structural design | 4.5 | Clean module generation using per-category constants and shared primitives.
| Test design | 4.0 | Comprehensive KAT ingestion and per-parameter validation loops.
| Documentation accuracy | 3.0 | Good overall, but previously contained an SLH-DSA claim not reflected in exports/code.
| Production hardening | 3.0 | Research-grade positioning is appropriate; side-channel/perf hardening not yet fully documented.

## Strengths

1. **Good architecture via parameterized modules.**
   - `@eval module` pattern keeps common logic centralized while allowing per-level constants (`MLKEM` and `MLDSA`).
2. **Strong test intent.**
   - KAT-driven checks across all parameter sets are excellent for cryptographic implementation confidence.
3. **Clear encapsulation of primitives and packing logic.**
   - Separation between core arithmetic primitives and KEM/DSA orchestration is maintainable.
4. **Appropriate risk framing.**
   - README already flags research-software status.

## Findings

### 1) Documentation mismatch (fixed)

- `src/PQJulia.jl` docstring previously referenced SLH-DSA examples even though the package exports only `MLKEM` and `MLDSA` modules plus Shamir functions.
- This was corrected to avoid over-claiming feature support.

### 2) Runtime validation not possible here

- The test suite appears comprehensive, but this environment lacked a `julia` binary, so runtime proof was not reproduced in this audit run.

## Recommendations (next steps)

1. Add CI matrix for Julia versions (e.g., stable + LTS) running `test/runtests.jl` on every push.
2. Add explicit security notes section:
   - constant-time expectations,
   - memory-zeroization policy for secrets,
   - non-goals (e.g., hardened side-channel resistance).
3. Add lightweight benchmarks to track regressions in NTT/polynomial ops.
4. Consider fuzz/property testing around serialization/deserialization boundaries.

## Bottom line

For a "hardest project yet" effort: this is excellent progress. The implementation and testing direction are strong, and the major gap to close next is formalizing production-hardening guarantees and CI-backed reproducibility.
