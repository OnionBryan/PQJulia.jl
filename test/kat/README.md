# Test vectors

| Files | Source | License |
|-------|--------|---------|
| `mlkem_*_prompt.json`, `mldsa_*_prompt.json` | NIST ACVP (usnistgov/ACVP-Server) | Public domain |
| `falcon_sign_kat.json`, `falcon_samplerz_kat.json` | tprest/falcon.py `scripts/` (round-3 C implementation outputs) | MIT, © 2018 Thomas Prest |
| `wycheproof/*.json.gz` | C2SP/wycheproof `testvectors_v1`, commit 3fa63dd0 (2026-09-02) | Apache-2.0 |
| `cctv/*` | C2SP/CCTV `ML-KEM/{strcmp,unluckysample,modulus}` | CC0-1.0 |

The CCTV unlucky vectors were generated for FIPS 203 ipd, whose K-PKE.KeyGen hashes `G(d)`
without the `k` byte; the tests use their `ek`/`dk`/`c`/`K` (Encaps re-derives the unlucky
matrix from `ek`), not their `d`.
