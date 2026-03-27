# PQJulia.jl

Julia implementations of NIST post-quantum cryptographic standards.

This is a doctoral research project implementing the three FIPS post-quantum standards in pure Julia.

## What's Implemented

| Standard | Algorithm | Levels |
|----------|-----------|--------|
| FIPS 203 | ML-KEM (Kyber) | 512 / 768 / 1024 |
| FIPS 204 | ML-DSA (Dilithium) | 44 / 65 / 87 |

Also includes Shamir (k,n)-threshold secret sharing.

Validated against NIST ACVP test vectors where available. ML-KEM and ML-DSA pass 750+ known answer tests. ML-KEM has been cross-validated against the pq-crystals C reference implementation.

## Usage

```julia
using PQJulia

# Key Encapsulation (ML-KEM-768)
pk, sk = MLKEM.Category3.kyber_kem_keypair()
ct, shared_secret = MLKEM.Category3.kyber_kem_enc(pk)
shared_secret2 = MLKEM.Category3.kyber_kem_dec(ct, sk)

# Digital Signatures (ML-DSA-65)
pk, sk = MLDSA.Category3.dilithium_keygen()
sig = MLDSA.Category3.dilithium_sign(msg, sk)
MLDSA.Category3.dilithium_verify(msg, sig, pk)

# Secret Sharing
shares = shamir_share(secret, 3, 5)
secret = shamir_reconstruct(shares[1:3], 3)
```

## Security Levels

| Level | ML-KEM | ML-DSA |
|-------|--------|--------|---------|
| 1 | Category1 (512) | Category2 (44) |
| 3 | Category3 (768) | Category3 (65) |
| 5 | Category5 (1024) | Category5 (87) |

## Notice

This is research software. It has not been audited for production use.

## References

- [FIPS 203](https://csrc.nist.gov/pubs/fips/203/final)
- [FIPS 204](https://csrc.nist.gov/pubs/fips/204/final)

## License

MIT
