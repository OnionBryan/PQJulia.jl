#!/usr/bin/env julia
"""
PQJulia.jl Property-Based Tests
================================
Random-input tests that verify structural properties hold beyond KAT vectors.
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PQJulia
using Test
using Random
using SHA

Random.seed!(42)
const NTRIALS = 100

println("=" ^ 70)
println("PQJulia.jl Property-Based Tests ($NTRIALS trials each)")
println("=" ^ 70)

@testset "Property Tests" begin

    # ═══════════════════════════════════════════════════════════════
    # ML-KEM: encaps/decaps round-trip for random keypairs
    # ═══════════════════════════════════════════════════════════════
    @testset "ML-KEM round-trip ($NTRIALS trials)" begin
        for (name, Cat) in [("512", MLKEM.Category1), ("768", MLKEM.Category3), ("1024", MLKEM.Category5)]
            @testset "ML-KEM-$name" begin
                for _ in 1:NTRIALS
                    pk, sk = Cat.kyber_kem_keypair()
                    ct, ss1 = Cat.kyber_kem_enc(pk)
                    ss2 = Cat.kyber_kem_dec(ct, sk)
                    @test ss1 == ss2
                    @test length(ss1) == 32
                end
            end
        end
    end

    # ═══════════════════════════════════════════════════════════════
    # ML-KEM: wrong sk must NOT produce matching shared secret
    # ═══════════════════════════════════════════════════════════════
    @testset "ML-KEM wrong-key rejection ($NTRIALS trials)" begin
        for _ in 1:NTRIALS
            pk1, sk1 = MLKEM.Category3.kyber_kem_keypair()
            _, sk2 = MLKEM.Category3.kyber_kem_keypair()
            ct, ss1 = MLKEM.Category3.kyber_kem_enc(pk1)
            ss_wrong = MLKEM.Category3.kyber_kem_dec(ct, sk2)
            @test ss1 != ss_wrong  # implicit rejection returns random-looking ss
        end
    end

    # ═══════════════════════════════════════════════════════════════
    # ML-DSA: sign/verify for random messages
    # ═══════════════════════════════════════════════════════════════
    @testset "ML-DSA sign/verify random messages ($NTRIALS trials)" begin
        for (name, Cat) in [("44", MLDSA.Category2), ("65", MLDSA.Category3), ("87", MLDSA.Category5)]
            @testset "ML-DSA-$name" begin
                pk, sk = Cat.dilithium_keygen()
                for _ in 1:NTRIALS
                    msg = rand(UInt8, rand(0:1000))
                    sig = Cat.dilithium_sign(msg, sk)
                    @test Cat.dilithium_verify(msg, sig, pk)
                end
            end
        end
    end

    # ═══════════════════════════════════════════════════════════════
    # ML-DSA: tampered message must NOT verify
    # ═══════════════════════════════════════════════════════════════
    @testset "ML-DSA tamper rejection ($NTRIALS trials)" begin
        pk, sk = MLDSA.Category3.dilithium_keygen()
        for _ in 1:NTRIALS
            msg = rand(UInt8, rand(1:500))
            sig = MLDSA.Category3.dilithium_sign(msg, sk)
            tampered = copy(msg)
            tampered[rand(1:length(tampered))] ⊻= rand(UInt8(1):UInt8(255))
            @test !MLDSA.Category3.dilithium_verify(tampered, sig, pk)
        end
    end

    # ═══════════════════════════════════════════════════════════════
    # ML-DSA: wrong key must NOT verify
    # ═══════════════════════════════════════════════════════════════
    @testset "ML-DSA wrong-key rejection ($NTRIALS trials)" begin
        for _ in 1:NTRIALS
            pk1, sk1 = MLDSA.Category3.dilithium_keygen()
            pk2, _ = MLDSA.Category3.dilithium_keygen()
            msg = rand(UInt8, 64)
            sig = MLDSA.Category3.dilithium_sign(msg, sk1)
            @test !MLDSA.Category3.dilithium_verify(msg, sig, pk2)
        end
    end

    # ═══════════════════════════════════════════════════════════════
    # ML-DSA: empty message
    # ═══════════════════════════════════════════════════════════════
    @testset "ML-DSA empty message" begin
        pk, sk = MLDSA.Category3.dilithium_keygen()
        sig = MLDSA.Category3.dilithium_sign(UInt8[], sk)
        @test MLDSA.Category3.dilithium_verify(UInt8[], sig, pk)
    end

    # ═══════════════════════════════════════════════════════════════
    # Shamir: round-trip for random 256-bit secrets
    # ═══════════════════════════════════════════════════════════════
    @testset "Shamir round-trip ($NTRIALS trials)" begin
        for _ in 1:NTRIALS
            secret = parse(BigInt, bytes2hex(rand(UInt8, 32)), base=16)
            shares = shamir_share(secret, 3, 5)
            # Pick 3 random shares
            idx = sort(shuffle(1:5)[1:3])
            recovered = shamir_reconstruct([shares[i] for i in idx], 3)
            @test recovered == secret
        end
    end

    # ═══════════════════════════════════════════════════════════════
    # Shamir: below-threshold MUST ALWAYS give wrong answer
    # (p = 2^521-1, probability of accidental match ≈ 2^-521)
    # ═══════════════════════════════════════════════════════════════
    @testset "Shamir below-threshold ($NTRIALS trials)" begin
        for _ in 1:NTRIALS
            secret = parse(BigInt, bytes2hex(rand(UInt8, 32)), base=16)
            shares = shamir_share(secret, 3, 5)
            wrong = shamir_reconstruct([shares[1], shares[2]], 2)
            @test wrong != secret  # must ALWAYS differ — any match is a catastrophic bug
        end
    end
end

println("\n" * "=" ^ 70)
println("ALL PROPERTY TESTS PASSED")
println("=" ^ 70)

# ═══════════════════════════════════════════════════════════════════════
# NEGATIVE / ADVERSARIAL TESTS
# Based on C2SP/CCTV, Wycheproof ML-KEM/ML-DSA, and mlkem-native
# ═══════════════════════════════════════════════════════════════════════

@testset "Negative & Edge Case Tests" begin

    # ─── ML-KEM: Invalid ciphertext → implicit rejection ─────────
    @testset "ML-KEM implicit rejection (corrupted ciphertext)" begin
        for (name, Cat) in [("512", MLKEM.Category1), ("768", MLKEM.Category3), ("1024", MLKEM.Category5)]
            @testset "ML-KEM-$name" begin
                pk, sk = Cat.kyber_kem_keypair()
                ct, ss_good = Cat.kyber_kem_enc(pk)

                for _ in 1:10
                    # Flip a random byte in the ciphertext
                    ct_bad = copy(ct)
                    pos = rand(1:length(ct_bad))
                    ct_bad[pos] ⊻= rand(UInt8(1):UInt8(255))

                    ss_bad = Cat.kyber_kem_dec(ct_bad, sk)
                    # Must NOT return the real shared secret (implicit rejection)
                    @test ss_bad != ss_good
                    # Must still return 32 bytes (not crash)
                    @test length(ss_bad) == 32
                end
            end
        end
    end

    # ─── ML-KEM: Completely random ciphertext → still returns 32 bytes ───
    @testset "ML-KEM random ciphertext decaps" begin
        pk, sk = MLKEM.Category3.kyber_kem_keypair()
        for _ in 1:20
            ct_rand = rand(UInt8, 1088)  # ML-KEM-768 ciphertext size
            ss = MLKEM.Category3.kyber_kem_dec(ct_rand, sk)
            @test length(ss) == 32  # Must not crash, must return pseudorandom
        end
    end

    # ─── ML-DSA: Bit-flipped signature → must reject ─────────────
    @testset "ML-DSA bit-flipped signature rejection" begin
        pk, sk = MLDSA.Category3.dilithium_keygen()
        msg = rand(UInt8, 100)
        sig = MLDSA.Category3.dilithium_sign(msg, sk)

        for _ in 1:20
            sig_bad = copy(sig)
            pos = rand(1:length(sig_bad))
            sig_bad[pos] ⊻= UInt8(1 << rand(0:7))  # flip one bit
            @test !MLDSA.Category3.dilithium_verify(msg, sig_bad, pk)
        end
    end

    # ─── ML-DSA: Wrong-length signature → must reject (not crash) ─
    @testset "ML-DSA wrong-length signature" begin
        pk, sk = MLDSA.Category3.dilithium_keygen()
        msg = Vector{UInt8}("test")
        sig = MLDSA.Category3.dilithium_sign(msg, sk)

        # Too short
        @test !MLDSA.Category3.dilithium_verify(msg, sig[1:end-1], pk)
        # Too long
        @test !MLDSA.Category3.dilithium_verify(msg, vcat(sig, UInt8[0x00]), pk)
        # Empty
        @test !MLDSA.Category3.dilithium_verify(msg, UInt8[], pk)
    end

    # ─── ML-DSA: Overlong context → must error ───────────────────
    @testset "ML-DSA context length validation" begin
        pk, sk = MLDSA.Category3.dilithium_keygen()
        msg = Vector{UInt8}("test")
        long_ctx = rand(UInt8, 256)  # >255 bytes

        @test_throws ErrorException MLDSA.Category3.dilithium_sign(msg, sk; context=long_ctx)
        sig = MLDSA.Category3.dilithium_sign(msg, sk)
        @test_throws ErrorException MLDSA.Category3.dilithium_verify(msg, sig, pk; context=long_ctx)

        # 255 bytes should be fine
        ok_ctx = rand(UInt8, 255)
        sig_ctx = MLDSA.Category3.dilithium_sign(msg, sk; context=ok_ctx)
        @test MLDSA.Category3.dilithium_verify(msg, sig_ctx, pk; context=ok_ctx)
    end

    # ─── ML-DSA: Context mismatch → must reject ──────────────────
    @testset "ML-DSA context mismatch rejection" begin
        pk, sk = MLDSA.Category3.dilithium_keygen()
        msg = Vector{UInt8}("test")
        ctx1 = Vector{UInt8}("context-A")
        ctx2 = Vector{UInt8}("context-B")

        sig = MLDSA.Category3.dilithium_sign(msg, sk; context=ctx1)
        @test MLDSA.Category3.dilithium_verify(msg, sig, pk; context=ctx1)
        @test !MLDSA.Category3.dilithium_verify(msg, sig, pk; context=ctx2)
        @test !MLDSA.Category3.dilithium_verify(msg, sig, pk)  # empty context
    end

    # ─── Shamir: Edge cases ──────────────────────────────────────
    @testset "Shamir edge cases" begin
        # k > n
        @test_throws ArgumentError shamir_share(big(42), 5, 3)

        # Secret = 0 (boundary)
        shares = shamir_share(big(0), 2, 3)
        @test shamir_reconstruct([shares[1], shares[3]], 2) == big(0)

        # Secret = 1 (boundary)
        shares = shamir_share(big(1), 2, 3)
        @test shamir_reconstruct([shares[2], shares[3]], 2) == big(1)

        # Large secret (512-bit)
        big_secret = parse(BigInt, bytes2hex(rand(UInt8, 64)), base=16)
        shares = shamir_share(big_secret, 3, 5)
        @test shamir_reconstruct([shares[1], shares[2], shares[5]], 3) == big_secret
    end
end

println("\n" * "=" ^ 70)
println("ALL PROPERTY + NEGATIVE TESTS PASSED")
println("=" ^ 70)

# ═══════════════════════════════════════════════════════════════════════
# FIXES FROM AUDIT ROUND (CCTV + Wycheproof + mlkem-native + Trail of Bits)
# ═══════════════════════════════════════════════════════════════════════

@testset "Audit-Driven Tests" begin

    # ─── Trail of Bits #1: Implicit rejection must be DETERMINISTIC ──
    @testset "ML-KEM implicit rejection determinism + formula" begin
        for (name, Cat) in [("512", MLKEM.Category1), ("768", MLKEM.Category3), ("1024", MLKEM.Category5)]
            @testset "ML-KEM-$name" begin
                pk, sk = Cat.kyber_kem_keypair()
                ct, ss_good = Cat.kyber_kem_enc(pk)
                ct_bad = copy(ct); ct_bad[1] ⊻= 0x01

                ss1 = Cat.kyber_kem_dec(ct_bad, sk)
                ss2 = Cat.kyber_kem_dec(ct_bad, sk)
                @test ss1 == ss2              # deterministic
                @test ss1 != ss_good          # not the real key
                @test length(ss1) == 32

                # Verify the actual FIPS 203 formula: ss = SHAKE256(z || ct)
                # z is the last 32 bytes of the KEM secret key
                z = sk[end-31:end]
                expected = SHA.shake256(vcat(z, ct_bad), UInt64(32))
                @test ss1 == expected         # must equal rkprf(z, ct_bad)
            end
        end
    end

    # ─── Trail of Bits #3: Encaps must be NON-deterministic ──────────
    @testset "ML-KEM encaps non-determinism" begin
        pk, _ = MLKEM.Category3.kyber_kem_keypair()
        ct1, ss1 = MLKEM.Category3.kyber_kem_enc(pk)
        ct2, ss2 = MLKEM.Category3.kyber_kem_enc(pk)
        @test ct1 != ct2
        @test ss1 != ss2
    end

    # ─── Trail of Bits #4: Deterministic signing is idempotent ───────
    @testset "ML-DSA deterministic signing idempotency" begin
        for (name, Cat) in [("44", MLDSA.Category2), ("65", MLDSA.Category3), ("87", MLDSA.Category5)]
            @testset "ML-DSA-$name" begin
                pk, sk = Cat.dilithium_keygen()
                msg = rand(UInt8, 64)
                sig1 = Cat.dilithium_sign(msg, sk)  # hedged=false default
                sig2 = Cat.dilithium_sign(msg, sk)
                @test sig1 == sig2
            end
        end
    end

    # ─── Trail of Bits #2: Shamir below-threshold must ALWAYS fail ───
    @testset "Shamir below-threshold is 100% failure" begin
        for _ in 1:100
            secret = parse(BigInt, bytes2hex(rand(UInt8, 32)), base=16)
            shares = shamir_share(secret, 3, 5)
            wrong = shamir_reconstruct([shares[1], shares[2]], 2)
            @test wrong != secret
        end
    end

    # ─── Wycheproof #2: Context sign-empty / verify-nonempty ─────────
    @testset "ML-DSA context injection rejection" begin
        pk, sk = MLDSA.Category3.dilithium_keygen()
        msg = Vector{UInt8}("test")
        ctx = Vector{UInt8}("bound-context")

        sig_noctx = MLDSA.Category3.dilithium_sign(msg, sk)  # empty context
        @test !MLDSA.Category3.dilithium_verify(msg, sig_noctx, pk; context=ctx)

        sig_ctx = MLDSA.Category3.dilithium_sign(msg, sk; context=ctx)
        @test !MLDSA.Category3.dilithium_verify(msg, sig_ctx, pk)  # empty verify
    end

    # ─── Wycheproof #3: Prehash invalid hash algorithm ──────────────
    @testset "ML-DSA prehash invalid hash algorithm" begin
        pk, sk = MLDSA.Category3.dilithium_keygen()
        msg = Vector{UInt8}("test")
        @test_throws ErrorException MLDSA.Category3.dilithium_sign_prehash(msg, sk, "SHA2-999")
        sig = MLDSA.Category3.dilithium_sign_prehash(msg, sk, "SHA2-256")
        @test_throws ErrorException MLDSA.Category3.dilithium_verify_prehash(msg, sig, pk, "SHA2-999")
    end

    # ─── Wycheproof #1: Wrong-length sig across all levels ───────────
    @testset "ML-DSA wrong-length signature (all levels)" begin
        for (name, Cat) in [("44", MLDSA.Category2), ("65", MLDSA.Category3), ("87", MLDSA.Category5)]
            @testset "ML-DSA-$name" begin
                pk, sk = Cat.dilithium_keygen()
                msg = Vector{UInt8}("test")
                sig = Cat.dilithium_sign(msg, sk)
                @test !Cat.dilithium_verify(msg, sig[1:end-1], pk)      # too short
                @test !Cat.dilithium_verify(msg, vcat(sig, UInt8[0]), pk) # too long
                @test !Cat.dilithium_verify(msg, UInt8[], pk)            # empty
            end
        end
    end
end
