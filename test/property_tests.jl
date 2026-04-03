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
    # Shamir: below-threshold gives wrong answer
    # ═══════════════════════════════════════════════════════════════
    @testset "Shamir below-threshold ($NTRIALS trials)" begin
        failures = 0
        for _ in 1:NTRIALS
            secret = parse(BigInt, bytes2hex(rand(UInt8, 32)), base=16)
            shares = shamir_share(secret, 3, 5)
            wrong = shamir_reconstruct([shares[1], shares[2]], 2)
            if wrong != secret
                failures += 1
            end
        end
        # With overwhelming probability, k-1 shares should not reconstruct
        @test failures >= NTRIALS * 0.99
    end
end

println("\n" * "=" ^ 70)
println("ALL PROPERTY TESTS PASSED")
println("=" ^ 70)
