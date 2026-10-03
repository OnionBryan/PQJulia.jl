using Test
using PQJulia
import JSON

println("=" ^ 70)
println("  PQJulia.jl — Test Suite")
println("=" ^ 70)

const KAT_DIR = joinpath(@__DIR__, "kat")
isdir(KAT_DIR) || error("KAT directory not found at $KAT_DIR — the ACVP vectors are required")
h(x) = hex2bytes(x)
groups(file, pred) = [g for g in JSON.parsefile(joinpath(KAT_DIR, file))["testGroups"] if pred(g)]

# ==================== ML-KEM (FIPS 203) ====================

@testset "ML-KEM Roundtrip (all levels)" begin
    for (name, Cat) in [("512", MLKEM.Category1), ("768", MLKEM.Category3), ("1024", MLKEM.Category5)]
        @testset "ML-KEM-$name" begin
            for _ in 1:3
                pk, sk = Cat.kyber_kem_keypair()
                ct, ss1 = Cat.kyber_kem_enc(pk)
                @test ss1 == Cat.kyber_kem_dec(ct, sk)
            end
        end
    end
end

for (ps, Cat) in [("ML-KEM-512", MLKEM.Category1), ("ML-KEM-768", MLKEM.Category3), ("ML-KEM-1024", MLKEM.Category5)]
    @testset "$ps ACVP" begin
        for g in groups("mlkem_keygen_prompt.json", g -> g["parameterSet"] == ps)
            @testset "KeyGen ($(length(g["tests"])))" begin
                for tc in g["tests"]
                    pk, sk = Cat.kyber_kem_keypair_derand(vcat(h(tc["d"]), h(tc["z"])))
                    @test pk == h(tc["ek"])
                    @test sk == h(tc["dk"])
                end
            end
        end
        for g in groups("mlkem_encapdecap_prompt.json", g -> g["parameterSet"] == ps)
            @testset "$(g["function"]) ($(length(g["tests"])))" begin
                for tc in g["tests"]
                    if g["function"] == "encapsulation"
                        ct, ss = Cat.kyber_kem_enc_derand(h(tc["ek"]), h(tc["m"]))
                        @test ct == h(tc["c"])
                        @test ss == h(tc["k"])
                    elseif g["function"] == "decapsulation"
                        @test Cat.kyber_kem_dec(h(tc["c"]), h(tc["dk"])) == h(tc["k"])
                    elseif g["function"] == "encapsulationKeyCheck"
                        @test Cat.kyber_ek_check(h(tc["ek"])) == tc["testPassed"]
                    elseif g["function"] == "decapsulationKeyCheck"
                        @test Cat.kyber_dk_check(h(tc["dk"])) == tc["testPassed"]
                    end
                end
            end
        end
    end
end

@testset "ML-KEM input validation" begin
    Cat = MLKEM.Category3
    pk, sk = Cat.kyber_kem_keypair()
    ct, _ = Cat.kyber_kem_enc(pk)
    bad_pk = copy(pk); bad_pk[1] = 0xff; bad_pk[2] |= 0x0f               # first coefficient = 4095 ≥ q
    @test !Cat.kyber_ek_check(bad_pk)
    @test_throws ArgumentError Cat.kyber_kem_enc(bad_pk)
    @test_throws ArgumentError Cat.kyber_kem_enc(pk[1:end-1])
    bad_sk = copy(sk); bad_sk[end-40] ⊻= 0x01                            # corrupt H(ek)
    @test_throws ArgumentError Cat.kyber_kem_dec(ct, bad_sk)
    @test_throws ArgumentError Cat.kyber_kem_dec(ct[1:end-1], sk)
end

# ==================== ML-DSA (FIPS 204) ====================

for (ps, Cat) in [("ML-DSA-44", MLDSA.Category2), ("ML-DSA-65", MLDSA.Category3), ("ML-DSA-87", MLDSA.Category5)]
    @testset "$ps ACVP" begin
        for g in groups("mldsa_keygen_prompt.json", g -> g["parameterSet"] == ps)
            @testset "KeyGen ($(length(g["tests"])))" begin
                for tc in g["tests"]
                    pk, sk = Cat.dilithium_keygen_derand(h(tc["seed"]))
                    @test pk == h(tc["pk"])
                    @test sk == h(tc["sk"])
                end
            end
        end

        for g in groups("mldsa_siggen_prompt.json", g -> g["parameterSet"] == ps)
            label = "SigGen $(g["signatureInterface"])/$(g["preHash"])$(g["externalMu"] ? "/μ" : "")" *
                    " $(g["deterministic"] ? "det" : "hedged") ($(length(g["tests"])))"
            @testset "$label" begin
                for tc in g["tests"]
                    rnd = g["deterministic"] ? zeros(UInt8, 32) : h(tc["rnd"])
                    sig = if g["signatureInterface"] == "internal"
                        g["externalMu"] ? Cat.dilithium_sign_internal_mu(h(tc["mu"]), h(tc["sk"]), rnd) :
                                          Cat.dilithium_sign_internal(h(tc["message"]), h(tc["sk"]), rnd)
                    elseif g["preHash"] == "preHash"
                        Cat.dilithium_sign_prehash_derand(h(tc["message"]), h(tc["sk"]), tc["hashAlg"], rnd; context=h(tc["context"]))
                    else
                        Cat.dilithium_sign_derand(h(tc["message"]), h(tc["sk"]), rnd; context=h(tc["context"]))
                    end
                    @test sig == h(tc["signature"])
                end
            end
        end

        for g in groups("mldsa_sigver_prompt.json", g -> g["parameterSet"] == ps)
            label = "SigVer $(g["signatureInterface"])/$(g["preHash"])$(g["externalMu"] ? "/μ" : "") ($(length(g["tests"])))"
            @testset "$label" begin
                for tc in g["tests"]
                    ok = if g["signatureInterface"] == "internal"
                        g["externalMu"] ? Cat.dilithium_verify_mu(h(tc["mu"]), h(tc["signature"]), h(tc["pk"])) :
                                          Cat.dilithium_verify_internal(h(tc["message"]), h(tc["signature"]), h(tc["pk"]))
                    elseif g["preHash"] == "preHash"
                        Cat.dilithium_verify_prehash(h(tc["message"]), h(tc["signature"]), h(tc["pk"]), tc["hashAlg"]; context=h(tc["context"]))
                    else
                        Cat.dilithium_verify(h(tc["message"]), h(tc["signature"]), h(tc["pk"]); context=h(tc["context"]))
                    end
                    @test ok == tc["testPassed"]
                end
            end
        end
    end
end

@testset "ML-DSA hedged default + interfaces agree" begin
    Cat = MLDSA.Category2
    pk, sk = Cat.dilithium_keygen()
    msg = Vector{UInt8}("interfaces"); ctx = Vector{UInt8}("ctx")
    s1 = Cat.dilithium_sign(msg, sk; context=ctx); s2 = Cat.dilithium_sign(msg, sk; context=ctx)
    @test s1 != s2                                                    # hedged: fresh rnd each call
    @test Cat.dilithium_verify(msg, s1, pk; context=ctx) && Cat.dilithium_verify(msg, s2, pk; context=ctx)
    @test Cat.dilithium_sign(msg, sk; hedged=false) == Cat.dilithium_sign(msg, sk; hedged=false)
    mprime = vcat(UInt8[0x00, UInt8(length(ctx))], ctx, msg)
    @test Cat.dilithium_sign_derand(msg, sk, zeros(UInt8, 32); context=ctx) ==
          Cat.dilithium_sign_internal(mprime, sk, zeros(UInt8, 32))
    @test Cat.dilithium_verify_internal(mprime, s1, pk)
    @test !Cat.dilithium_verify(msg, s1, pk)                         # wrong context
    @test !Cat.dilithium_verify(msg, s1, pk[1:end-1])                # malformed pk → false, not an exception

    # Regression: MakeHint must see HighBits(w), not w. This (key, message) hits w0 = −γ2 with
    # HighBits = 0, where the pre-0.2 signer emitted a signature that did not verify.
    pk, sk = Cat.dilithium_keygen_derand(zeros(UInt8, 32))
    m = collect(reinterpret(UInt8, [18297]))
    @test Cat.dilithium_verify(m, Cat.dilithium_sign_derand(m, sk, zeros(UInt8, 32)), pk)
end

# ==================== FN-DSA (Falcon) ====================

include("falcon_kat.jl")
include("falcon_certify.jl")
include("falcon_mrm.jl")
include("falcon_fxp.jl")
include("vectors_extra.jl")

# ==================== X25519 and X-Wing ====================

include("x25519_xwing.jl")
include("wipe.jl")

@testset "FN-DSA API ($(M.IDENTIFIER))" for M in (FNDSA.Falcon512, FNDSA.Falcon1024)
    pk, sk = M.falcon_keygen()
    @test (length(pk), length(sk)) == (M.PK_BYTES, M.SK_BYTES)
    msg = Vector{UInt8}("FN-DSA")
    sig = M.falcon_sign(msg, sk)
    @test length(sig) == M.SIG_BYTES
    @test M.falcon_verify(msg, sig, pk)
    ek = M.falcon_expand_sk(sk)
    @test all(M.falcon_verify(UInt8[msg; i], M.falcon_sign(UInt8[msg; i], ek), pk) for i in 0x00:0x04)
    @test !M.falcon_verify(UInt8[msg; 0x00], sig, pk)
    bad = copy(sig); bad[60] ⊻= 0x01
    @test !M.falcon_verify(msg, bad, pk)
    @test !M.falcon_verify(msg, sig[1:end-1], pk)
    @test !M.falcon_verify(msg, sig, pk[1:end-1])
    pk2, _ = M.falcon_keygen()
    @test !M.falcon_verify(msg, sig, pk2)
    @test_throws ArgumentError M.falcon_sign(msg, sk[1:end-1])
end

# ==================== Shamir ====================

@testset "Shamir Secret Sharing" begin
    @testset "Basic roundtrip" begin
        for secret in [0, 1, 42, 1000, big(2)^126]
            shares = shamir_share(secret, 3, 5)
            @test shamir_reconstruct(shares, 3) == secret
        end
    end
    @testset "Any k-subset" begin
        shares = shamir_share(big(123456789), 3, 7)
        for combo in [[1,2,3], [1,4,7], [3,5,7]]
            @test shamir_reconstruct([shares[i] for i in combo], 3) == big(123456789)
        end
    end
    @testset "Byte-level API" begin
        secret_bytes = rand(UInt8, 32)
        shares = shamir_share_bytes(secret_bytes, 3, 5)
        @test shamir_reconstruct_bytes(shares, 3, 32) == secret_bytes
    end
end


# ==================== Property tests ====================

include("property_tests.jl")

println("\n" * "=" ^ 70)
println("  All PQJulia.jl tests complete!")
println("=" ^ 70)
