using Test
import JSON

# Load the package
include(joinpath(@__DIR__, "..", "src", "PQJulia.jl"))
using .PQJulia

println("=" ^ 70)
println("  PQJulia.jl — Comprehensive Test Suite")
println("=" ^ 70)

# ==================== ML-KEM ====================

@testset "ML-KEM Roundtrip (all levels)" begin
    for (name, Cat) in [("512", MLKEM.Category1), ("768", MLKEM.Category3), ("1024", MLKEM.Category5)]
        @testset "ML-KEM-$name" begin
            for _ in 1:3
                pk, sk = Cat.kyber_kem_keypair()
                ct, ss1 = Cat.kyber_kem_enc(pk)
                ss2 = Cat.kyber_kem_dec(ct, sk)
                @test ss1 == ss2
            end
        end
    end
end

kat_dir = joinpath(@__DIR__, "kat")
if isdir(kat_dir)
    keygen_kem = JSON.parsefile(joinpath(kat_dir, "mlkem_keygen_prompt.json"))
    encap_kem = JSON.parsefile(joinpath(kat_dir, "mlkem_encapdecap_prompt.json"))

    for (ps, Cat) in [("ML-KEM-512", MLKEM.Category1), ("ML-KEM-768", MLKEM.Category3), ("ML-KEM-1024", MLKEM.Category5)]
        gkg = [g for g in keygen_kem["testGroups"] if g["parameterSet"] == ps][1]
        @testset "$ps KAT KeyGen ($(length(gkg["tests"])))" begin
            for tc in gkg["tests"]
                pk, sk = Cat.kyber_kem_keypair_derand(vcat(hex2bytes(tc["d"]), hex2bytes(tc["z"])))
                @test pk == hex2bytes(tc["ek"])
                @test sk == hex2bytes(tc["dk"])
            end
        end

        genc = [g for g in encap_kem["testGroups"] if g["parameterSet"] == ps && g["function"] == "encapsulation"][1]
        @testset "$ps KAT Encaps ($(length(genc["tests"])))" begin
            for tc in genc["tests"]
                ct, ss = Cat.kyber_kem_enc_derand(hex2bytes(tc["ek"]), hex2bytes(tc["m"]))
                @test ct == hex2bytes(tc["c"])
                @test ss == hex2bytes(tc["k"])
            end
        end

        gdec = [g for g in encap_kem["testGroups"] if g["parameterSet"] == ps && g["function"] == "decapsulation"][1]
        @testset "$ps KAT Decaps ($(length(gdec["tests"])))" begin
            for tc in gdec["tests"]
                k = Cat.kyber_kem_dec(hex2bytes(tc["c"]), hex2bytes(tc["dk"]))
                @test k == hex2bytes(tc["k"])
            end
        end
    end

    # ==================== ML-DSA ====================

    keygen_dsa = JSON.parsefile(joinpath(kat_dir, "mldsa_keygen_prompt.json"))
    siggen_dsa = JSON.parsefile(joinpath(kat_dir, "mldsa_siggen_prompt.json"))
    sigver_dsa = JSON.parsefile(joinpath(kat_dir, "mldsa_sigver_prompt.json"))

    for (ps, Cat) in [("ML-DSA-44", MLDSA.Category2), ("ML-DSA-65", MLDSA.Category3), ("ML-DSA-87", MLDSA.Category5)]
        # KeyGen KAT
        gkg = [g for g in keygen_dsa["testGroups"] if g["parameterSet"] == ps][1]
        @testset "$ps KAT KeyGen ($(length(gkg["tests"])))" begin
            for tc in gkg["tests"]
                pk, sk = Cat.dilithium_keygen_derand(hex2bytes(tc["seed"]))
                @test pk == hex2bytes(tc["pk"])
                @test sk == hex2bytes(tc["sk"])
            end
        end

        # SigGen pure KAT
        local gsg = nothing
        for g in siggen_dsa["testGroups"]
            if g["parameterSet"] == ps && g["deterministic"] == true &&
               get(g, "preHash", "") == "pure" && g["signatureInterface"] == "external"
                gsg = g; break
            end
        end
        if gsg !== nothing
            @testset "$ps KAT SigGen pure ($(length(gsg["tests"])))" begin
                for tc in gsg["tests"]
                    sig = Cat.dilithium_sign_derand(hex2bytes(tc["message"]), hex2bytes(tc["sk"]), zeros(UInt8, 32); context=hex2bytes(tc["context"]))
                    @test sig == hex2bytes(tc["signature"])
                end
            end
        end

        # SigGen preHash KAT
        local gph = nothing
        for g in siggen_dsa["testGroups"]
            if g["parameterSet"] == ps && g["deterministic"] == true &&
               get(g, "preHash", "") == "preHash" && g["signatureInterface"] == "external"
                gph = g; break
            end
        end
        if gph !== nothing
            @testset "$ps KAT SigGen preHash ($(length(gph["tests"])))" begin
                for tc in gph["tests"]
                    sig = Cat.dilithium_sign_prehash(hex2bytes(tc["message"]), hex2bytes(tc["sk"]), tc["hashAlg"]; context=hex2bytes(tc["context"]))
                    @test sig == hex2bytes(tc["signature"])
                end
            end
        end

        # SigVer pure KAT
        local gvp = nothing
        for g in sigver_dsa["testGroups"]
            if g["parameterSet"] == ps && get(g, "preHash", "") == "pure" && g["signatureInterface"] == "external"
                gvp = g; break
            end
        end
        if gvp !== nothing
            @testset "$ps KAT SigVer pure ($(length(gvp["tests"])))" begin
                for tc in gvp["tests"]
                    result = Cat.dilithium_verify(hex2bytes(tc["message"]), hex2bytes(tc["signature"]), hex2bytes(tc["pk"]); context=hex2bytes(tc["context"]))
                    @test result == tc["testPassed"]
                end
            end
        end

        # SigVer preHash KAT
        local gvph = nothing
        for g in sigver_dsa["testGroups"]
            if g["parameterSet"] == ps && get(g, "preHash", "") == "preHash" && g["signatureInterface"] == "external"
                gvph = g; break
            end
        end
        if gvph !== nothing
            @testset "$ps KAT SigVer preHash ($(length(gvph["tests"])))" begin
                for tc in gvph["tests"]
                    result = Cat.dilithium_verify_prehash(hex2bytes(tc["message"]), hex2bytes(tc["signature"]), hex2bytes(tc["pk"]), tc["hashAlg"]; context=hex2bytes(tc["context"]))
                    @test result == tc["testPassed"]
                end
            end
        end
    end
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
end

println("\n" * "=" ^ 70)
println("  All PQJulia.jl tests complete!")
println("=" ^ 70)
