# test/vectors_extra.jl — third-party edge-case suites (test/kat/README.md has provenance):
# Wycheproof ML-DSA/ML-KEM, CCTV ML-KEM strcmp/unlucky/modulus, CCTV accumulated ML-DSA.
using Test, JSON, SHA, CodecZlib
using PQJulia
isdefined(Main, :ShakeReader) || include("shake_reader.jl")

const WP = joinpath(@__DIR__, "kat", "wycheproof")
const CV = joinpath(@__DIR__, "kat", "cctv")
hx(x) = hex2bytes(x)
hx(::Nothing) = nothing
wp(name) = JSON.parse(read(GzipDecompressorStream(open(joinpath(WP, name * ".json.gz"))), String))["testGroups"]
throws(f) = try f(); false catch e; e isa ArgumentError || e isa ErrorException end

@testset "Wycheproof ML-DSA-$lv" for (lv, C) in [(44, MLDSA.Category2), (65, MLDSA.Category3), (87, MLDSA.Category5)]
    @testset "verify" begin
        for g in wp("mldsa_$(lv)_verify_test"), t in g["tests"]
            ok = C.dilithium_verify(hx(t["msg"]), hx(t["sig"]), hx(g["publicKey"]); context=hx(get(t, "ctx", "")))
            @test ok == (t["result"] == "valid")
        end
    end
    @testset "sign ($kind)" for kind in ("seed", "noseed")
        for g in wp("mldsa_$(lv)_sign_$(kind)_test")
            sk = if kind == "seed"
                s = hx(g["privateSeed"]); s !== nothing && length(s) == 32 ? C.dilithium_keygen_derand(s)[2] : nothing
            else
                hx(g["privateKey"])
            end
            for t in g["tests"]
                ctx = hx(get(t, "ctx", ""))
                sign() = haskey(t, "msg") ? C.dilithium_sign_derand(hx(t["msg"]), sk, zeros(UInt8, 32); context=ctx) :
                                            C.dilithium_sign_internal_mu(hx(t["mu"]), sk, zeros(UInt8, 32))
                if t["result"] == "invalid"
                    @test sk === nothing || throws(sign)
                elseif "Randomized" in t["flags"]
                    @test C.dilithium_verify(hx(t["msg"]), hx(t["sig"]), hx(g["publicKey"]); context=ctx)
                else
                    @test sign() == hx(t["sig"])
                end
            end
        end
    end
end

@testset "Wycheproof ML-KEM-$lv" for (lv, C) in [(512, MLKEM.Category1), (768, MLKEM.Category3), (1024, MLKEM.Category5)]
    @testset "decaps" begin
        for g in wp("mlkem_$(lv)_test"), t in g["tests"]
            run() = (kp = C.kyber_kem_keypair_derand(hx(t["seed"])); kp[1] == hx(t["ek"]) && C.kyber_kem_dec(hx(t["c"]), kp[2]) == hx(t["K"]))
            @test t["result"] == "valid" ? run() : (try !run() catch; true end)
        end
    end
    @testset "encaps" begin
        for g in wp("mlkem_$(lv)_encaps_test"), t in g["tests"]
            enc() = C.kyber_kem_enc_derand(hx(t["ek"]), hx(t["m"]))
            @test t["result"] == "valid" ? enc() == (hx(t["c"]), hx(t["K"])) : throws(enc)
        end
    end
    @testset "keygen" begin
        for g in wp("mlkem_$(lv)_keygen_seed_test"), t in g["tests"]
            kg() = C.kyber_kem_keypair_derand(hx(t["seed"]))
            @test t["result"] == "valid" ? kg() == (hx(t["ek"]), hx(t["dk"])) : throws(kg)
        end
    end
    @testset "decaps (expanded dk)" begin
        for g in wp("mlkem_$(lv)_semi_expanded_decaps_test"), t in g["tests"]
            dec() = C.kyber_kem_dec(hx(t["c"]), hx(t["dk"]))
            @test t["result"] == "valid" ? dec() == hx(t["K"]) : throws(dec)
        end
    end
end

kv(file) = Dict(m[1] => m[2] for m in eachmatch(r"^(\w+) = ([0-9a-f]+)$"m, read(file, String)))
@testset "CCTV ML-KEM-$lv" for (lv, C) in [(512, MLKEM.Category1), (768, MLKEM.Category3), (1024, MLKEM.Category5)]
    u = kv(joinpath(CV, "unlucky-ML-KEM-$lv.txt"))             # Encaps re-samples the unlucky Â from ek
    @test C.kyber_kem_enc_derand(hx(u["ek"]), hx(u["m"])) == (hx(u["c"]), hx(u["K"]))
    @test C.kyber_kem_dec(hx(u["c"]), hx(u["dk"])) == hx(u["K"])
    s = kv(joinpath(CV, "strcmp-ML-KEM-$lv.txt"))
    @test C.kyber_kem_dec(hx(s["c"]), hx(s["dk"])) == hx(s["K"])
    bad = [hx(l) for l in eachline(GzipDecompressorStream(open(joinpath(CV, "modulus-ML-KEM-$lv.txt.gz")))) if !isempty(l)]
    @test length(bad) > 700
    @test all(ek -> !C.kyber_ek_check(ek) && throws(() -> C.kyber_kem_enc(ek)), bad)
end

# CCTV accumulated ML-DSA: seeds from SHAKE128(""), pk and deterministic sig of "" absorbed into SHAKE128.
const ACC_MLDSA_10K = Dict(44 => "e7fd21f6a59bcba60d65adc44404bb29a7c00e5d8d3ec06a732c00a306a7d143",
                           65 => "5ff5e196f0b830c3b10a9eb5358e7c98a3a20136cb677f3ae3b90175c3ace329",
                           87 => "80a8cf39317f7d0be0e24972c51ac152bd2a3e09bc0c32ce29dd82c4e7385e60")
const ACC_N = get(ENV, "PQJULIA_LONG_TESTS", "") == "1" ? 10_000 : 100
const ACC_MLDSA_100 = Dict(44 => "d51148e1f9f4fa1a723a6cf42e25f2a99eb5c1b378b3d2dbbd561b1203beeae4",
                          65 => "8358a1843220194417cadbc2651295cd8fc65125b5a5c1a239a16dc8b57ca199",
                          87 => "8c3ad714777622b8f21ce31bb35f71394f23bc0fcf3c78ace5d608990f3b061b")
# 100 keys by default; PQJULIA_LONG_TESTS=1 runs the 10 000-key CI tier.
@testset "CCTV accumulated ML-DSA-$lv ($ACC_N keys)" for (lv, C) in [(44, MLDSA.Category2), (65, MLDSA.Category3), (87, MLDSA.Category5)]
    s = ShakeReader(UInt8[]; bits=128); a = SHA.SHAKE_128_CTX(); ok = true
    for _ in 1:ACC_N
        pk, sk = C.dilithium_keygen_derand(take!(s, 32))
        sig = C.dilithium_sign_derand(UInt8[], sk, zeros(UInt8, 32))
        SHA.update!(a, pk); SHA.update!(a, sig)
        ok &= C.dilithium_verify(UInt8[], sig, pk)
    end
    d = Vector{UInt8}(undef, 32); SHA.digest!(a, UInt64(32), pointer(d))
    @test ok
    @test bytes2hex(d) == (ACC_N == 100 ? ACC_MLDSA_100 : ACC_MLDSA_10K)[lv]
end
