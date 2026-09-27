# X25519 (RFC 7748) and X-Wing (draft-connolly-cfrg-xwing-kem-11).
using Random, JSON, CodecZlib
const XC = PQJulia.X25519
hb(s) = hex2bytes(s)

@testset "X25519 (RFC 7748)" begin
    @testset "§5.2 vectors" begin
        for (k, u, out) in (
            ("a546e36bf0527c9d3b16154b82465edd62144c0ac1fc5a18506a2244ba449ac4",
             "e6db6867583030db3594c1a424b15f7c726624ec26b3353b10a903a6d0ab1c4c",
             "c3da55379de9c6908e94ea4df28d084f32eccf03491c71f754b4075577a28552"),
            ("4b66e9d4d1b4673c5ad22691957d6af5c11b6421e0ea01d42ca4169e7918ba0d",
             "e5210f12786811d3f4b7959d0538ae2c31dbe7106fc03c3efc4cd549c715a493",
             "95cbde9476e8907d7aade45cb4b873f88b595a68799fa152e6f8f7647aac7957"))
            @test bytes2hex(XC.x25519(hb(k), hb(u))) == out
            @test bytes2hex(XC.x25519_ref(hb(k), hb(u))) == out
        end
    end

    @testset "§5.2 iterated" begin
        iter(n) = (k = copy(XC.BASE); u = copy(XC.BASE); for _ in 1:n; k, u = XC.x25519(k, u), k; end; bytes2hex(k))
        @test iter(1) == "422c8e7a6227d7bca1350b3e2bb7279f7897b87bb6854b783c60e80311ae3079"
        @test iter(1000) == "684cf59ba83309552800ef566f2f4d3c1c3887c49360e3875f2eb94d99532c51"
        if get(ENV, "PQJULIA_LONG_TESTS", "") == "1"
            @test iter(1_000_000) == "7c3911e0ab2586fd864497297e575e6f3bc601c0883c30df5f4dd2d24f665424"
        end
    end

    @testset "§6.1 Diffie–Hellman" begin
        a = hb("77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a")
        b = hb("5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb")
        A = XC.x25519_base(a); B = XC.x25519_base(b)
        @test bytes2hex(A) == "8520f0098930a754748b7ddcb43ef75a0dbf3a0d26381af4eba4a98eaa9b4e6a"
        @test bytes2hex(B) == "de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f"
        @test bytes2hex(XC.x25519(a, B)) == bytes2hex(XC.x25519(b, A)) ==
              "4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742"
    end

    @testset "Wycheproof (518 vectors, twist / low-order / non-canonical inputs)" begin
        g = JSON.parse(read(GzipDecompressorStream(open(joinpath(@__DIR__, "kat", "wycheproof", "x25519_test.json.gz"))), String))["testGroups"]
        tests = reduce(vcat, [t["tests"] for t in g])
        @test length(tests) == 518
        @test all(XC.x25519(hb(t["private"]), hb(t["public"])) == hb(t["shared"]) for t in tests)
    end

    @testset "independent Edwards-model oracle (paper5-haskell)" begin
        vs = JSON.parsefile(joinpath(@__DIR__, "kat", "x25519_edwards_oracle.json"))["tests"]
        @test length(vs) == 164
        @test all(bytes2hex(XC.x25519(hb(v["scalar"]), hb(v["u"]))) == v["out"] for v in vs)
    end

    @testset "limb engine == BigInt reference on random inputs" begin
        rng = MersenneTwister(7748)
        @test all(1:300) do _
            k = rand(rng, UInt8, 32); u = rand(rng, UInt8, 32)
            XC.x25519(k, u) == XC.x25519_ref(k, u)
        end
        # u ∈ [p, 2²⁵⁵) is reduced, and the top bit is ignored.
        p = big(2)^255 - 19
        for x in (p, p + 1, p + 18, big(2)^255 - 1)
            u = XC.tobytes(x); k = rand(rng, UInt8, 32)
            @test XC.x25519(k, u) == XC.x25519_ref(k, u) == XC.x25519(k, XC.tobytes(mod(x, p)))
        end
        k = rand(rng, UInt8, 32); u = rand(rng, UInt8, 32); u2 = copy(u); u2[32] ⊻= 0x80
        @test XC.x25519(k, u) == XC.x25519(k, u2)
        @test_throws ArgumentError XC.x25519(k[1:31], u)
    end
end

@testset "X-Wing (draft-connolly-cfrg-xwing-kem-11)" begin
    @testset "Appendix C vectors" begin
        vs = JSON.parsefile(joinpath(@__DIR__, "kat", "xwing_draft11.json"))["tests"]
        @test length(vs) == 3
        for v in vs
            pk, sk = XWing.xwing_keypair_derand(hb(v["seed"]))
            @test bytes2hex(sk) == v["sk"] && bytes2hex(pk) == v["pk"]
            ct, ss = XWing.xwing_encaps_derand(pk, hb(v["eseed"]))
            @test bytes2hex(ct) == v["ct"] && bytes2hex(ss) == v["ss"]
            @test bytes2hex(XWing.xwing_decaps(hb(v["ct"]), hb(v["sk"]))) == v["ss"]
        end
    end

    @testset "round trip and rejection" begin
        pk, sk = XWing.xwing_keypair()
        @test (length(pk), length(sk)) == (XWing.PK_BYTES, XWing.SK_BYTES)
        ct, ss = XWing.xwing_encaps(pk)
        @test length(ct) == XWing.CT_BYTES && XWing.xwing_decaps(ct, sk) == ss
        badM = copy(ct); badM[10] ⊻= 0x01                  # ML-KEM part: implicit rejection
        badX = copy(ct); badX[1100] ⊻= 0x01                # X25519 part: different ss_X
        @test XWing.xwing_decaps(badM, sk) != ss && XWing.xwing_decaps(badX, sk) != ss
        @test XWing.xwing_decaps(ct, XWing.xwing_keypair()[2]) != ss
        @test_throws ArgumentError XWing.xwing_encaps(pk[1:end-1])
        @test_throws ArgumentError XWing.xwing_decaps(ct[1:end-1], sk)
        badpk = copy(pk); badpk[1:2] .= 0xff                # ML-KEM modulus check (FIPS 203 §7.2)
        @test_throws ArgumentError XWing.xwing_encaps(badpk)
    end
end
