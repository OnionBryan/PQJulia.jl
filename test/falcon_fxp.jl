# Float-free key generation (FNDSA.FalconFxp, after ePrint 2023/290 and pornin/ntrugen).
using Random, SHA
const FX = FNDSA.FalconFxp

@testset "Falcon fixed-point keygen" begin
    @testset "root table equals ntrugen GM_TAB" begin
        # SHA-256 of ntrugen ng_fxp.c GM_TAB[0..1023] as little-endian (re, im) u64 pairs.
        b = reinterpret(UInt8, reduce(vcat, [[reinterpret(UInt64, r), reinterpret(UInt64, i)] for (r, i) in FX.GM_TAB]))
        @test bytes2hex(sha256(collect(b))) == "e0746da83816f05bd5a3bbe0867bba11ee2bab9cd781780ce4d205d6c51b03d9"
    end

    @testset "fxr division is round-to-nearest of x·2³²/y" begin
        rng = MersenneTwister(290)
        for _ in 1:2000
            x = rand(rng, -(Int64(1) << 50):(Int64(1) << 50)); y = rand(rng, -(Int64(1) << 45):(Int64(1) << 45))
            y == 0 && continue
            r = big(x) << 32 // y
            @test FX.fxr_div(x, y) == sign(r) * floor(BigInt, abs(r) + 1 // 2)
        end
    end

    @testset "fixed-point FFT product is the negacyclic product, n=$n" for n in (4, 16, 256, 1024)
        rng = MersenneTwister(n)
        a = rand(rng, -40:40, n); b = rand(rng, -40:40, n)
        ra = FX.vect_FFT!(FX.fxr_of.(a)); rb = FX.vect_FFT!(FX.fxr_of.(b))
        @test [FX.fxr_round(x) for x in FX.vect_iFFT!(FX.vect_mul_fft!(ra, rb))] == FNDSA.FalconNTRUGen.negamul(a, b)
    end

    @testset "fixed-point GS check agrees with the exact GS norm" begin
        rng = FNDSA.FalconChaCha.ChaCha20(collect(UInt8, 1:56))
        tested = 0; agree = 0; near = 0; passed = 0
        for _ in 1:300
            f = FX.gauss_sample_poly(64, rng); g = FX.gauss_sample_poly(64, rng)
            tested += 1
            r = FNDSA.FalconCertify.exact_gs_norm2(f, g) / ((117 // 100)^2 * 12289)
            passed += r < 1
            if abs(r - 1) < 1 // 10^6
                near += 1
            else
                agree += FX.gs_ok(f, g) == (r < 1)
            end
        end
        @test agree == tested - near && near <= 1 && passed > 20
    end

    @testset "keys (n=$n)" for n in (4, 64, 512)
        sk = Fal.keygen(n; fixedpoint=true)
        B(v) = BigInt.(v)
        @test FNDSA.FalconNTRUGen.ntru_check(B(sk.f), B(sk.g), B(sk.F), B(sk.G))
        @test all(x -> -128 < x < 128, [sk.F; sk.G])
        @test isodd(sum(sk.f)) && isodd(sum(sk.g))
        @test Fal.certify(sk).ok
        ek = Fal.expand_sk(sk)
        @test all(Fal.verify(UInt8[i], Fal.sign(UInt8[i], ek), Fal.encode_pk(sk.h, n)) for i in 1:3)
    end

    @testset "public API, with the FIPS 206 bound" begin
        pk, sk = FNDSA.Falcon512.falcon_keygen(fixedpoint=true, fips206=true)
        @test FNDSA.Falcon512.falcon_certify(sk; fips206=true).ok
        @test FNDSA.Falcon512.falcon_verify(UInt8[7], FNDSA.Falcon512.falcon_sign(UInt8[7], sk), pk)
        @test_throws ArgumentError Fal.keygen(2; fixedpoint=true)
    end
end
