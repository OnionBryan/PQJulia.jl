# Exact (float-free) key certificate: FNDSA.FalconCertify.
using Random

const FCe = FNDSA.FalconCertify
const Fal = FNDSA.Falcon

# NTRU-complete (f, g, F, G) from f, g, or nothing if f is not invertible mod q.
function ntru_key(f, g)
    Fal.is_invertible(f) || return nothing
    F, G = try FNDSA.FalconNTRUGen.ntru_solve(BigInt.(f), BigInt.(g)) catch; return nothing end
    Fal.reduce_FG!(BigInt.(f), BigInt.(g), F, G)
    Fal.secret_key(f, g, F, G)
end

@testset "Falcon exact key certificate" begin
    @testset "multi-modular leaves == ℚ(x) tower, n=$n" for n in (2, 4, 8, 16)
        sk = Fal.keygen(n)
        d = FCe.exact_leaves(sk.f, sk.g, sk.F, sk.G)
        @test d == FCe.exact_leaves_tower(sk.f, sk.g, sk.F, sk.G)
        @test prod(d) == big(12289)^n
    end

    @testset "float ffLDL leaves agree with exact, n=$n" for n in (64, 512)
        sk = Fal.keygen(n)
        c = Fal.certify(sk)
        σ = Fal.params(n).σ
        exact = [σ / sqrt(Float64(x)) for x in c.leaves]
        @test maximum(abs.(Fal.leaves(Fal.sign_setup(sk).T) .- exact) ./ exact) < 1e-12
    end

    @testset "KAT keys certify" begin
        kat = JSON.parsefile(joinpath(@__DIR__, "kat", "falcon_sign_kat.json"))
        for t in kat["tests"][1:12:end]                          # one key per n = 2..1024
            @test Fal.certify((Int.(t[k]) for k in ("f", "g", "F", "G"))...).ok
        end
    end

    @testset "public API" begin
        pk, sk = FNDSA.Falcon512.falcon_keygen(certified=true)
        c = FNDSA.Falcon512.falcon_certify(sk)
        @test c.ok && c.leaves_ok && c.gs_ok && length(c.leaves) == 512
        @test FNDSA.Falcon512.falcon_verify(UInt8[1, 2, 3], FNDSA.Falcon512.falcon_sign(UInt8[1, 2, 3], sk), pk)
        @test_throws ArgumentError FNDSA.Falcon512.falcon_certify(sk[1:end-1])
    end

    @testset "FIPS 206: GS bound 0.9999·1.17√q, NTT-form public key" begin
        sk = Fal.keygen(64; fips206=true)
        c = Fal.certify(sk; fips206=true)
        @test c.ok && c.gs_norm2 <= (9999 // 10000 * 117 // 100)^2 * 12289
        # The GS decision is exact: a factor 2⁻²⁰⁰ either side of √(gs²/q) flips it.
        r = c.gs_norm2 / 12289
        lo = isqrt(numerator(r) * big(2)^400 ÷ denominator(r)) // big(2)^200
        p = Fal.params(64)
        cert(k) = FCe.certify(sk.f, sk.g, sk.F, sk.G; σ=p.σ, σmin=p.σmin, σmax=FNDSA.FalconSampler.MAX_SIGMA, gs_factor=k)
        @test !cert(lo).gs_ok && cert(lo + 1 // big(2)^200).gs_ok
        ĥ = Fal.pk_ntt(sk.h)
        @test mod.(ĥ .* Fal.ntt_fwd(mod.(sk.f, 12289)), 12289) == Fal.ntt_fwd(mod.(sk.g, 12289))
        ek = Fal.expand_sk(Fal.encode_sk(sk))
        for i in 1:5
            m = UInt8[i]; sig = Fal.sign_poly(ek.sk, ek.gs, m)
            @test Fal.verify_poly_ntt(ĥ, 64, m, sig) && Fal.verify_poly(sk.h, 64, m, sig)
            @test !Fal.verify_poly_ntt(ĥ, 64, UInt8[i, 0], sig)
        end
    end

    @testset "oversized basis is rejected, exactly and by the signer" begin
        rng = MersenneTwister(2026)
        sk = nothing
        while sk === nothing
            sk = ntru_key(rand(rng, -40:40, 16), rand(rng, -40:40, 16))
        end
        c = Fal.certify(sk)
        @test !c.ok && !c.gs_ok && !c.leaves_ok
        @test prod(c.leaves) == big(12289)^16
        @test minimum(c.leaves) < c.bounds[1] || maximum(c.leaves) > c.bounds[2]
        @test_throws ArgumentError Fal.sign_setup(sk)
    end
end
