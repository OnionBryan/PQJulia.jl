# test/falcon_fpr.jl — the integer-emulated signer (falcon_fpr.jl) against hardware binary64:
# each operation bit for bit, SamplerZ on the reference vectors, and whole signatures.
using Test, JSON, Random
using PQJulia
const FP = PQJulia.FNDSA.FalconFpr

@testset "Falcon integer floating point" begin
    rng = Xoshiro(2026)
    rnd() = (k = rand(rng, 1:12);
             k == 1 ? rand(rng, (0.0, -0.0)) :
             k == 2 ? Float64(rand(rng, -100000:100000)) :
             k == 3 ? rand(rng, (-1, 1)) * (1.0 + rand(rng, 0:7) * eps()) :
             rand(rng, (-1.0, 1.0)) * ldexp(1.0 + rand(rng), rand(rng, -60:60)))
    bits(x) = reinterpret(UInt64, x)
    eq(a, b) = bits(a) == bits(b) || (a == 0 && b == 0)
    bad = 0
    for _ in 1:100_000
        a, b = rnd(), rnd(); A, B = FP.fpr(a), FP.fpr(b)
        bad += !eq(Float64(FP.add(A, B)), a + b) + !eq(Float64(FP.sub(A, B)), a - b) +
               !eq(Float64(FP.mul(A, B)), a * b) + (b != 0 && !eq(Float64(FP.div(A, B)), a / b)) +
               (bits(Float64(sqrt(FP.fpr(abs(a))))) != bits(sqrt(abs(a)))) +
               (floor(FP.fpr(a)) != floor(Int, a)) + (FP.rint(FP.fpr(a)) != round(Int, a))
        i = rand(rng, -(2^62):2^62) >> rand(rng, 0:62)
        bad += bits(Float64(FP.fpr_of(i))) != bits(Float64(i))
    end
    @test bad == 0
    @test all(k -> FP.rint(FP.fpr(k + 0.5)) == round(Int, k + 0.5), -1000:1000)   # ties to even

    sz = JSON.parsefile(joinpath(@__DIR__, "kat", "falcon_samplerz_kat.json"))["tests"]
    @test count(sz) do t
        FP.samplerz(FP.fpr(parse(Float64, t["mu"])), FP.fpr(parse(Float64, t["sigma"])),
                    FP.fpr(parse(Float64, t["sigmin"])), FNDSA.FalconSampler.KATSource(t["octets"])) != t["z"]
    end == 0
end

@testset "Falcon integer signer = Float64 signer" begin
    F = FNDSA.Falcon
    for n in (4, 16, 64, 256, 512, 1024), trial in 1:(n >= 512 ? 2 : 4)
        sk = F.keygen(n)
        ek = F.expand_sk(sk); gs = F.sign_setup(sk)
        @test FP.leaves(ek.gs.T) == F.leaves(gs.T)
        msg = rand(UInt8, 33)
        seed = rand(UInt8, 48)
        rb() = (r = Xoshiro(reinterpret(UInt64, seed[1:8])[1]); k -> rand(r, UInt8, k))
        s_int = F.sign_poly(sk, ek.gs, msg; randombytes = rb())
        s_f64 = F.sign_poly(sk, gs, msg; randombytes = rb())
        @test s_int == s_f64
    end
end
