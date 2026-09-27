# test/falcon_kat.jl — Falcon against the round-3 C implementation (test/kat/falcon_*_kat.json,
# vendored from tprest/falcon.py, MIT). Signing KAT: from the same (f,g,F,G), message and
# SHAKE256("external") randomness stream, the padded signature must match byte-for-byte.
using Test, JSON, SHA
using PQJulia
const FAL = PQJulia.FNDSA.Falcon
const FS = PQJulia.FNDSA.FalconSampler

# One SHAKE256 stream is shared across all vectors, as in the reference.
include("shake_reader.jl")

@testset "Falcon KATs (round-3 C reference)" begin
    kat_dir = joinpath(@__DIR__, "kat")

    sz = JSON.parsefile(joinpath(kat_dir, "falcon_samplerz_kat.json"))["tests"]
    @testset "SamplerZ bit-exact ($(length(sz)) vectors)" begin
        bad = count(sz) do t
            FS.samplerz(parse(Float64, t["mu"]), parse(Float64, t["sigma"]), parse(Float64, t["sigmin"]),
                        FS.KATSource(t["octets"])) != t["z"]
        end
        @test bad == 0
    end

    kat = JSON.parsefile(joinpath(kat_dir, "falcon_sign_kat.json"))
    msg = Vector{UInt8}(kat["message"]); vecs = kat["tests"]
    seed = Vector{UInt8}(kat["shake_seed"])
    @test take!(ShakeReader(seed), 1000) == SHA.shake256(seed, UInt64(1000))
    stream = ShakeReader(seed)

    for n in sort(unique(t["n"] for t in vecs))
        vn = [t for t in vecs if t["n"] == n]
        signed = 0; verified = 0; tamper = 0; skround = 0
        for t in vn
            sk = FAL.secret_key(t["f"], t["g"], t["F"], t["G"])
            ref = hex2bytes(t["sig"]); pk = FAL.encode_pk(sk.h, n)
            skround += FAL.decode_sk(FAL.encode_sk(sk)) == sk
            take!(stream, 8 * t["read_bytes"])
            sig = FAL.sign(msg, FAL.expand_sk(sk); randombytes = k -> take!(stream, k))
            signed += sig == ref
            verified += FAL.verify(msg, ref, pk)
            bad = copy(ref); bad[end-3] ⊻= 0x10
            tamper += !FAL.verify(msg, bad, pk) && !FAL.verify(UInt8[msg; 0x00], ref, pk)
        end
        @testset "n=$n ($(length(vn)) vectors)" begin
            @test signed == length(vn)
            @test verified == length(vn)
            @test tamper == length(vn)
            @test skround == length(vn)
        end
    end
end
