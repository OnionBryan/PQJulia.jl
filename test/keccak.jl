# Keccak (src/keccak.jl) against SHA.jl: every input length through three blocks of each rate,
# outputs at and across the rate boundaries, non-Vector inputs, and the FIPS 202 empty-string values.
using SHA

@testset "Keccak vs SHA.jl" begin
    K = PQJulia.Keccak
    @test bytes2hex(K.sha3_256(UInt8[])) == "a7ffc6f8bf1ed76651c14756a061d662f580ff4de43b49fa82d80a4b80f8434a"
    @test bytes2hex(K.shake128(UInt8[], 16)) == "7f9c2ba4e88f827d616045507605853e"
    bad = 0
    for len in 0:3 * 168 + 1
        d = rand(UInt8, len)
        K.sha3_256(d) == SHA.sha3_256(d) || (bad += 1)
        K.sha3_512(d) == SHA.sha3_512(d) || (bad += 1)
        K.sha3_224(d) == SHA.sha3_224(d) || (bad += 1)
        K.sha3_384(d) == SHA.sha3_384(d) || (bad += 1)
        for out in (0, 1, 32, 64, 135, 136, 137, 167, 168, 169, 672, 840)
            K.shake128(d, out) == SHA.shake128(d, UInt64(out)) || (bad += 1)
            K.shake256(d, out) == SHA.shake256(d, UInt64(out)) || (bad += 1)
        end
        v = @view vcat(UInt8[0xaa], d)[2:end]
        K.shake256(v, 100) == SHA.shake256(d, UInt64(100)) || (bad += 1)
    end
    @test bad == 0
    d = rand(UInt8, 33)
    @test K.shake256(d, 10_000) == SHA.shake256(d, UInt64(10_000))
end
