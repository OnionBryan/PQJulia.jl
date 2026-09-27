# FALCON-MRM, message recovery (ePrint 2026/420).
using Random
const MRM = FNDSA.FalconMRM

@testset "FALCON-MRM" begin
    @testset "sizes (ePrint 2026/420 Tables 1–2)" begin
        @test (MRM.m1_len(512), MRM.m1_len(1024)) == (478, 962)
        @test (MRM.sig_bytes(512), MRM.sig_bytes(1024)) == (1251, 2479)
        @test (MRM.max_bits(478), MRM.max_bits(962)) == (6492, 13067)
        @test FNDSA.Falcon512.MRM_MAX_BITS == floor(Int, 478 * log2(12289)) - 1
    end

    @testset "§1.4 encode/decode" begin
        rng = MersenneTwister(420)
        for k in (478, 962, 12, 2, 14)
            μ = rand(rng, Bool, MRM.max_bits(k))
            z = MRM.encode(μ, k)
            @test length(z) == k && all(x -> 0 <= x < 12289, z) && MRM.decode(z, length(μ)) == μ
        end
        @test MRM.encode(rand(rng, Bool, 6491), 478) === nothing
        @test MRM.decode(fill(12288, 12), 163) === nothing          # (q¹² − 1) ≥ 2¹⁶³
        @test MRM.decode(fill(12288, 2), 27) === nothing            # q² − 1 ≥ 2²⁷
        @test MRM.decode([12289; zeros(Int, 11)], 163) === nothing
    end

    @testset "sign / recover ($(M.IDENTIFIER))" for M in (FNDSA.Falcon512, FNDSA.Falcon1024)
        rng = MersenneTwister(2026)
        pk, sk = M.falcon_keygen()
        ek = M.falcon_expand_sk(sk)
        μ = rand(rng, Bool, M.MRM_MAX_BITS)
        M1 = MRM.encode(μ, M.MRM_M1_LEN); M2 = Vector{UInt8}("header")
        sig = M.falcon_mrm_sign(M1, M2, ek)
        @test length(sig) == M.MRM_SIG_BYTES
        rec = M.falcon_mrm_verify(M2, sig, pk)
        @test rec == M1 && MRM.decode(rec, length(μ)) == μ
        @test M.falcon_mrm_verify(UInt8[M2; 0x00], sig, pk) === nothing
        bad = copy(sig); bad[100] ⊻= 0x04
        @test M.falcon_mrm_verify(M2, bad, pk) === nothing
        @test M.falcon_mrm_verify(M2, sig[1:end-1], pk) === nothing
        @test M.falcon_mrm_verify(M2, sig, M.falcon_keygen()[1]) === nothing
        @test M.falcon_mrm_sign(M1, M2, ek) != sig                    # fresh ρ per signature
        @test_throws ArgumentError M.falcon_mrm_sign(M1[1:end-1], M2, ek)
        @test_throws ArgumentError M.falcon_mrm_sign([M1[1:end-1]; 12289], M2, ek)
        # A plain Falcon signature is not an MRM signature and vice versa.
        @test !M.falcon_verify(M2, sig, pk)
        @test M.falcon_mrm_verify(M2, M.falcon_sign(M2, ek), pk) === nothing
    end
end
