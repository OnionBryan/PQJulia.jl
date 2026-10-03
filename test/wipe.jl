# Secret-buffer wiping: wipe! zeroes and survives compilation; no API call touches caller inputs.
using InteractiveUtils
const W = PQJulia.Wipe

@testset "wipe!" begin
    for T in (UInt8, Int16, Int32, Int64, Float64, ComplexF64)
        x = rand(T, 37) .+ one(T)
        @test all(iszero, W.wipe!(x))
    end
    v = [rand(Int32, 8) for _ in 1:3]; a = rand(UInt8, 5)
    W.wipe!(v, a)
    @test all(p -> all(iszero, p), v) && all(iszero, a)
    ir = sprint(io -> code_llvm(io, W.wipe!, (Vector{UInt8},); debuginfo=:none))
    @test occursin("asm sideeffect", ir)                      # the zeroing is not a dead store
end

@testset "API leaves caller inputs intact" begin
    unchanged(f, xs...) = (cs = deepcopy(xs); f(xs...); all(map(==, cs, xs)))
    for C in (MLKEM.Category1, MLKEM.Category3, MLKEM.Category5)
        coins = rand(UInt8, 64); m = rand(UInt8, 32)
        pk, sk = C.kyber_kem_keypair_derand(copy(coins))
        @test unchanged(C.kyber_kem_keypair_derand, coins)
        @test unchanged(C.kyber_kem_enc_derand, pk, m)
        ct, ss = C.kyber_kem_enc_derand(pk, m)
        @test unchanged(C.kyber_kem_dec, ct, sk)
        @test C.kyber_kem_dec(ct, sk) == ss == C.kyber_kem_enc_derand(pk, m)[2]
    end
    for C in (MLDSA.Category2, MLDSA.Category3, MLDSA.Category5)
        xi = rand(UInt8, 32); rnd = rand(UInt8, 32); msg = rand(UInt8, 50)
        @test unchanged(C.dilithium_keygen_derand, xi)
        pk, sk = C.dilithium_keygen_derand(xi)
        @test unchanged((s, r) -> C.dilithium_sign_derand(msg, s, r), sk, rnd)
        @test C.dilithium_sign_derand(msg, sk, rnd) == C.dilithium_sign_derand(msg, sk, rnd)
        @test C.dilithium_verify(msg, C.dilithium_sign(msg, sk), pk)
    end
    seed = rand(UInt8, 32); eseed = rand(UInt8, 64)
    @test unchanged(XWing.xwing_keypair_derand, seed)
    pk, sk = XWing.xwing_keypair_derand(seed)
    @test unchanged(e -> XWing.xwing_encaps_derand(pk, e), eseed)
    ct, ss = XWing.xwing_encaps_derand(pk, eseed)
    @test unchanged(s -> XWing.xwing_decaps(ct, s), sk)
    @test XWing.xwing_decaps(ct, sk) == ss
end

@testset "Falcon key wiping" begin
    F = FNDSA.Falcon512
    pk, sk = F.falcon_keygen(); msg = rand(UInt8, 40)
    s0 = copy(sk)
    @test F.falcon_verify(msg, F.falcon_sign(msg, sk), pk) && sk == s0
    @test F.falcon_verify(msg, F.falcon_sign(msg, sk), pk)   # the byte key still signs
    ek = F.falcon_expand_sk(sk)
    @test F.falcon_verify(msg, F.falcon_sign(msg, ek), pk)
    F.falcon_wipe!(ek)
    @test all(iszero, ek.sk.f) && all(iszero, ek.sk.F) && all(iszero, ek.gs.b00)
    @test all(iszero, FNDSA.Falcon.leaves(ek.gs.T))
end
