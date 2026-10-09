# test/falcon_hardening.jl — checks on the integer signer's ffLDL leaves (falcon_fpr.jl): the range
# [σmin, σmax] is re-checked before every signature, and each leaf σ/√d is evaluated twice at key
# expansion to detect a fault (Kaihara et al., ePrint 2026/2046). Signatures are unchanged.
using Test, JSON, SHA, InteractiveUtils
using PQJulia

@testset "Falcon signer hardening" begin
    FP = PQJulia.FNDSA.FalconFpr; FAL = PQJulia.FNDSA.Falcon; MRM = PQJulia.FNDSA.FalconMRM
    σmax = PQJulia.FNDSA.FalconSampler.MAX_SIGMA
    leaflist(T, out=FP.Leaf[]) = T isa FP.Leaf ? push!(out, T) : (leaflist(T.t0, out); leaflist(T.t1, out))
    # Fixed randomness: one SHAKE256 output read in order.
    stream(tag) = (buf = SHA.shake256(Vector{UInt8}(tag), UInt64(8192)); pos = Ref(0);
                   k -> (v = buf[pos[]+1:pos[]+k]; pos[] += k; v))
    kat = JSON.parsefile(joinpath(@__DIR__, "kat", "falcon_sign_kat.json"))["tests"]
    katkey(n) = (t = first(x for x in kat if x["n"] == n); FAL.secret_key(t["f"], t["g"], t["F"], t["G"]))
    msg = Vector{UInt8}("PQJulia hardening")
    m1(n) = [mod(37i, 12289) for i in 1:MRM.m1_len(n)]
    # SHA-256 of the Falcon and FALCON-MRM signatures produced before these checks were added.
    pinned = Dict(512  => ("2427199e3cf5e18e3f477624b429a11ba518f810d6a1f2be0baae0a9ab2503bd",
                           "242b26ce3796ff06c1978bcf442650fb889b6c1c79106e3199cd69cdfabb3cf4"),
                  1024 => ("b7a17ccc8fa58fc12d8f9ab18ab3249164f259555de8b8bd526fea0bbd00c56b",
                           "310cd575e8ef152eb81876423dc032d3ec1f237dc336eb6e1a1bc942db5000c2"))

    @testset "leaf range test at the bounds" begin
        σmin = FAL.params(512).σmin
        bad(v) = FP.leafbad(FP.Leaf(FP.fpr(v)), FP.fpr(σmin), FP.fpr(σmax))
        @test bad(σmin) == bad(σmax) == bad((σmin + σmax) / 2) == 0
        @test all(v -> bad(v) == 1, (prevfloat(σmin), nextfloat(σmax), 0.0, -0.0, -1.5, Inf, -Inf, NaN, -NaN))
    end

    @testset "untouched key signs as before (n=$n)" for n in (512, 1024)
        sk = katkey(n); pk = FAL.encode_pk(sk.h, n)
        ek = FAL.expand_sk(sk)
        sig = FAL.sign(msg, ek; randombytes = stream("sign $n"))
        @test bytes2hex(sha256(sig)) == pinned[n][1]
        @test FAL.verify(msg, sig, pk)
        # The Float64 reference signer (no re-check) gives the same bytes.
        @test FAL.encode_sig(FAL.sign_poly(sk, FAL.sign_setup(sk), msg; randombytes = stream("sign $n")), n) == sig
        ms = MRM.sign(m1(n), msg, ek; randombytes = stream("mrm $n"))
        @test bytes2hex(sha256(ms)) == pinned[n][2]
        @test MRM.verify(msg, ms, sk.h) == m1(n)
    end

    @testset "corrupted expanded key is refused ($(M.IDENTIFIER))" for M in (FNDSA.Falcon512, FNDSA.Falcon1024)
        n = M.N; sk = katkey(n); pk = FAL.encode_pk(sk.h, n); σmin = FAL.params(n).σmin
        ek = FAL.expand_sk(sk)
        L = leaflist(ek.gs.T)
        @test length(L) == n
        ref = FAL.sign(msg, ek; randombytes = stream("ref"))
        draws = Ref(0); counted = k -> (draws[] += 1; rand(UInt8, k))
        for (i, v) in ((1, prevfloat(σmin)), (n ÷ 2, nextfloat(σmax)), (n, 0.0), (3, NaN), (n - 5, -1.5))
            old = L[i].σ; L[i].σ = FP.fpr(v)
            @test_throws ArgumentError M.falcon_sign(msg, ek)
            @test_throws ArgumentError M.falcon_mrm_sign(m1(n), msg, ek)
            @test_throws ArgumentError FAL.sign(msg, ek; randombytes = counted)
            @test_throws ArgumentError MRM.sign(m1(n), msg, ek; randombytes = counted)
            L[i].σ = old
        end
        @test draws[] == 0                                  # refused before any randomness is drawn
        @test FAL.sign(msg, ek; randombytes = stream("ref")) == ref
        @test M.falcon_verify(msg, M.falcon_sign(msg, ek), pk)
        @test M.falcon_mrm_verify(msg, M.falcon_mrm_sign(m1(n), msg, ek), pk) == m1(n)
        # A wiped key has zero leaves; signing with it used to resample forever. A bounded randomness
        # source turns a missing check into a test failure instead of a hang.
        M.falcon_wipe!(ek)
        budget(k) = (c = Ref(0); m -> (c[] += 1; c[] > k && error("randomness budget exhausted"); rand(UInt8, m)))
        @test_throws ArgumentError FAL.sign(msg, ek; randombytes = budget(10_000))
        @test_throws ArgumentError MRM.sign(m1(n), msg, ek; randombytes = budget(10_000))
    end

    @testset "duplicate leaf evaluation" begin
        σ = FP.fpr(FAL.params(512).σ)
        for d in FP.fpr.([8300.0, 9329.071102002914, 12289.0, 16900.0])
            v, e = FP.leafsigma(σ, d)
            @test e == 0 && Float64(v) === Float64(σ) / sqrt(Float64(d))
            # A fault in either evaluation, modelled as a flipped bit in its operand (bits 1–62; the
            # sign bit does not enter sqrt), changes that evaluation's result and is reported.
            for k in 1:62
                dk = FP.Fpr(d.b ⊻ (UInt64(1) << k))
                v2, e2 = FP.leafsigma(σ, d, dk)
                v1, e1 = FP.leafsigma(σ, dk, d)
                @test v2 == v && e2 == v.b ⊻ FP.div(σ, sqrt(dk)).b != 0
                @test e1 == v1.b ⊻ v.b != 0
            end
        end
        # A mismatch at one leaf, injected through setup's evaluator, aborts key expansion after
        # every leaf has been computed; the unperturbed evaluator gives the expand_sk tree.
        sk = katkey(512); p = FAL.params(512)
        for at in (1, 300, 512)
            calls = Ref(0)
            glitch(s, d) = (calls[] += 1; FP.leafsigma(s, d, calls[] == at ? FP.Fpr(d.b ⊻ 0x10) : d))
            @test_throws ErrorException FP.setup(sk, p.σ, p.σmin; leaf = glitch)
            @test calls[] == 512
        end
        @test leaflist(FP.setup(sk, p.σ, p.σmin; leaf = FP.leafsigma).T) .|> (l -> l.σ) ==
              leaflist(FAL.expand_sk(sk).gs.T) .|> (l -> l.σ)
        # Both evaluations survive compilation: ffldl passes the second operand through the asm
        # barrier, and leafsigma keeps two square roots (inlined loops, r = 2^53 on entry, or calls).
        ir = sprint(io -> code_llvm(io, FP.ffldl, (Vector{FP.CF}, Vector{FP.CF}, Vector{FP.CF}, FP.Fpr,
                                                   Base.RefValue{UInt64}, typeof(FP.leafsigma)); debuginfo=:none))
        @test count("asm sideeffect \"\", \"=r,0\"", ir) == 2
        ir = sprint(io -> code_llvm(io, FP.leafsigma, (FP.Fpr, FP.Fpr, FP.Fpr); debuginfo=:none))
        @test count(r"phi i64 [^\n]*\[ 9007199254740992, ", ir) + count(r"call .*@j(ulia)?_sqrt_\d+\(", ir) == 2
    end

    if Sys.ARCH === :x86_64
        @testset "leaf range test has no conditional branch (x86-64)" begin
            asm = sprint(io -> code_native(io, FP.leafbad, (FP.Leaf, FP.Fpr, FP.Fpr); debuginfo=:none, syntax=:intel))
            @test !occursin(r"^\s+j(?!mp\b)[a-z]+\s"m, asm)
        end
    end
end
