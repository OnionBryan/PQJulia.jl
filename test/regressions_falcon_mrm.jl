# test/regressions_falcon_mrm.jl — regressions from the defect sweep for FALCON-MRM signing and the
# byte-key Falcon API (src/falcon/falcon_mrm.jl, src/fndsa.jl). Each testset failed at 2e9a475.
#   F13, F23, F31  MRM signing left the sampler seed, the ChaCha20 state and keystream, the targets
#                  t0, t1 (which give f and F from the public c) and rejected candidates in memory
#   F10(a)         a byte key's ExpandedKey was not wiped when signing threw; M1 is now checked
#                  before the key is decoded
#   F36            encode(μ, k) returned k − 1 elements for odd k
#   inputs         offset-axes vectors, negative k, decode with a bit count encode cannot produce
# Watching package-internal buffers needs wrappers around package methods; those run in a child
# process, so this process keeps the package's own methods.
using Test, PQJulia

@testset "Falcon MRM / byte API regressions (sweep)" begin
    F = FNDSA.Falcon512; MRM = FNDSA.FalconMRM; FAL = FNDSA.Falcon
    pk, sk = F.falcon_keygen(); ek = F.falcon_expand_sk(sk)
    M1 = [mod(37i, 12289) for i in 1:F.MRM_M1_LEN]; M2 = Vector{UInt8}("PQJulia regression")
    bits(B, s) = [isodd((i * 2654435761 + s) >> 7) for i in 1:B]
    # The same values on axes 4:3+length(v), as an OffsetArray would have.
    off(v) = view([zeros(eltype(v), 3); v], Base.IdentityUnitRange(4:3+length(v)))

    @testset "MRM sign wipes its sampler seeds (F13, F23, F31)" begin
        for _ in 1:3
            bufs = Vector{UInt8}[]
            rb = k -> (b = rand(UInt8, k); push!(bufs, b); b)      # keeps every buffer it hands out
            sig = MRM.sign(M1, M2, ek; randombytes = rb)
            @test F.falcon_mrm_verify(M2, sig, pk) == M1
            seeds = filter(b -> length(b) == FAL.SEED_LEN, bufs)
            @test !isempty(seeds) && all(s -> all(iszero, s), seeds)
        end
        # A seed the generator refuses (55 bytes) is wiped before the error propagates.
        short = UInt8[]
        rb = k -> k == FAL.SEED_LEN ? append!(short, rand(UInt8, k - 1)) : rand(UInt8, k)
        @test_throws ArgumentError MRM.sign(M1, M2, ek; randombytes = rb)
        @test length(short) == FAL.SEED_LEN - 1 && all(iszero, short)
    end

    @testset "encode / decode on odd k (F36)" begin
        for k in (1, 3, 13, 25, 477, 479, 961, 963)
            μ = bits(MRM.max_bits(k), k)
            z = MRM.encode(μ, k)
            @test length(z) == k && z[end] == 0 && all(x -> 0 <= x < 12289, z)
            @test z[1:end-1] == MRM.encode(μ, k - 1)              # the odd element carries no bits
            @test MRM.decode(z, length(μ)) == μ
            z[end] = 1
            @test MRM.decode(z, length(μ)) === nothing              # and must be 0
        end
        @test MRM.encode(Bool[], 1) == [0] && MRM.decode([0], 0) == Bool[]
        # Even k is unchanged: SHA3-256 of the encodings produced at 2e9a475.
        for (k, h) in (478 => "439db5d69634dd387994a212d3da9ba0394ba38dfad1114f903a5a0be6484740",
                       962 => "a5b03e4f8edb56007db99fb6f16012d14ad2217ee1cc899842c0cad6f55c4869",
                       14  => "2edd1f5bb9fb66f6e3099be6d5d419c80d1ad87fa3ee24c38969a6992dfdefde")
            z = MRM.encode(bits(MRM.max_bits(k), k), k)
            @test bytes2hex(PQJulia.Keccak.sha3_256(Vector{UInt8}(join(z, ",")))) == h
        end
    end

    @testset "encode / decode accept only what encode produces" begin
        μ = bits(F.MRM_MAX_BITS, 0); z = MRM.encode(μ, F.MRM_M1_LEN)
        @test MRM.decode(z, 163) === nothing                        # was: the first 163 bits
        @test MRM.decode([z; 0; 0], length(μ)) === nothing          # was: extra elements ignored
        @test MRM.decode([z; 0], length(μ)) == μ                    # ℤ_q^479 holds the same bits
        @test_throws ArgumentError MRM.encode(Bool[], -1)           # was: an empty vector
        @test_throws ArgumentError MRM.max_bits(-12)                # was: -163
    end

    @testset "offset-axes inputs" begin
        sig = F.falcon_mrm_sign(off(M1), off(M2), ek)
        @test F.falcon_mrm_verify(M2, sig, pk) == M1
        @test F.falcon_mrm_verify(off(M2), off(sig), off(pk)) == M1
        @test MRM.verify(M2, off(sig), off(FAL.decode_pk(pk, F.N))) == M1
        @test F.falcon_mrm_verify(M2, F.falcon_mrm_sign(M1, M2, off(sk)), pk) == M1
        s = F.falcon_sign(off(M2), off(sk))
        @test F.falcon_verify(M2, s, pk) && F.falcon_verify(off(M2), off(s), off(pk))
        @test F.falcon_verify(M2, F.falcon_sign(off(M2), ek), pk)
        e2 = F.falcon_expand_sk(off(sk))
        @test F.falcon_verify(M2, F.falcon_sign(M2, e2), pk)
        F.falcon_wipe!(e2)
        @test_throws ArgumentError F.falcon_certify(off(zeros(UInt8, F.SK_BYTES)))   # was: BoundsError
        μ = bits(F.MRM_MAX_BITS, 1)
        @test MRM.decode(off(MRM.encode(off(μ), F.MRM_M1_LEN)), length(μ)) == μ
    end
    F.falcon_wipe!(ek)

    # Child process: wrappers that record (and, for the forced rejections, alter) package internals.
    child = raw"""
    redirect_stderr(stdout)
    using PQJulia
    const FAL = FNDSA.Falcon; const F = FNDSA.Falcon512; const MRM = FNDSA.FalconMRM
    const FP = FNDSA.FalconFpr; const CC = FNDSA.FalconChaCha; const FE = FNDSA.FalconEncoding
    check(name, ok) = println("CHECK ", name, " ", ok === true)
    threw(f) = try f(); false catch; true end
    zerocf(v) = all(x -> x.re.b == 0 && x.im.b == 0, v)
    wiped(ek) = all(iszero, ek.sk.f) && all(iszero, ek.sk.g) && all(iszero, ek.sk.F) && all(iszero, ek.sk.G) &&
                zerocf(ek.gs.b00) && zerocf(ek.gs.b01) && zerocf(ek.gs.b10) && zerocf(ek.gs.b11) &&
                all(iszero, FP.leaves(ek.gs.T))
    struct Unreadable <: AbstractVector{UInt8} end            # signing fails when it reads the message
    Base.size(::Unreadable) = (3,)
    Base.getindex(::Unreadable, ::Int) = error("unreadable message")
    pk, sk = F.falcon_keygen()
    M1 = [mod(37i, 12289) for i in 1:F.MRM_M1_LEN]; M2 = UInt8[1, 2, 3]

    # F10(a): record every key the byte-key API decodes and the ExpandedKey it builds.
    const EKS = Any[]; const SKB = Vector{UInt8}[]
    FAL.expand_sk(skb::Vector{UInt8}) =
        (push!(SKB, skb); ek = invoke(FAL.expand_sk, Tuple{AbstractVector{UInt8}}, skb); push!(EKS, ek); ek)
    check("mrm_bad_M1_rejected_before_decoding", threw(() -> F.falcon_mrm_sign([M1[1:end-1]; 12289], M2, sk)) && isempty(EKS))
    empty!(EKS)
    check("mrm_throw_wipes_expanded_key", threw(() -> F.falcon_mrm_sign(M1, Unreadable(), sk)) && length(EKS) == 1 && wiped(EKS[1]))
    empty!(EKS)
    check("sign_throw_wipes_expanded_key", threw(() -> F.falcon_sign(Unreadable(), sk)) && length(EKS) == 1 && wiped(EKS[1]))
    empty!(EKS); empty!(SKB)
    sig = try F.falcon_sign(M2, view([zeros(UInt8, 3); sk], Base.IdentityUnitRange(4:3+length(sk)))) catch; nothing end
    check("offset_sk_copy_wiped", sig !== nothing && F.falcon_verify(M2, sig, pk) && length(SKB) == 1 &&
                                  all(iszero, SKB[1]) && wiped(EKS[1]))

    # F13, F23, F31: record every seed, generator, target, candidate and encoding of one MRM
    # signature, with three forced rejections: by norm (attempt 1), at e2 once e1 is encoded
    # (attempt 2), at e1 (attempt 3); attempt 4 signs.
    const SEEDS = Vector{UInt8}[]; const GENS = CC.ChaCha20[]; const TGTS = Any[]; const CANDS = Any[]
    const ENCS = Vector{UInt8}[]; const NPRE = Ref(0); const NCOMP = Ref(0); const FAILSEED = Ref(0)
    function rb(k)
        k == FAL.SEED_LEN && length(SEEDS) + 1 == FAILSEED[] && error("randomness source failed")
        b = rand(UInt8, k); k == FAL.SEED_LEN && push!(SEEDS, b); b
    end
    CC.ChaCha20(seed::Vector{UInt8}) = (r = invoke(CC.ChaCha20, Tuple{AbstractVector{UInt8}}, seed); push!(GENS, r); r)
    FP.target(gs::FP.Setup, c::Vector{Int}) =
        (t = invoke(FP.target, Tuple{FP.Setup, AbstractVector{<:Integer}}, gs, c); push!(TGTS, t); t)
    function FAL.preimage(gs::FP.Setup, c::Vector{Int}, t0::Vector{FP.CF}, t1::Vector{FP.CF}, σmin::Float64, r::CC.ChaCha20)
        s = invoke(FAL.preimage, Tuple{FP.Setup, Any, Any, Any, Any, Any}, gs, c, t0, t1, σmin, r)
        push!(CANDS, s)
        (NPRE[] += 1) == 1 && (s[1] .*= 64; s[2] .*= 64)
        s
    end
    function FE.compress(s::Vector{Int}, slen::Int)
        (NCOMP[] += 1) in (2, 3) && return nothing
        e = invoke(FE.compress, Tuple{AbstractVector{<:Integer}, Int}, s, slen); push!(ENCS, e); e
    end
    ek = F.falcon_expand_sk(sk)
    sig = MRM.sign(M1, M2, ek; randombytes = rb)
    check("mrm_signature_valid", F.falcon_mrm_verify(M2, sig, pk) == M1)
    check("mrm_four_attempts", length(SEEDS) == length(GENS) == length(TGTS) == length(CANDS) == 4)
    check("mrm_seeds_wiped", all(s -> all(iszero, s), SEEDS))
    check("mrm_chacha_wiped", all(r -> all(iszero, r.s) && r.ctr == 0 && all(iszero, r.buf), GENS))
    check("mrm_targets_wiped", all(t -> zerocf(t[1]) && zerocf(t[2]), TGTS))
    check("mrm_rejected_candidates_wiped", all(s -> all(iszero, s[1]) && all(iszero, s[2]), CANDS[1:3]))
    check("mrm_rejected_e1_wiped", !isempty(ENCS) && all(iszero, ENCS[1]))
    # The randomness source fails on attempt 2's seed, after attempt 1 was rejected.
    foreach(empty!, (SEEDS, GENS, TGTS, CANDS, ENCS)); NPRE[] = 0; NCOMP[] = 0; FAILSEED[] = 2
    check("mrm_throw_after_rejection", threw(() -> MRM.sign(M1, M2, ek; randombytes = rb)) &&
                                       length(SEEDS) == 1 && length(TGTS) == 2)
    check("mrm_throw_wipes_targets", all(t -> zerocf(t[1]) && zerocf(t[2]), TGTS))
    check("mrm_throw_wipes_attempt1", all(iszero, SEEDS[1]) && all(iszero, GENS[1].buf) &&
                                      all(iszero, CANDS[1][1]) && all(iszero, CANDS[1][2]))
    """
    @testset "wipes observed in a child process (F10(a), F13, F23, F31)" begin
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $child`
        out = read(ignorestatus(cmd), String)
        got = Dict(m[1] => m[2] == "true" for m in eachmatch(r"^CHECK (\S+) (true|false)$"m, out))
        names = ("mrm_bad_M1_rejected_before_decoding", "mrm_throw_wipes_expanded_key",
                 "sign_throw_wipes_expanded_key", "offset_sk_copy_wiped", "mrm_signature_valid",
                 "mrm_four_attempts", "mrm_seeds_wiped", "mrm_chacha_wiped", "mrm_targets_wiped",
                 "mrm_rejected_candidates_wiped", "mrm_rejected_e1_wiped", "mrm_throw_after_rejection",
                 "mrm_throw_wipes_targets", "mrm_throw_wipes_attempt1")
        for name in names
            @test (name, get(got, name, false)) == (name, true)
        end
        all(name -> get(got, name, false), names) || println("child process output:\n", out)
    end
end
