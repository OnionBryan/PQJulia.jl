# test/regressions_falcon.jl — regressions for the Falcon core defects of the 2026 sweep:
#   F34          Base.vcat on byte data (its assertion is a runtime dispatch on Julia 1.12)
#   F11, F30     byte-key decoding left f, g, F, G mod q, their NTTs and products unwiped
#   F12, F24, F32 keygen left its seed, the ChaCha20 generator, f, g and their scratch unwiped
#   F27, F33, F10(b) a key refused by expansion, Falcon.sign(msg, sk::NamedTuple) and byte-key
#                signing that throws left the decoded key or the expanded key unwiped
#   F35          keygen's bare catch around the NTRU solver swallowed every exception
# The "no copy is left" checks run the operations in a child process with the GC disabled, so that
# nothing is freed or reused, and scan that process's memory for each secret value; the scanner
# holds each value only XOR-masked. Run alone with
#   julia --project=. -e 'using Test, PQJulia; include("test/regressions_falcon.jl")'
using Test, PQJulia, Random

module FalconSweepRegressions

using PQJulia, Random
const FAL = PQJulia.FNDSA.Falcon; const NG = PQJulia.FNDSA.FalconNTRUGen; const FX = PQJulia.FNDSA.FalconFxp
const FP = PQJulia.FNDSA.FalconFpr; const FC = PQJulia.FNDSA.FalconCertify; const FF = PQJulia.FNDSA.FalconFFT
const FS = PQJulia.FNDSA.FalconSampler; const CC = PQJulia.FNDSA.FalconChaCha; const KEC = PQJulia.Keccak
const q = 12289

# Deterministic randomness: counter-mode SHAKE256, a fresh vector per call.
det(tag) = (c = Ref(0); k -> (c[] += 1; KEC.shake256(Vector{UInt8}(string(tag, "/", c[])), UInt64(k))))

# A valid encoding of a degree-16 key whose ffLDL leaves lie outside [σmin, σmax]: it decodes, but
# key expansion refuses it (as in test/falcon_certify.jl).
function refused_key16()
    rng = MersenneTwister(1)
    for _ in 1:100_000
        f = rand(rng, -15:15, 16); g = rand(rng, -15:15, 16)
        FAL.is_invertible(f) || continue
        FG = try NG.ntru_solve(BigInt.(f), BigInt.(g)) catch e; e isa InterruptException && rethrow(); nothing end
        FG === nothing && continue
        F, G = FG
        FAL.reduce_FG!(BigInt.(f), BigInt.(g), F, G)
        (NG.ntru_check(BigInt.(f), BigInt.(g), F, G) && FAL.fits(F, 8) && FAL.fits(G, 8)) || continue
        skb = FAL.encode_sk(FAL.secret_key(f, g, F, G))
        k = FAL.decode_sk(skb); k === nothing && continue
        refused = try FAL.wipe!(FAL.expand_sk(k)); false catch e; e isa ArgumentError || rethrow(); true end
        refused && return skb, k
    end
    error("no refused degree-16 key found")
end

# ── masked patterns ─────────────────────────────────────────────────────────
# A pattern is a run of 8-byte words P, stored as P ⊻ M with M random. All patterns share their
# first mask word M0, so the scanner tests each memory word against one table.
struct Pat
    name::String
    pm::Vector{UInt64}
    m::Vector{UInt64}
end
words(v::AbstractVector, r) = collect(reinterpret(UInt64, collect(v[r])))
function mkpat(name, P::Vector{UInt64}, M0::UInt64)
    m = rand(RandomDevice(), UInt64, length(P)); m[1] = M0
    Pat(name, P .⊻ m, m)
end
# A 16-coefficient window with a negative entry (so that v mod q differs from v), as words.
function window(v::AbstractVector{<:Integer}, len = 16)
    for s in 1:len:length(v)-len+1
        any(<(0), v[s:s+len-1]) && return s:s+len-1
    end
    1:len
end

# Writes the job for the child: secret keys and keygen seeds (seeds masked), and the patterns.
function write_job(path)
    M0 = rand(RandomDevice(), UInt64)
    pats = Pat[]; lines = String[]
    add(name, P) = push!(pats, mkpat(name, P, M0))
    key(name, b) = push!(lines, "KEY $name $(bytes2hex(b))")
    function seed(name, s)
        mk = rand(RandomDevice(), UInt8, length(s))
        push!(lines, "SEED $name $(bytes2hex(s .⊻ mk)) $(bytes2hex(mk))")
    end
    add("control (materialized by the scanner)", rand(UInt64, 16))

    # F11/F30: Falcon512.falcon_sign(msg, skb) decodes the key for one call.
    pk, skb = FAL.keypair(512; randombytes = det("decode512")); key("decode512", skb); key("decode512pk", pk)
    k = FAL.decode_sk(skb); r = window(k.f)
    add("decode: f", words(k.f, r)); add("decode: G", words(k.G, window(k.G)))
    add("decode: f mod q", words(mod.(k.f, q), r)); add("decode: g mod q", words(mod.(k.g, q), window(k.g)))
    add("decode: F mod q", words(mod.(k.F, q), window(k.F))); add("decode: G mod q", words(mod.(k.G, q), window(k.G)))
    add("decode: NTT(f)", words(FAL.ntt_fwd(mod.(k.f, q)), 201:216))
    add("decode: NTT(g)", words(FAL.ntt_fwd(mod.(k.g, q)), 201:216))
    add("decode: NTT(F)", words(FAL.ntt_fwd(mod.(k.F, q)), 201:216))
    add("decode: f*G", words(FAL.negamul_int(k.f, k.G), 201:216))
    add("decode: g*F", words(FAL.negamul_int(k.g, k.F), 201:216))
    add("decode: g*F mod q", words(FAL.poly_mul_modq(k.g, k.F), 201:216))
    add("decode: FFT basis b01 (expanded key)", words(FP.setup(k, FAL.params(512).σ, FAL.params(512).σmin).b01, 101:108))

    # Certificate of a decoded key: the Int copies and the mod-p transforms of f, g, F, G.
    _, skb = FAL.keypair(64; randombytes = det("cert64")); key("cert64", skb)
    k = FAL.decode_sk(skb); p1 = FC.ntt_primes(64, 1)[1]
    add("certify: g", words(k.g, window(k.g))); add("certify: -f", words(.-k.f, window(.-k.f)))
    add("certify: G", words(k.G, window(k.G)))
    add("certify: NTT mod p of g", words(FC.nttm(FC.ModRing(64, p1), k.g), 1:16))

    # F27/F33/F10(b): a key that decodes but that expansion refuses.
    skb, k = refused_key16(); key("refused16", skb)
    add("refused key: f", words(k.f, 1:16)); add("refused key: G", words(k.G, 1:16))

    # F33: Falcon.sign(msg, sk::NamedTuple) expands the key for one call.
    _, skb = FAL.keypair(64; randombytes = det("nt64")); key("nt64", skb)
    k = FAL.decode_sk(skb); gs = FP.setup(k, FAL.params(64).σ, FAL.params(64).σmin)
    add("sign(NamedTuple): FFT basis b01", words(gs.b01, 1:8)); add("sign(NamedTuple): ffLDL root l10", words(gs.T.l10, 1:8))

    # F10(b): Falcon.sign(msg, skb) where signing throws after the key is expanded.
    _, skb = FAL.keypair(64; randombytes = det("throw64")); key("throw64", skb)
    k = FAL.decode_sk(skb)
    add("sign throws: f", words(k.f, window(k.f)))
    add("sign throws: FFT basis b01", words(FP.setup(k, FAL.params(64).σ, FAL.params(64).σmin).b01, 1:8))

    # F12/F24/F32: keygen from a known seed; the child runs the same keypair call.
    s = det("kg64")(56); seed("kg64", s)
    k = FAL.keygen(64; randombytes = _ -> copy(s))
    r0 = CC.ChaCha20(s); raw = [FS.samplerz(0.0, FAL.SIGMA_FG, FAL.SIGMA_FG - 0.001, r0) for _ in 1:4096]
    add("keygen: seed (and ChaCha20 state)", words(s, 1:56)); add("keygen: f", words(k.f, window(k.f)))
    add("keygen: g", words(k.g, window(k.g))); add("keygen: f mod q", words(mod.(k.f, q), window(k.f)))
    add("keygen: NTT(f)", words(FAL.ntt_fwd(mod.(k.f, q)), 1:16))
    add("keygen: Float64 f", words(Float64.(k.f), 1:16)); add("keygen: FFT(f)", words(FF.fft(Float64.(k.f)), 1:8))
    add("keygen: raw samples of gen_poly", words(raw, 1:16))
    s = det("kgfxp64")(56); seed("kgfxp64", s)
    k = FAL.keygen(64; fixedpoint = true, randombytes = _ -> copy(s))
    add("keygen fixedpoint: seed", words(s, 1:56)); add("keygen fixedpoint: f", words(k.f, window(k.f)))
    add("keygen fixedpoint: g", words(k.g, window(k.g)))
    add("keygen fixedpoint: keystream chunk", words(CC.randbytes(CC.ChaCha20(s), 512), 1:64))
    s = det("kgcert64")(56); seed("kgcert64", s)
    k = FAL.keygen(64; certified = true, randombytes = _ -> copy(s))
    add("keygen certified: seed", words(s, 1:56)); add("keygen certified: F", words(k.F, window(k.F)))
    add("keygen certified: G", words(k.G, window(k.G)))

    push!(lines, "M0 $(string(M0, base = 16))")
    for p in pats
        push!(lines, "PAT $(bytes2hex(reinterpret(UInt8, p.pm))) $(bytes2hex(reinterpret(UInt8, p.m))) $(p.name)")
    end
    write(path, join(lines, "\n") * "\n")
    [p.name for p in pats]
end

# The child: runs every scenario with the GC disabled, scans its writable anonymous memory, then
# checks keygen's handling of solver exceptions with an injected ntru_solve method.
const CHILD = raw"""
using PQJulia
const FAL = PQJulia.FNDSA.Falcon; const NG = PQJulia.FNDSA.FalconNTRUGen; const F5 = FNDSA.Falcon512
KEYS = Dict{String,String}(); SEEDS = Dict{String,Tuple{Vector{UInt8},Vector{UInt8}}}()
PATS = Tuple{String,Vector{UInt64},Vector{UInt64}}[]; M0 = UInt64(0)
for line in eachline(ARGS[1])
    t = split(line, ' '; limit = 4)
    if t[1] == "KEY"; KEYS[t[2]] = t[3]
    elseif t[1] == "SEED"; SEEDS[t[2]] = (hex2bytes(t[3]), hex2bytes(split(t[4])[1]))
    elseif t[1] == "M0"; global M0 = parse(UInt64, t[2]; base = 16)
    elseif t[1] == "PAT"
        push!(PATS, (String(t[4]), collect(reinterpret(UInt64, hex2bytes(t[2]))), collect(reinterpret(UInt64, hex2bytes(t[3])))))
    end
end
key(name) = hex2bytes(KEYS[name])
# randombytes returning the keygen seed: the only clear copy is the returned vector.
function seedfn(name)
    sm, mk = SEEDS[name]
    k -> (s = copy(sm); s .⊻= mk; s)
end
outcome(f) = try f(); "returned" catch e; string(nameof(typeof(e))) end
msg = Vector{UInt8}("PQJulia sweep regression")

GC.enable(false)
println("RESULT decode512-signs ", F5.falcon_verify(msg, F5.falcon_sign(msg, key("decode512")), key("decode512pk")))
let k = FAL.decode_sk(key("cert64")); FAL.certify(k); FAL.wipe_sk!(k) end
println("RESULT refused16 ", outcome(() -> FAL.expand_sk(key("refused16"))))
println("RESULT refused16 ", outcome(() -> FAL.sign(msg, key("refused16"))))
let k = FAL.decode_sk(key("nt64")); FAL.sign(msg, k); FAL.wipe_sk!(k) end
nosampler(n) = n == FAL.SALT_LEN ? zeros(UInt8, n) : error("no sampler seed")
println("RESULT throw64 ", outcome(() -> FAL.sign(msg, key("throw64"); randombytes = nosampler)))
for (name, kw) in (("kg64", (;)), ("kgfxp64", (; fixedpoint = true)), ("kgcert64", (; certified = true)))
    local pk, sk = FAL.keypair(64; randombytes = seedfn(name), kw...)
    wipe!(sk)
end

# Read buffer (1 MiB): bytes for readbytes!, the same memory as words for the comparisons.
const BUF8 = Vector{UInt8}(undef, 1 << 20)
const BUF = unsafe_wrap(Array, Ptr{UInt64}(pointer(BUF8)), length(BUF8) ÷ 8)
function regions()
    out = Tuple{UInt64,UInt64}[]
    for line in eachline("/proc/self/maps")
        p = split(line)
        startswith(p[2], "rw") || continue
        nm = length(p) >= 6 ? p[6] : ""
        (startswith(nm, "/") || nm in ("[vvar]", "[vvar_vclock]", "[vsyscall]")) && continue   # file-backed: no heap
        lo, hi = (parse(UInt64, s; base = 16) for s in split(p[1], '-'))
        push!(out, (lo, hi))
    end
    out
end
function scan(pats, m0)
    heads = Dict{UInt64,Vector{Int}}(); bloom = falses(1 << 16)
    for (i, (_, pm, _)) in enumerate(pats)
        push!(get!(heads, pm[1], Int[]), i); bloom[(pm[1] & 0xffff) + 1] = true
    end
    hits = [Set{UInt64}() for _ in pats]; L = maximum(p -> length(p[2]), pats)
    lo_buf = UInt64(pointer(BUF)); hi_buf = lo_buf + sizeof(BUF)
    io = open("/proc/self/mem")
    for (lo, hi) in regions()
        addr = lo
        while addr < hi
            nb = Int(min(hi - addr, UInt64(sizeof(BUF))))
            got = try seek(io, Int64(addr)); readbytes!(io, BUF8, nb) catch; 0 end
            nw = got ÷ 8
            @inbounds for j in 1:nw
                x = BUF[j] ⊻ m0
                bloom[(x & 0xffff) + 1] || continue
                ids = get(heads, x, nothing); ids === nothing && continue
                a = addr + 8(j - 1)
                lo_buf <= a < hi_buf && continue
                for i in ids
                    _, pm, m = pats[i]
                    j + length(pm) - 1 <= nw || continue
                    ok = true
                    for t in 2:length(pm); BUF[j+t-1] ⊻ m[t] == pm[t] || (ok = false; break); end
                    ok && push!(hits[i], a)
                end
            end
            nb <= 8L && break
            addr += UInt64(nb - 8L)                       # overlap: a pattern across chunks is seen
        end
    end
    close(io)
    length.(hits)
end
ctl = PATS[1][2] .⊻ PATS[1][3]                            # positive control, in the clear
counts = GC.@preserve ctl BUF8 scan(PATS, M0)
for ((name, _, _), c) in zip(PATS, counts); println("HIT ", c, " ", name); end
GC.enable(true)

# F35: an exception other than the solver's failure signal propagates out of keygen; the failure
# signal itself still makes keygen draw a new candidate.
const MODE = Ref(:interrupt); const CALLS = Ref(0)
NG.ntru_solve(f::Vector{BigInt}, g::Vector{BigInt}) =
    MODE[] === :interrupt ? throw(InterruptException()) :
    (CALLS[] += 1; CALLS[] <= 2 ? throw(NG.NTRUSolveFailure("injected")) : invoke(NG.ntru_solve, Tuple{Any,Any}, f, g))
# Enough tries that candidates reach the solver; before the fix every one was swallowed and keygen
# ended with "keygen failed in 200 tries".
for fixedpoint in (false, true)
    println("RESULT interrupt-fixedpoint=$fixedpoint ", outcome(() -> FAL.keygen(16; fixedpoint, maxtries = 200)))
end
MODE[] = :failure
for fixedpoint in (false, true)
    CALLS[] = 0
    r = try string(FAL.keygen(16; fixedpoint, maxtries = 200).n == 16 && CALLS[] > 2) catch e; string(nameof(typeof(e))) end
    println("RESULT resample-fixedpoint=$fixedpoint ", r)
end
println("DONE")
"""

function run_child()
    dir = mktempdir()
    try
        job = joinpath(dir, "job.txt"); script = joinpath(dir, "child.jl")
        names = write_job(job); write(script, CHILD)
        proj = Base.active_project()
        cmd = proj === nothing ? `$(Base.julia_cmd()) --startup-file=no $script $job` :
                                 `$(Base.julia_cmd()) --startup-file=no --project=$proj $script $job`
        err = IOBuffer()
        out = read(pipeline(ignorestatus(cmd); stderr = err), String)
        occursin("DONE", out) || println(stderr, "child process output:\n", out, "\n", String(take!(err)))
        hits = Dict{String,Int}(); results = Dict{String,String}()
        for line in split(out, '\n')
            t = split(line, ' '; limit = 3)
            length(t) == 3 || continue
            t[1] == "HIT" && (hits[t[3]] = parse(Int, t[2]))
            t[1] == "RESULT" && (results[t[2]] = haskey(results, t[2]) ? results[t[2]] * "," * t[3] : String(t[3]))
        end
        names, hits, results, occursin("DONE", out)
    finally
        rm(dir; recursive = true, force = true)
    end
end

end # module FalconSweepRegressions

@testset "Falcon regressions (sweep)" begin
    R = FalconSweepRegressions
    FAL = R.FAL; NG = R.NG; FX = R.FX; CC = R.CC
    msg = Vector{UInt8}("PQJulia sweep regression")

    @testset "F34: the Falcon core concatenates bytes without Base.vcat" begin
        dir = joinpath(pkgdir(PQJulia), "src", "falcon")
        for f in ("falcon.jl", "falcon_ntrugen.jl", "falcon_fxp.jl", "falcon_certify.jl", "falcon_sampler.jl",
                  "falcon_chacha.jl", "falcon_fpr.jl", "falcon_fft.jl", "falcon_encoding.jl")
            @test !occursin("vcat(", read(joinpath(dir, f), String))
        end
        salt = collect(UInt8, 1:40)
        @test FAL.encode_sig((; salt, enc = UInt8[7, 8]), 512) == UInt8[0x39; salt; 7; 8]
        @test FAL.hash_to_point(msg, salt, 64) == FAL.hash_to_point(msg, view(salt, 1:40), 64)
    end

    @testset "F12/F24/F32: keygen wipes its seed on every exit" begin
        rec = Vector{UInt8}[]
        rb = k -> (b = rand(RandomDevice(), UInt8, k); push!(rec, b); b)
        for (n, kw) in ((64, (;)), (64, (; fixedpoint = true)), (64, (; certified = true)), (512, (;)))
            empty!(rec)
            pk, sk = FAL.keypair(n; randombytes = rb, kw...)
            @test length(rec) == 1 && length(rec[1]) == FAL.SEED_LEN
            @test all(iszero, rec[1])
            @test FAL.verify(msg, FAL.sign(msg, sk), pk)
        end
        empty!(rec)
        @test_throws ErrorException FAL.keygen(64; randombytes = rb, maxtries = 0)
        @test length(rec) == 1 && all(iszero, rec[1])
        empty!(rec)
        @test_throws ArgumentError FAL.keygen(64; randombytes = k -> rb(k - 1))   # 55-byte seed
        @test length(rec) == 1 && all(iszero, rec[1])
    end

    @testset "signing wipes its sampler seed when the generator cannot be built" begin
        pk, skb = FAL.keypair(64; randombytes = R.det("seedwipe"))
        ek = FAL.expand_sk(skb)
        seeds = Vector{UInt8}[]
        rb = k -> (b = rand(RandomDevice(), UInt8, k == FAL.SEED_LEN ? k - 1 : k); k == FAL.SEED_LEN && push!(seeds, b); b)
        @test_throws ArgumentError FAL.sign(msg, ek; randombytes = rb)
        @test length(seeds) == 1 && all(iszero, seeds[1])
        empty!(seeds)
        rb2 = k -> (b = rand(RandomDevice(), UInt8, k); k == FAL.SEED_LEN && push!(seeds, b); b)
        @test FAL.verify(msg, FAL.sign(msg, ek; randombytes = rb2), pk)
        @test !isempty(seeds) && all(s -> all(iszero, s), seeds)
        FAL.wipe!(ek)
    end

    @testset "F33: Falcon.sign(msg, sk::NamedTuple) leaves the caller's key usable" begin
        pk, skb = FAL.keypair(64; randombytes = R.det("ntkey"))
        k = FAL.decode_sk(skb); k0 = deepcopy(k)
        @test FAL.verify(msg, FAL.sign(msg, k), pk)
        @test k == k0
        @test FAL.verify(msg, FAL.sign(msg, k), pk)
        FAL.wipe_sk!(k)
    end

    @testset "F35: the NTRU solver signals a failed solve explicitly" begin
        @test NG.NTRUSolveFailure <: Exception
        @test_throws NG.NTRUSolveFailure NG.ntru_solve(BigInt[2, 0], BigInt[2, 0])          # resultants 4, 4
        @test_throws NG.NTRUSolveFailure FX.ntru_solve(BigInt[2, 0, 0, 0], BigInt[2, 0, 0, 0])
        @test_throws NG.NTRUSolveFailure FX.to_fxr(big(2)^80, 0)                          # outside the fxr range
        @test FX.to_fxr(big(3), 1) == Int64(3) << 31
        @test occursin("non-coprime", sprint(showerror, NG.NTRUSolveFailure("non-coprime resultants")))
    end

    @testset "review fixes: ranges, offset axes, wrong-degree keys, NTRUSolveFailure()" begin
        FF = R.FF; FC = R.FC; FE = PQJulia.FNDSA.FalconEncoding
        @test FF.fft(1:8) == FF.fft(collect(1:8)) && FF.fft(0.0:7.0) == FF.fft(collect(0.0:7.0))
        P = FC.ModRing(8, FC.ntt_primes(8, 1)[1])
        @test FC.nttm(P, 0:7) == FC.nttm(P, collect(0:7))
        @test NG.NTRUSolveFailure() isa NG.NTRUSolveFailure
        v = view([7, 1, 2, 3, 4], Base.IdentityUnitRange(2:5)); w = view([9, 5, 6, 7, 8], Base.IdentityUnitRange(2:5))
        @test FAL.poly_div_modq(v, [1, 0, 0, 0]) == [1, 2, 3, 4]
        @test FAL.poly_mul_modq(v, w) == FAL.poly_mul_modq([1, 2, 3, 4], [5, 6, 7, 8])
        @test FAL.is_invertible(v) == FAL.is_invertible([1, 2, 3, 4])
        @test_throws ArgumentError FAL.negamul_int(v, w)
        pk, skb = FAL.keypair(64; randombytes = R.det("degree"))
        k = FAL.decode_sk(skb); sig = FAL.sign_poly(k, FAL.expand_sk(k).gs, msg)
        @test FAL.verify_poly(k.h, 64, msg, sig)
        @test !FAL.verify_poly([k.h; k.h], 64, msg, sig) && !FAL.verify_poly(k.h[1:32], 64, msg, sig)
        FAL.wipe_sk!(k)
        # A rejected compression leaves no partial encoding in the output buffer.
        out = fill(0xff, 2)
        @test FE.compress(fill(300, 4), 16, out) === nothing && all(iszero, out)
        @test FE.compress([1, -2, 3], 32, fill(0xff, 4)) == FE.compress([1, -2, 3], 32)
    end

    # Signing and decoding wipes observed through recording methods in a child process: more
    # specific methods that forward with `invoke`, so the test process keeps the package's own.
    @testset "signing and decoding wipes (child process)" begin
        spy = raw"""
        using PQJulia
        const FAL = PQJulia.FNDSA.Falcon; const FP = PQJulia.FNDSA.FalconFpr; const KEC = PQJulia.Keccak
        det(tag) = (c = Ref(0); k -> (c[] += 1; KEC.shake256(Vector{UInt8}(string(tag, "/", c[])), UInt64(k))))
        check(name, ok) = println("CHECK ", name, " ", ok)
        const TGT = Any[]; const PRE = Any[]; const DEC = Any[]; const FFS = Any[]; const KZ = Any[]
        FP.target(gs::FP.Setup, c::Vector{Int}) =
            (t = invoke(FP.target, Tuple{FP.Setup, AbstractVector{<:Integer}}, gs, c); push!(TGT, t); t)
        FP.preimage(gs::FP.Setup, c::Vector{Int}, t0, t1, σmin::Float64, r) =
            (s = invoke(FP.preimage, Tuple{FP.Setup, AbstractVector{<:Integer}, Any, Any, Float64, Any}, gs, c, t0, t1, σmin, r); push!(PRE, s); s)
        FAL.decode_signed(bytes::SubArray, n::Int, bits::Int) =
            (s = invoke(FAL.decode_signed, NTuple{3, Any}, bytes, n, bits); push!(DEC, s); s)
        FAL.ffsample(t0::Vector{ComplexF64}, t1::Vector{ComplexF64}, T, σmin::Float64, r) =
            (z = invoke(FAL.ffsample, NTuple{5, Any}, t0, t1, T, σmin, r); push!(FFS, z); z)
        FAL.klein_sample(B::Matrix{Float64}, Bgs, gsnorm2, target, σ, σmin, r) =
            (z = invoke(FAL.klein_sample, NTuple{7, Any}, B, Bgs, gsnorm2, target, σ, σmin, r); push!(KZ, z); z)
        nz(v) = any(x -> x != zero(x), v)
        fprnz(v) = any(x -> x.re.b != 0 || x.im.b != 0, v)
        function main()
            msg = Vector{UInt8}("spy")
            # Signing throws after the target is computed: t0, t1 are wiped.
            pk, skb = FAL.keypair(64; randombytes = det("t64")); ek = FAL.expand_sk(skb)
            r = try FAL.sign(msg, ek; randombytes = k -> k == FAL.SEED_LEN ? error("no seed") : zeros(UInt8, k)); "returned" catch e; string(typeof(e)) end
            check("throw-after-target", r == "ErrorException" && !isempty(TGT))
            check("throw-target-zero", !fprnz(TGT[end][1]) && !fprnz(TGT[end][2]))
            # Degree 2 rejects often: every rejected (s0, s1) is wiped.
            pk2, skb2 = FAL.keypair(2; randombytes = det("n2")); ek2 = FAL.expand_sk(skb2)
            rej = 0; left = 0
            for i in 1:200
                empty!(PRE); FAL.sign(Vector{UInt8}("m$i"), ek2)
                for (s0, s1) in PRE[1:end-1]; rej += 1; (nz(s0) || nz(s1)) && (left += 1); end
            end
            check("rejections-seen", rej > 0)
            check("rejected-candidates-zero", left == 0)
            # A key whose F field does not decode: the decoded f and g are wiped.
            pk5, skb5 = FAL.keypair(512; randombytes = det("p512")); bad = copy(skb5); bad[end] = 0x80
            empty!(DEC)
            check("undecodable-rejected", FAL.decode_sk(bad) === nothing && length(DEC) == 3 && DEC[3] === nothing)
            check("undecodable-fields-zero", !nz(DEC[1]) && !nz(DEC[2]))
            # Float64 reference signer: every level's ffSampling output is wiped.
            k64 = FAL.decode_sk(skb); gsf = FAL.sign_setup(k64)
            empty!(FFS); sig = FAL.sign_poly(k64, gsf, msg; randombytes = det("f64"))
            check("float-signer-verifies", FAL.verify_poly(k64.h, 64, msg, sig))
            check("float-ffsample-zero", !isempty(FFS) && all(z -> !nz(z[1]) && !nz(z[2]), FFS))
            # Klein reference signer: every lattice sample z is wiped.
            pk16, skb16 = FAL.keypair(16; randombytes = det("k16")); k16 = FAL.decode_sk(skb16)
            gk = FAL.sign_setup_klein(k16)
            empty!(KZ); sigk = FAL.sign_poly_klein(k16, gk, msg; randombytes = det("kl"))
            check("klein-verifies", FAL.verify_poly(k16.h, 16, msg, sigk))
            check("klein-z-zero", !isempty(KZ) && all(z -> !nz(z), KZ))
            println("DONE")
        end
        main()
        """
        dir = mktempdir()
        try
            script = joinpath(dir, "spy.jl"); write(script, spy)
            proj = Base.active_project()
            cmd = proj === nothing ? `$(Base.julia_cmd()) --startup-file=no $script` :
                                     `$(Base.julia_cmd()) --startup-file=no --project=$proj $script`
            err = IOBuffer()
            out = read(pipeline(ignorestatus(cmd); stderr = err), String)
            occursin("DONE", out) || println(stderr, "spy child output:\n", out, "\n", String(take!(err)))
            got = Dict(String(t[2]) => t[3] == "true" for t in (split(l, ' ') for l in split(out, '\n')) if length(t) == 3 && t[1] == "CHECK")
            for name in ("throw-after-target", "throw-target-zero", "rejections-seen", "rejected-candidates-zero",
                         "undecodable-rejected", "undecodable-fields-zero", "float-signer-verifies",
                         "float-ffsample-zero", "klein-verifies", "klein-z-zero")
                @test (name => get(got, name, false)) == (name => true)
            end
        finally
            rm(dir; recursive = true, force = true)
        end
    end

    if Sys.islinux() && isfile("/proc/self/mem")
        names, hits, results, done = R.run_child()
        @testset "child process ran every scenario" begin
            @test done
            @test get(results, "decode512-signs", "") == "true"
            @test get(results, "refused16", "") == "ArgumentError,ArgumentError"
            @test get(results, "throw64", "") == "ErrorException"
        end
        @testset "F11/F30, F12/F24/F32, F27/F33/F10(b): no copy of the key is left in memory" begin
            @test get(hits, names[1], 0) >= 1                         # the scanner finds a planted value
            for name in names[2:end]
                @test (name => get(hits, name, -1)) == (name => 0)
            end
        end
        @testset "F35: keygen lets every exception but the solver's failure through" begin
            @test get(results, "interrupt-fixedpoint=false", "") == "InterruptException"
            @test get(results, "interrupt-fixedpoint=true", "") == "InterruptException"
            @test get(results, "resample-fixedpoint=false", "") == "true"
            @test get(results, "resample-fixedpoint=true", "") == "true"
        end
    else
        @test_skip "memory scan needs Linux /proc/self/mem"
    end
end
