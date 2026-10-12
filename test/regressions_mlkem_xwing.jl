# Regression tests for the ML-KEM, X-Wing and X25519 findings of the 2026-10 sweep:
#   F5/F15  the basemul scratch (Â∘ŝ, t̂∘r̂, ŝ∘NTT(u)) was never wiped;
#   F6/F15  decapsulation never wiped the re-encryption of m′, and no wipe ran on an exception;
#   F7      Base.vcat on byte data failed the runtime-dispatch gate (test/timing/dispatch.jl);
#   F8      xwing_encaps_derand threw MethodError for a Memory or range eseed;
#   F9      ek/dk checks, X25519 and X-Wing indexed AbstractVector inputs from 1;
#   F16     X-Wing encapsulation left ekX, m, ssX (and xwing_encaps its eseed) on the error path.
# Included from runtests.jl (`using Test, PQJulia`), or on its own:
#   julia --project=. -e 'using Test, PQJulia; include("test/regressions_mlkem_xwing.jl")'
# The wipe checks spy on internal helpers by redefining them, so they run in a child process and
# nothing leaks into the other testsets.

# A minimal 0-based AbstractVector, the shape of an OffsetArrays.OffsetVector (not a dependency).
struct ZeroBasedBytes <: AbstractVector{UInt8}
    d::Vector{UInt8}
end
Base.size(v::ZeroBasedBytes) = size(v.d)
Base.axes(v::ZeroBasedBytes) = (Base.IdentityUnitRange(0:length(v.d)-1),)
Base.@propagate_inbounds Base.getindex(v::ZeroBasedBytes, i::Int) =
    (@boundscheck checkbounds(v, i); v.d[i+1])

# Child process: spies record the buffers each step allocates, then the parent asserts on the
# printed `CHECK <true|false> <name>` lines. Every spy reproduces the original behaviour, and the
# round trips are checked, so a spy cannot hide a wrong result.
const MLKEM_XWING_WIPE_CHILD = raw"""
using PQJulia
const M = PQJulia.MLKEM; const X = PQJulia.X25519; const XW = PQJulia.XWing
const LEVELS = (M.Category1, M.Category3, M.Category5)
check(name, ok) = println("CHECK ", ok === true, " ", name)
zeroed(bufs) = all(b -> all(iszero, b), bufs)

# The third argument of every polynomial add/sub on these paths is secret: the basemul scratch,
# e, e1, e2, the message polynomial, or ŝᵀu.
const ADDENDS = Vector{Vector{Int16}}()
function M.kyber_poly_add!(r::Vector{Int16}, a::Vector{Int16}, b::Vector{Int16})
    push!(ADDENDS, b)
    @inbounds for i in 1:256; r[i] = a[i] + b[i]; end
    r
end
function M.kyber_poly_sub!(r::Vector{Int16}, a::Vector{Int16}, b::Vector{Int16})
    push!(ADDENDS, b)
    @inbounds for i in 1:256; r[i] = a[i] - b[i]; end
    r
end
# The second argument of the FO comparison is the re-encryption of m′.
const CMPS = Vector{Vector{UInt8}}()
function M.kyber_verify(a::Vector{UInt8}, b::Vector{UInt8})
    push!(CMPS, b)
    d = 0x00
    for i in eachindex(a, b); d |= a[i] ⊻ b[i]; end
    Int(d != 0x00)
end
# G's input and output: d ‖ k and ρ ‖ σ, or m ‖ H(ek) and (K, r).
const GCALLS = Tuple{Vector{UInt8},Vector{UInt8}}[]
function M.kyber_hash_g(data::Vector{UInt8})
    out = PQJulia.Keccak.sha3_512(data); push!(GCALLS, (data, out)); out
end
# X25519 scalars are secret; outputs are secret unless u is the base point.
const XCALLS = Tuple{Vector{UInt8},Vector{UInt8},Vector{UInt8}}[]
function X.x25519(k::Vector{UInt8}, u::Vector{UInt8})
    r = invoke(X.x25519, Tuple{AbstractVector{UInt8}, AbstractVector{UInt8}}, k, u)
    push!(XCALLS, (k, u, r)); r
end
x25519_clean() = all(((k, u, r),) -> all(iszero, k) && (u === X.BASE || all(iszero, r)), XCALLS)
const EXPANDED = Any[]
function XW.expand(sk::Vector{UInt8})
    k = invoke(XW.expand, Tuple{AbstractVector{UInt8}}, sk); push!(EXPANDED, k); k
end
const ESEEDS = Vector{Vector{UInt8}}()
function XW.xwing_encaps_derand(pk::Vector{UInt8}, eseed::Vector{UInt8})
    push!(ESEEDS, eseed)
    invoke(XW.xwing_encaps_derand, Tuple{AbstractVector{UInt8}, AbstractVector{UInt8}}, pk, eseed)
end
reset() = foreach(empty!, (ADDENDS, CMPS, GCALLS, XCALLS, EXPANDED, ESEEDS))

# ── F5/F15, F6/F15: every secret buffer is zero after keygen, encaps and decaps ──
const KEYS = Dict{Any,Any}()
for C in LEVELS
    id = C.IDENTIFIER
    reset(); pk, sk = C.kyber_kem_keypair_derand(collect(UInt8, 1:64)); KEYS[C] = (pk, sk)
    check("$id keygen: basemul scratch and e zeroed", !isempty(ADDENDS) && zeroed(ADDENDS))
    check("$id keygen: G input and output zeroed", length(GCALLS) == 1 && zeroed(GCALLS[1]))
    reset(); ct, ss = C.kyber_kem_enc_derand(pk, fill(0x42, 32))
    check("$id encaps: basemul scratch, e1, e2, message polynomial zeroed", !isempty(ADDENDS) && zeroed(ADDENDS))
    check("$id encaps: G input and output zeroed", length(GCALLS) == 1 && zeroed(GCALLS[1]))
    bad = copy(ct); bad[end] ⊻= 0xff                         # decrypts to m′ ≠ m: rejected
    for (label, c, want) in (("valid ct", ct, true), ("rejected ct", bad, false))
        reset(); out = C.kyber_kem_dec(c, sk)
        check("$id decaps ($label): shared secret correct", (out == ss) == want)
        check("$id decaps ($label): basemul scratch ŝ∘NTT(u) and re-encryption temporaries zeroed",
              !isempty(ADDENDS) && zeroed(ADDENDS))
        check("$id decaps ($label): re-encryption cmp zeroed", length(CMPS) == 1 && zeroed(CMPS))
        check("$id decaps ($label): m′ ‖ H(ek) and (K′, r′) zeroed", length(GCALLS) == 1 && zeroed(GCALLS[1]))
    end
end
reset(); xpk, xsk = XW.xwing_keypair_derand(collect(UInt8, 1:32))
check("X-Wing keypair: expanded skM, skX and ML-KEM temporaries zeroed",
      zeroed((EXPANDED[1].skM, EXPANDED[1].skX)) && zeroed(ADDENDS) && x25519_clean())
reset(); xct, xss = XW.xwing_encaps_derand(xpk, collect(UInt8, 65:128))
check("X-Wing encaps: ekX, ssX and ML-KEM temporaries zeroed", !isempty(XCALLS) && x25519_clean() && zeroed(ADDENDS))
reset(); xout = XW.xwing_decaps(xct, xsk)
check("X-Wing decaps: shared secret correct", xout == xss)
check("X-Wing decaps: skM, skX, ssX, cmp and ML-KEM temporaries zeroed",
      zeroed((EXPANDED[1].skM, EXPANDED[1].skX)) && x25519_clean() && zeroed(ADDENDS) && zeroed(CMPS))

# ── F16: an encapsulation key that fails the FIPS 203 check leaves no secret behind ──
xbad = copy(xpk); xbad[1] = 0xff; xbad[2] = 0xff             # first ML-KEM coefficient 4095 ≥ q
reset(); eseed = collect(UInt8, 101:164)
threw = try XW.xwing_encaps_derand(xbad, eseed); false catch e; e isa ArgumentError end
check("X-Wing encaps_derand(bad pk): ArgumentError", threw)
check("X-Wing encaps_derand(bad pk): no ekX or ssX left", x25519_clean())
check("X-Wing encaps_derand(bad pk): caller's eseed untouched", eseed == collect(UInt8, 101:164))
reset()
threw = try XW.xwing_encaps(xbad); false catch e; e isa ArgumentError end
check("X-Wing encaps(bad pk): ArgumentError", threw)
check("X-Wing encaps(bad pk): eseed zeroed", length(ESEEDS) == 1 && zeroed(ESEEDS))

# ── F6/F15, F16: an exception after the secrets exist still runs every wipe ──
function M.kyber_verify(a::Vector{UInt8}, b::Vector{UInt8})
    push!(CMPS, b); error("injected")
end
injected(e) = e isa ErrorException && e.msg == "injected"
for C in LEVELS
    id = C.IDENTIFIER; pk, sk = KEYS[C]
    ct, _ = C.kyber_kem_enc_derand(pk, fill(0x42, 32))
    reset(); threw = try C.kyber_kem_dec(ct, sk); false catch e; injected(e) end
    check("$id decaps, exception after re-encryption: propagated", threw)
    check("$id decaps, exception after re-encryption: cmp zeroed", length(CMPS) == 1 && zeroed(CMPS))
    check("$id decaps, exception after re-encryption: m′ ‖ H(ek) and (K′, r′) zeroed",
          length(GCALLS) == 1 && zeroed(GCALLS[1]))
    check("$id decaps, exception after re-encryption: polynomial temporaries zeroed", zeroed(ADDENDS))
end
reset(); threw = try XW.xwing_decaps(xct, xsk); false catch e; injected(e) end
check("X-Wing decaps, exception inside ML-KEM: propagated", threw)
check("X-Wing decaps, exception inside ML-KEM: skM, skX zeroed",
      length(EXPANDED) == 1 && zeroed((EXPANDED[1].skM, EXPANDED[1].skX)) && x25519_clean())

const MS = Vector{Vector{UInt8}}()
function M.Category3.kyber_kem_enc_derand(pk::Vector{UInt8}, coins::Vector{UInt8})
    push!(MS, coins); error("injected")
end
reset(); empty!(MS)
threw = try XW.xwing_encaps_derand(xpk, collect(UInt8, 65:128)); false catch e; injected(e) end
check("X-Wing encaps_derand, exception inside ML-KEM: propagated", threw)
check("X-Wing encaps_derand, exception inside ML-KEM: m zeroed", length(MS) == 1 && zeroed(MS))
check("X-Wing encaps_derand, exception inside ML-KEM: ekX, ssX zeroed", length(XCALLS) == 2 && x25519_clean())
reset(); empty!(MS)
threw = try XW.xwing_encaps(xpk); false catch e; injected(e) end
check("X-Wing encaps, exception inside ML-KEM: eseed, m, ekX, ssX zeroed",
      threw && length(ESEEDS) == 1 && zeroed(ESEEDS) && zeroed(MS) && x25519_clean())
empty!(MS)
threw = try M.Category3.kyber_kem_enc(KEYS[M.Category3][1]); false catch e; injected(e) end
check("ML-KEM-768 kem_enc, exception in encapsulation: coins zeroed", threw && length(MS) == 1 && zeroed(MS))
println("DONE")
"""

@testset "ML-KEM / X-Wing regressions (sweep)" begin
    zb(v) = ZeroBasedBytes(collect(UInt8, v))
    mem(v) = (m = Memory{UInt8}(undef, length(v)); copyto!(m, v); m)
    LEVELS = (MLKEM.Category1, MLKEM.Category3, MLKEM.Category5)

    @testset "F7: no Base.vcat on byte data (runtime-dispatch gate)" begin
        src = joinpath(pkgdir(PQJulia), "src")
        files = [f for f in readdir(src) if (startswith(f, "kyber_") && endswith(f, ".jl")) ||
                                            f in ("xwing.jl", "x25519.jl")]
        @test length(files) == 4
        for f in files
            @test (f, occursin("vcat(", read(joinpath(src, f), String))) == (f, false)
        end
    end

    @testset "F5/F15: basemul scratch is caller-owned and zeroed" begin
        for C in LEVELS
            k = C.KYBER_K
            a = [Int16.(rand(-1664:1664, 256)) for _ in 1:k]
            b = [Int16.(rand(-1664:1664, 256)) for _ in 1:k]
            ref = MLKEM.kyber_poly_basemul_montgomery!(zeros(Int16, 256), a[1], b[1])
            for i in 2:k
                ref .+= MLKEM.kyber_poly_basemul_montgomery!(zeros(Int16, 256), a[i], b[i])
            end
            MLKEM.kyber_poly_reduce!(ref)
            tmp = fill(Int16(0x5a5a), 256); r = zeros(Int16, 256)
            C.kyber_polyvec_basemul_acc!(r, a, b, k, tmp)
            @test r == ref
            @test all(iszero, tmp)
        end
    end

    @testset "F5/F6/F15/F16: secrets zeroed on every exit (child process)" begin
        script = tempname() * ".jl"
        write(script, MLKEM_XWING_WIPE_CHILD)
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) $script`
        out = IOBuffer(); err = IOBuffer()
        proc = run(pipeline(ignorestatus(cmd); stdout=out, stderr=err))
        rm(script; force=true)
        lines = split(String(take!(out)), '\n')
        ok = success(proc) && "DONE" in lines
        ok || println(stderr, String(take!(err)))
        @test ok
        checks = [(String(p[3]), p[2] == "true") for p in (split(l, ' '; limit=3) for l in lines)
                  if p[1] == "CHECK"]
        @test length(checks) == 64
        for (name, passed) in checks
            @test (name, passed) == (name, true)
        end
    end

    @testset "F9: inputs with other axes give the same result as a Vector" begin
        for C in LEVELS
            pk, sk = C.kyber_kem_keypair_derand(collect(UInt8, 1:64))
            bad_first = copy(pk); bad_first[1] = 0xff; bad_first[2] |= 0x0f   # coefficient 0 = 4095 ≥ q
            bad_last = copy(pk); bad_last[C.KYBER_POLYVECBYTES] = 0xff        # last coefficient ≥ q
            for (ek, want) in ((pk, true), (bad_first, false), (bad_last, false))
                @test C.kyber_ek_check(ek) == want
                @test C.kyber_ek_check(zb(ek)) == want
                @test C.kyber_ek_check(mem(ek)) == want
            end
            bad_sk = copy(sk); bad_sk[end-40] ⊻= 0x01                         # corrupt H(ek)
            for (dk, want) in ((sk, true), (bad_sk, false))
                @test C.kyber_dk_check(zb(dk)) == want
                @test C.kyber_dk_check(mem(dk)) == want                       # MethodError before
                @test C.kyber_dk_check(view(repeat(dk, inner=2), 1:2:2length(dk))) == want
            end
            @test !C.kyber_ek_check(zb(pk[1:end-1])) && !C.kyber_dk_check(zb(sk[2:end]))
        end
        k = collect(UInt8, 1:32); u = collect(UInt8, 101:132)
        @test X25519.x25519(zb(k), zb(u)) == X25519.x25519(k, u) == X25519.x25519_ref(zb(k), zb(u))
        @test X25519.x25519(mem(k), view(u, :)) == X25519.x25519(k, u)
        @test X25519.x25519_base(zb(k)) == X25519.x25519_base(k)
        @test_throws ArgumentError X25519.x25519(zb(k[1:31]), u)

        seed = collect(UInt8, 1:32); eseed = collect(UInt8, 65:128)
        pk, sk = XWing.xwing_keypair_derand(seed)
        @test XWing.xwing_keypair_derand(zb(seed)) == (pk, sk)
        ct, ss = XWing.xwing_encaps_derand(pk, eseed)
        @test XWing.xwing_encaps_derand(zb(pk), zb(eseed)) == (ct, ss)
        @test XWing.xwing_decaps(zb(ct), zb(sk)) == ss
        @test XWing.xwing_decaps(zb(ct), sk) == XWing.xwing_decaps(ct, zb(sk)) == ss
        ct2, ss2 = XWing.xwing_encaps(zb(pk))
        @test XWing.xwing_decaps(ct2, sk) == ss2
        bad = copy(pk); bad[1] = 0xff; bad[2] = 0xff
        @test_throws ArgumentError XWing.xwing_encaps(zb(bad))
        @test_throws ArgumentError XWing.xwing_decaps(zb(ct[1:end-1]), sk)
    end

    @testset "F8: X-Wing accepts any AbstractVector eseed" begin
        pk, sk = XWing.xwing_keypair_derand(collect(UInt8, 1:32))
        e = collect(0x00:0x3f)
        ref = XWing.xwing_encaps_derand(pk, e)
        @test XWing.xwing_decaps(ref[1], sk) == ref[2]
        for es in (view(e, :), mem(e), 0x00:0x3f, zb(e), view(repeat(e, inner=2), 1:2:128),
                   reinterpret(UInt8, reinterpret(UInt16, e)))
            before = collect(es)
            @test XWing.xwing_encaps_derand(pk, es) == ref
            @test collect(es) == before                          # the caller's eseed is not wiped
        end
        @test XWing.xwing_encaps_derand(mem(pk), mem(e)) == ref
        @test XWing.xwing_decaps(mem(ref[1]), mem(sk)) == ref[2]
        @test_throws ArgumentError XWing.xwing_encaps_derand(pk, 0x00:0x3e)
        bad = copy(pk); bad[1] = 0xff; bad[2] = 0xff
        eb = copy(e)
        @test_throws ArgumentError XWing.xwing_encaps_derand(bad, eb)
        @test eb == e
    end
end
