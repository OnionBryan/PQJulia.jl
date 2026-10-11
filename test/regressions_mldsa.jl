# ML-DSA regressions from the defect sweep:
#   F1          Base.vcat on byte strings: test/timing/dispatch.jl reported runtime dispatch on
#               every ML-DSA keygen and signing path (the assertion inside vcat, Julia 1.12.1)
#   F2          each signing attempt's w1 packing, H input, c̃, c, NTT(c) and SampleInBall stream
#               were never wiped
#   F3/F14/F26  a secret key with s1/s2 outside [-η, η], or any other exception while signing,
#               left the decoded K, s1, s2 and t0 unwiped
#   F4          dilithium_sign / dilithium_sign_prehash left the hedged rnd unwiped when they threw
# Included from runtests.jl, which has `using Test, PQJulia`.

@testset "ML-DSA regressions (sweep)" begin

@testset "F1: no Base.vcat in the ML-DSA sources" begin
    for f in ("mldsa.jl", "dilithium_core.jl", "dilithium_level.jl")
        hits = String[]
        for (i, l) in enumerate(eachline(joinpath(pkgdir(PQJulia), "src", f)))
            occursin(r"\bvcat\b", first(split(l, '#'; limit=2))) && push!(hits, "src/$f:$i: $(strip(l))")
        end
        @test isempty(hits)
    end
end

@testset "F1: hash inputs built by catbytes are the vcat ones" begin
    KC = PQJulia.Keccak
    rho = collect(UInt8, 101:132); seed = collect(UInt8, 1:64)
    # ExpandA: SHAKE128(ρ ‖ nonce_le16), nonce with both bytes set; rejection sampling of 23-bit words
    stream = KC.shake128(vcat(rho, UInt8[0x05, 0x03]), 168 * 20)
    ref = Int32[]; p = 1
    while length(ref) < 256
        t = (UInt32(stream[p]) | UInt32(stream[p+1]) << 8 | UInt32(stream[p+2]) << 16) & 0x7fffff
        p += 3
        t < MLDSA.Q && push!(ref, Int32(t))
    end
    @test MLDSA.poly_uniform!(zeros(Int32, 256), rho, 0x0305) == ref
    msg = collect(UInt8, 1:50); ctx = UInt8[7, 8, 9]; tr = collect(UInt8, 201:255)
    for C in (MLDSA.Category2, MLDSA.Category3, MLDSA.Category5)
        @testset "$(C.IDENTIFIER)" begin
            # ExpandMask / ExpandS: SHAKE256(ρ′ ‖ nonce_le16), nonces above 255
            @test C.poly_uniform_gamma1!(zeros(Int32, 256), seed, 0x1234) ==
                  C.polyz_unpack!(zeros(Int32, 256), KC.shake256(vcat(seed, UInt8[0x34, 0x12]), C.POLYZ_PACKED))
            buf = KC.shake256(vcat(seed, UInt8[0x02, 0x01]), 136 * 8); eta = zeros(Int32, 256)
            @test C.rej_eta!(eta, 0, 256, buf, length(buf)) == 256
            @test C.poly_uniform_eta!(zeros(Int32, 256), seed, 0x0102) == eta
            # M′ and μ
            @test C.pure_mprime(msg, ctx) == vcat(UInt8[0x00, 0x03], ctx, msg)
            @test C.prehash_mprime(msg, "SHA3-256", ctx) ==
                  vcat(UInt8[0x01, 0x03], ctx, C.HASH_OIDS["SHA3-256"], KC.sha3_256(msg))
            @test C.mu_of(tr, msg) == KC.shake256(vcat(tr, msg), 64)
        end
    end
end

# The wipe checks watch package internals through spies: methods more specific than the
# package's (or Base's) that forward through `invoke`, so behaviour is unchanged, and keep
# references to the buffers they see. Adding them would invalidate and change this process's
# method tables, so they run in a child process that reports one line per check.
child = raw"""
using PQJulia, Random
const W = PQJulia.Wipe
const KC = PQJulia.Keccak
const V = SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int}}, true}

const SPY = Ref(false)
const SHAKE = Vector{UInt8}[]     # SHAKE256 inputs and outputs
const COPIES = Vector{Int32}[]    # Int32 vectors copied, and their copies (y, NTT(y), c, NTT(c))
const POLYS = Vector{Int32}[]     # the buffers s1, s2 and t0 are decoded into
const WIPED = Vector{UInt8}[]     # what each byte buffer held when it was wiped
const RNDS = Vector{UInt8}[]      # bytes drawn from RandomDevice (keygen ξ, hedged rnd)
const CHAL_LEN = Ref(0)           # length of the challenge hash input μ ‖ w1Encode(w1)
const NCHAL = Ref(0)              # challenge hashes so far
const FAULT = Ref(0)              # > 0: the FAULT[]-th challenge hash throws Injected()
struct Injected <: Exception end

function KC.shake256(data::Vector{UInt8}, outlen::UInt64)
    if SPY[]
        push!(SHAKE, data)
        if length(data) == CHAL_LEN[]
            NCHAL[] += 1
            NCHAL[] == FAULT[] && throw(Injected())
        end
    end
    out = invoke(KC.shake256, Tuple{AbstractVector{UInt8}, Integer}, data, outlen)
    SPY[] && push!(SHAKE, out)
    out
end
function Base.copy(a::Vector{Int32})
    r = invoke(copy, Tuple{Array}, a)
    SPY[] && push!(COPIES, a, r)
    r
end
function W.wipe!(x::Vector{UInt8})
    SPY[] && push!(WIPED, invoke(copy, Tuple{Array}, x))
    invoke(W.wipe!, Tuple{Array}, x)
end
function Random.rand(rd::RandomDevice, ::Type{UInt8}, n::Int)
    r = invoke(rand, Tuple{Random.AbstractRNG, Type{UInt8}, Integer}, rd, UInt8, n)
    SPY[] && push!(RNDS, r)
    r
end
for C in (MLDSA.Category2, MLDSA.Category3, MLDSA.Category5)
    @eval begin
        $C.polyeta_unpack!(r::Vector{Int32}, a::V) =
            (invoke($C.polyeta_unpack!, Tuple{Vector{Int32}, AbstractVector{UInt8}}, r, a); SPY[] && push!(POLYS, r); r)
        $C.polyt0_unpack!(r::Vector{Int32}, a::V) =
            (invoke($C.polyt0_unpack!, Tuple{Vector{Int32}, AbstractVector{UInt8}}, r, a); SPY[] && push!(POLYS, r); r)
    end
end

# Run f with the spies on; return its value, or the exception it threw.
function spied(f)
    foreach(empty!, (SHAKE, COPIES, POLYS, WIPED, RNDS)); NCHAL[] = 0
    SPY[] = true
    try
        f()
    catch e
        e
    finally
        SPY[] = false
    end
end
allzero(bufs) = all(b -> all(iszero, b), bufs)
report(id, C, ok) = (println("CHECK\t", id, "\t", C.IDENTIFIER, "\t", ok); flush(stdout))

for C in (MLDSA.Category2, MLDSA.Category3, MLDSA.Category5)
    CHAL_LEN[] = C.CRHBYTES + C.K * C.POLYW1_PACKED
    pk, sk = C.dilithium_keygen_derand(fill(0x5a, 32))
    K = sk[33:64]
    msg = collect(UInt8, 1:40); rnd = fill(0x33, 32)
    npolys = C.L + 2 * C.K

    # F2. External-μ signing computes no public μ, so every SHAKE256 input and output of the
    # call (ρ″ and its input, each attempt's y stream, H input, c̃ and SampleInBall stream) is
    # secret until zeroed, and so is every Int32 vector it copies (y, NTT(y), c, NTT(c)).
    shake_ok = copies_ok = true; seen = 0; mu_multi = UInt8[]
    for m in 0x00:0x3f
        mu = fill(m, 64)
        sig = spied(() -> C.dilithium_sign_internal_mu(mu, sk, rnd))
        shake_ok &= sig isa Vector{UInt8} && C.dilithium_verify_mu(mu, sig, pk) && NCHAL[] >= 1 && allzero(SHAKE)
        copies_ok &= !isempty(COPIES) && allzero(COPIES)
        NCHAL[] >= 3 && isempty(mu_multi) && (mu_multi = mu)
        seen += 1
        seen >= 8 && !isempty(mu_multi) && break
    end
    report("F2.rejected", C, !isempty(mu_multi))
    report("F2.shake", C, shake_ok)
    report("F2.copies", C, copies_ok)

    # F3/F14/F26. A key whose last s2 byte decodes outside [-η, η], through every signing entry point.
    bad = copy(sk); bad[2 * 32 + 64 + (C.L + C.K) * C.POLYETA_PACKED] = 0xff
    entries = (s -> C.dilithium_sign(msg, s),
               s -> C.dilithium_sign_derand(msg, s, rnd),
               s -> C.dilithium_sign_internal(msg, s, rnd),
               s -> C.dilithium_sign_internal_mu(fill(0x01, 64), s, rnd),
               s -> C.dilithium_sign_prehash(msg, s, "SHA2-256"),
               s -> C.dilithium_sign_prehash_derand(msg, s, "SHA2-256", rnd))
    polys_ok = key_ok = true
    for f in entries
        e = spied(() -> f(bad))
        polys_ok &= e isa ArgumentError && length(POLYS) == npolys && allzero(POLYS)
        key_ok &= e isa ArgumentError && any(==(K), WIPED)
    end
    report("F3.polys", C, polys_ok)
    report("F3.key", C, key_ok)

    # F3. An exception in the middle of signing (thrown here by the second attempt's challenge
    # hash) still wipes everything the signer holds.
    FAULT[] = 2
    e = spied(() -> C.dilithium_sign_internal_mu(mu_multi, sk, rnd))
    FAULT[] = 0
    report("F3.fault", C, e isa Injected && NCHAL[] == 2 && length(POLYS) == npolys && allzero(POLYS) &&
                          allzero(SHAKE) && allzero(COPIES) && any(==(K), WIPED))

    # F4. The hedged rnd is wiped whether signing returns or throws.
    errs = (() -> C.dilithium_sign(msg, bad),
            () -> C.dilithium_sign(msg, sk[1:end-1]),
            () -> C.dilithium_sign(msg, sk; context=zeros(UInt8, 256)),
            () -> C.dilithium_sign_prehash(msg, bad, "SHA2-256"),
            () -> C.dilithium_sign_prehash(msg, sk, "MD5"),
            () -> C.dilithium_sign_prehash(msg, sk, "SHA2-256"; context=zeros(UInt8, 256)))
    ok = true
    for f in errs
        e = spied(f)
        ok &= e isa ArgumentError && length(RNDS) == 1 && all(iszero, RNDS[1])
    end
    report("F4.throw", C, ok)
    ok = true
    for f in (() -> C.dilithium_sign(msg, sk), () -> C.dilithium_sign_prehash(msg, sk, "SHA2-256"))
        sig = spied(f)
        ok &= sig isa Vector{UInt8} && length(RNDS) == 1 && all(iszero, RNDS[1])
    end
    report("F4.return", C, ok)
    kp = spied(() -> C.dilithium_keygen())
    report("keygen.xi", C, kp isa Tuple && length(RNDS) == 1 && all(iszero, RNDS[1]))
end
"""

checks = ("F2.rejected" => "a signature with rejected attempts was observed",
          "F2.shake"    => "every attempt's H input (μ ‖ w1Encode), c̃, SampleInBall stream, y stream and ρ″ zeroed",
          "F2.copies"   => "every attempt's c, NTT(c), y and NTT(y) zeroed",
          "F3.polys"    => "out-of-range key, 6 signing entry points: decoded s1, s2, t0 zeroed",
          "F3.key"      => "out-of-range key, 6 signing entry points: decoded K zeroed",
          "F3.fault"    => "exception in the second attempt: K, s1, s2, t0, ρ″, y and challenge buffers zeroed",
          "F4.throw"    => "hedged rnd zeroed when dilithium_sign / dilithium_sign_prehash throw (6 error paths)",
          "F4.return"   => "hedged rnd zeroed after a signature",
          "keygen.xi"   => "keygen seed ξ zeroed")

@testset "F2, F3, F4: secret buffers wiped on every attempt and every exit (child process)" begin
    proj = Base.active_project()
    projflag = proj === nothing ? String[] : ["--project=$proj"]
    cmd = addenv(`$(Base.julia_cmd()) --startup-file=no $projflag -e $child`,
                 "JULIA_LOAD_PATH" => join(LOAD_PATH, Sys.iswindows() ? ';' : ':'))
    out = IOBuffer(); err = IOBuffer()
    proc = run(pipeline(ignorestatus(cmd); stdout=out, stderr=err))
    got = Dict{Tuple{String, String}, Bool}()
    for l in eachline(seekstart(out))
        f = split(l, '\t')
        length(f) == 4 && f[1] == "CHECK" && (got[(f[2], f[3])] = f[4] == "true")
    end
    success(proc) || print(stderr, String(take!(err)))
    @test success(proc)
    for lv in ("ML-DSA-44", "ML-DSA-65", "ML-DSA-87"), (id, what) in checks
        @test (lv, what, get(got, (id, lv), false)) == (lv, what, true)
    end
end

end # ML-DSA regressions (sweep)
