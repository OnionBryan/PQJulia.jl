# Defects found by the adversarial review of the sweep fixes: a Falcon secret key whose reads fail
# part-way (a key store, an offset-axes or 1-based non-Vector key) leaves no copy and no decoded
# field behind, and Shamir reconstruction leaves no plain-data copy of Int shares. Observed through
# recording methods (forwarding with `invoke`) in a child process.
@testset "Review regressions: failed key reads, Int shares (child process)" begin
    child = raw"""
    using PQJulia
    const FAL = FNDSA.Falcon; const F = FNDSA.Falcon512
    check(name, ok) = println("CHECK ", name, " ", ok)
    # A key store whose read of byte `lim + 1` fails; axes start at off + 1.
    struct FailingKey <: AbstractVector{UInt8}; a::Vector{UInt8}; lim::Int; off::Int; end
    Base.size(k::FailingKey) = size(k.a)
    Base.axes(k::FailingKey) = k.off == 0 ? (Base.OneTo(length(k.a)),) : (Base.IdentityUnitRange(k.off+1:k.off+length(k.a)),)
    Base.getindex(k::FailingKey, i::Int) = i - k.off > k.lim ? error("key store read failed") : k.a[i - k.off]
    const SPY = Ref(false); const DEST = Vector{UInt8}[]; const DEC = Any[]; const SUB = Any[]
    Base.copyto!(dest::Vector{UInt8}, src::FailingKey) =
        (SPY[] && push!(DEST, dest); invoke(copyto!, Tuple{AbstractArray, AbstractArray}, dest, src))
    FAL.decode_signed(bytes::SubArray{UInt8, 1, FailingKey}, n::Int, bits::Int) =
        (r = invoke(FAL.decode_signed, NTuple{3, Any}, bytes, n, bits); SPY[] && push!(DEC, r); r)
    # Decoding from a plain Vector (the wrappers' copy): records f and g, and with INJECT set throws
    # when F (8-bit fields) is reached, after f and g are decoded.
    const INJECT = Ref(false); struct Injected <: Exception end
    FAL.decode_signed(bytes::SubArray{UInt8, 1, Vector{UInt8}, Tuple{UnitRange{Int}}, true}, n::Int, bits::Int) =
        (INJECT[] && bits == 8 && throw(Injected());
         r = invoke(FAL.decode_signed, NTuple{3, Any}, bytes, n, bits); SPY[] && push!(DEC, r); r)
    Base.getindex(v::Vector{Tuple{Int, Int}}, r::UnitRange{Int}) =
        (c = invoke(getindex, Tuple{Array, AbstractUnitRange}, v, r); SPY[] && push!(SUB, c); c)
    zero_or_nothing(x) = x === nothing || all(iszero, x)
    function spied(f)
        foreach(empty!, (DEST, DEC, SUB)); SPY[] = true
        try f() catch e; e finally SPY[] = false end
    end
    function main()
        pk, sk = F.falcon_keygen()
        calls = ("sign" => k -> F.falcon_sign(UInt8[1, 2, 3], k),
                 "mrm_sign" => k -> F.falcon_mrm_sign(zeros(Int, F.MRM_M1_LEN), UInt8[], k),
                 "expand_sk" => k -> F.falcon_expand_sk(k),
                 "certify" => k -> F.falcon_certify(k))
        flen = 1 + 2 * (512 * 6 ÷ 8)                   # header, f and g: the read fails inside F
        for (off, lim, tag) in ((10, length(sk) - 1, "offset"), (0, flen + 10, "onebased"))
            for (name, call) in calls
                e = spied(() -> call(FailingKey(sk, lim, off)))
                check("$tag-$name-threw", e isa ErrorException)
                check("$tag-$name-no-copy", !isempty(DEST) && all(d -> all(iszero, d), DEST))
            end
        end
        # The copy succeeds and decoding fails after f and g: the copy, f and g are wiped.
        INJECT[] = true
        for (name, call) in calls
            e = spied(() -> call(FailingKey(sk, length(sk), 10)))
            check("inject-$name-threw", e isa Injected)
            check("inject-$name-no-copy", length(DEST) == 1 && all(iszero, DEST[1]))
            check("inject-$name-no-field", length(DEC) == 2 && all(d -> d !== nothing && all(iszero, d), DEC))
        end
        e = spied(() -> FAL.decode_sk(copy(sk)))
        check("inject-decode_sk-no-field", e isa Injected && length(DEC) == 2 && all(d -> all(iszero, d), DEC))
        INJECT[] = false
        # The math-level decoder on a 1-based key store: f and g are decoded, then reading F fails.
        e = spied(() -> FAL.decode_sk(FailingKey(sk, flen + 10, 0)))
        check("decode_sk-threw", e isa ErrorException && length(DEC) == 2)
        check("decode_sk-fields-zero", all(zero_or_nothing, DEC))
        # Int shares: reconstruction makes no copy of the k shares it uses.
        p = 9223372036854775783                        # the largest prime below 2^63
        secret = 1234567890123
        ish = Tuple{Int, Int}[(x, Int(y)) for (x, y) in shamir_share(secret, 3, 5, p)]
        r = spied(() -> shamir_reconstruct(ish, 3, p))
        check("shamir-int-correct", r == secret)
        check("shamir-int-no-copy", all(c -> all(t -> t == (0, 0), c), SUB))
        println("DONE")
    end
    main()
    """
    dir = mktempdir()
    try
        script = joinpath(dir, "review.jl"); write(script, child)
        proj = Base.active_project()
        cmd = proj === nothing ? `$(Base.julia_cmd()) --startup-file=no $script` :
                                 `$(Base.julia_cmd()) --startup-file=no --project=$proj $script`
        err = IOBuffer()
        out = read(pipeline(ignorestatus(cmd); stderr = err), String)
        occursin("DONE", out) || println(stderr, "review child output:\n", out, "\n", String(take!(err)))
        got = Dict(String(t[2]) => t[3] == "true" for t in (split(l, ' ') for l in split(out, '\n')) if length(t) == 3 && t[1] == "CHECK")
        @test occursin("DONE", out)
        names = String[]
        for tag in ("offset", "onebased"), name in ("sign", "mrm_sign", "expand_sk", "certify"), what in ("threw", "no-copy")
            push!(names, "$tag-$name-$what")
        end
        for name in ("sign", "mrm_sign", "expand_sk", "certify"), what in ("threw", "no-copy", "no-field")
            push!(names, "inject-$name-$what")
        end
        push!(names, "inject-decode_sk-no-field")
        append!(names, ["decode_sk-threw", "decode_sk-fields-zero", "shamir-int-correct", "shamir-int-no-copy"])
        for name in names
            @test (name => get(got, name, false)) == (name => true)
        end
    finally
        rm(dir; recursive = true, force = true)
    end
end
