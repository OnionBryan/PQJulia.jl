# Widened Vararg specializations and dynamic calls on secret paths. A method `f(xs...)` is
# compiled for a widened signature (`Vararg{T}` of unknown length) unless it is written
# `xs::Vararg{T, N} where {N}`, and an argument of function type that a method only passes on
# (`randombytes`) is not specialized unless it is typed `::R where {R}`; in both cases the code
# that runs calls through dynamic dispatch. JET's report_opt analyses the concrete call and does
# not report either, so test/timing/dispatch.jl cannot see them. Here the optimized code of every
# path in that gate is walked through each PQJulia method it invokes (and each Base method invoked
# with a PQJulia function or closure as an argument: map, foreach, ntuple). The walk fails on
#   - an invoked PQJulia specialization whose signature ends in a Vararg,
#   - a dynamic call (`:call`, `Core.kwcall`, `Core._apply_iterate`) of a PQJulia function,
#   - a `Core._apply_iterate` splat in a PQJulia method,
# and each path must reach the function that does its work, so a walk that stops early fails.
using InteractiveUtils

module SpecializationCheck
using PQJulia

inpkg(m::Module) = m === PQJulia || (parentmodule(m) !== m && inpkg(parentmodule(m)))

# The value a statement argument stands for, when it is known at compile time.
function constval(ci, x)
    x isa GlobalRef && return isdefined(x.mod, x.name) ? getglobal(x.mod, x.name) : nothing
    x isa QuoteNode && return x.value
    if x isa Core.SSAValue
        t = ci.ssavaluetypes[x.id]
        return t isa Core.Const ? t.val : nothing
    end
    x isa Union{Core.Argument, Core.SlotNumber, Expr} && return nothing
    x
end

# `own(m)`: whether module m holds code under test (the package; the probes below add their own).
function walk(f, types; maxdepth = 16, own = inpkg)
    ownf(g) = own(parentmodule(typeof(g)))
    ownfuntype(T) = T isa DataType && T <: Function && own(T.name.module)
    bad = Set{String}(); seen = Set{Any}(); reached = Set{Symbol}()
    function visit(sig, depth, ownpkg)
        (sig in seen || depth > maxdepth) && return
        push!(seen, sig)
        cis = try Base.code_typed_by_type(sig) catch; return end
        for (ci, _) in cis
            ci isa Core.CodeInfo || continue
            for st in ci.code
                if Meta.isexpr(st, :invoke)
                    mi = st.args[1]; mi isa Core.CodeInstance && (mi = mi.def)
                    mi isa Core.MethodInstance || continue
                    m = mi.def; m isa Method || continue
                    ps = Base.unwrap_unionall(mi.specTypes).parameters
                    if own(m.module)
                        push!(reached, m.name)
                        !isempty(ps) && Base.isvarargtype(ps[end]) && push!(bad, "widened: " * string(mi.specTypes))
                        visit(mi.specTypes, depth + 1, true)
                    elseif any(ownfuntype, ps)
                        visit(mi.specTypes, depth + 1, false)       # map/foreach/ntuple over a package function
                    end
                elseif Meta.isexpr(st, :call)
                    g = constval(ci, st.args[1])
                    if g === Core._apply_iterate || g === Core.kwcall
                        ownpkg && g === Core._apply_iterate && push!(bad, "splat in $(sig)")
                        for x in st.args[2:end]
                            h = constval(ci, x)
                            h isa Function && ownf(h) && push!(bad, "dynamic call of $h in $(sig)")
                        end
                    elseif g isa Function && !(g isa Core.Builtin) && ownf(g)
                        push!(bad, "dynamic call of $g in $(sig)")
                    elseif g === nothing && ownpkg
                        push!(bad, "dynamic call through $(st.args[1]) in $(sig)")   # callee unknown until run time
                    end
                end
            end
        end
    end
    visit(Base.signature_type(f, types), 0, true)
    sort!(collect(bad)), reached
end

# (name, function, argument types, a function the walk must reach)
function paths()
    V = Vector{UInt8}
    ps = Any[]
    for (lv, C) in ((512, MLKEM.Category1), (768, MLKEM.Category3), (1024, MLKEM.Category5))
        push!(ps, ("ML-KEM-$lv keygen", C.kyber_kem_keypair_derand, (V,), :kyber_indcpa_keypair_derand),
                  ("ML-KEM-$lv encaps", C.kyber_kem_enc_derand, (V, V), :kyber_indcpa_enc),
                  ("ML-KEM-$lv decaps", C.kyber_kem_dec, (V, V), :kyber_indcpa_dec))
    end
    for (lv, C) in ((44, MLDSA.Category2), (65, MLDSA.Category3), (87, MLDSA.Category5))
        push!(ps, ("ML-DSA-$lv keygen", C.dilithium_keygen_derand, (V,), :poly_uniform_eta!),
                  ("ML-DSA-$lv sign", C.dilithium_sign_derand, (V, V, V), :compute_challenge!),
                  ("ML-DSA-$lv sign (hedged)", C.dilithium_sign, (V, V), :compute_challenge!))
    end
    for F in (FNDSA.Falcon512, FNDSA.Falcon1024)
        EK = typeof(F.falcon_expand_sk(F.falcon_keygen()[2]))
        push!(ps, ("$(F.IDENTIFIER) expand_sk", F.falcon_expand_sk, (V,), :ffldl),
                  ("$(F.IDENTIFIER) sign", F.falcon_sign, (V, EK), :samplerz),
                  ("$(F.IDENTIFIER) sign (byte key)", F.falcon_sign, (V, V), :samplerz),
                  ("$(F.IDENTIFIER) MRM sign", F.falcon_mrm_sign, (Vector{Int}, V, EK), :samplerz),
                  ("$(F.IDENTIFIER) keygen", F.falcon_keygen, (), :keygen_candidate))
    end
    push!(ps, ("X25519", X25519.x25519, (V, V), :fe_mul),
              ("X-Wing keygen", XWing.xwing_keypair_derand, (V,), :kyber_indcpa_keypair_derand),
              ("X-Wing decaps", XWing.xwing_decaps, (V, V), :kyber_indcpa_dec),
              ("X-Wing encaps", XWing.xwing_encaps_derand, (V, V), :kyber_indcpa_enc),
              ("Shamir share bytes", shamir_share_bytes, (V, Int, Int), :shamir_share),
              ("Shamir reconstruct bytes", shamir_reconstruct_bytes, (Vector{Tuple{Int, BigInt}}, Int, Int), :shamir_reconstruct))
    ps
end
end # module SpecializationCheck

@testset "No widened Vararg specializations or dynamic calls on secret paths" begin
    for (name, f, T, must) in SpecializationCheck.paths()
        bad, reached = SpecializationCheck.walk(f, T)
        @test (name => bad) == (name => String[])
        @test (name => must in reached) == (name => true)
    end
    # wipe! and catbytes are straight-line code at any argument count (no tuple splat).
    for (g, T) in ((PQJulia.Wipe.wipe!, ntuple(i -> isodd(i) ? Vector{UInt8} : Vector{Int32}, 40)),
                   (PQJulia.Bytes.catbytes, ntuple(_ -> Vector{UInt8}, 40)))
        ci, _ = only(Base.code_typed_by_type(Tuple{typeof(g), T...}))
        @test !any(st -> Meta.isexpr(st, :call) && occursin("_apply_iterate", string(st)), ci.code)
    end
    # The walk detects a widened helper and a dynamic call, so its silence above means something.
    @eval SpecializationCheck begin
        @noinline widecat(xs::Union{UInt8, AbstractVector{UInt8}}...) = PQJulia.Bytes.catbytes(xs...)   # as when not inlined
        probe_widened(a, b) = widecat(a, 0x01, b)
        const HOLDER = Ref{Any}(PQJulia.Bytes.catbytes)
        probe_dynamic(r) = HOLDER[](r, r)                     # callee known only at run time
    end
    V = Vector{UInt8}
    probe(m) = SpecializationCheck.inpkg(m) || m === SpecializationCheck
    @test !isempty(SpecializationCheck.walk(SpecializationCheck.probe_widened, (V, V); own = probe)[1])
    @test !isempty(SpecializationCheck.walk(SpecializationCheck.probe_dynamic, (V,); own = probe)[1])
    # catbytes allocates only its result, as vcat does.
    a = rand(UInt8, 32); b = rand(UInt8, 64); c = rand(UInt8, 16)
    cat3(a, b, c) = PQJulia.Bytes.catbytes(a, b, c)
    cat3(a, b, c)
    @test @allocated(cat3(a, b, c)) == @allocated(vcat(a, b, c))
end
