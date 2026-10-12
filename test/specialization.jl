# Widened Vararg specializations on secret paths. A method `f(xs...)` is compiled for a widened
# signature (`Vararg{T}` of unknown length) unless it is written `xs::Vararg{T, N} where {N}`; the
# caller then splats its arguments through runtime dispatch. JET's report_opt analyses the concrete
# call and does not report it, so test/timing/dispatch.jl cannot see this. Here the optimized code
# of every path in that gate is walked through each PQJulia method it invokes, and any invoked
# PQJulia specialization whose signature ends in a Vararg fails the test.
using InteractiveUtils

module SpecializationCheck
using PQJulia

inpkg(m::Module) = m === PQJulia || (parentmodule(m) !== m && inpkg(parentmodule(m)))

function widened(f, types; maxdepth = 12)
    bad = Set{String}(); seen = Set{Any}()
    function walk(sig, depth)
        (sig in seen || depth > maxdepth) && return
        push!(seen, sig)
        cis = try Base.code_typed_by_type(sig) catch; return end
        for (ci, _) in cis
            ci isa Core.CodeInfo || continue
            for st in ci.code
                Meta.isexpr(st, :invoke) || continue
                mi = st.args[1]; mi isa Core.CodeInstance && (mi = mi.def)
                mi isa Core.MethodInstance || continue
                m = mi.def; m isa Method && inpkg(m.module) || continue
                s = Base.unwrap_unionall(mi.specTypes)
                ps = s.parameters
                !isempty(ps) && Base.isvarargtype(ps[end]) && push!(bad, string(mi.specTypes))
                walk(mi.specTypes, depth + 1)
            end
        end
    end
    walk(Base.signature_type(f, types), 0)
    sort!(collect(bad))
end

function paths()
    V = Vector{UInt8}
    ps = Any[]
    for (lv, C) in ((512, MLKEM.Category1), (768, MLKEM.Category3), (1024, MLKEM.Category5))
        push!(ps, ("ML-KEM-$lv keygen", C.kyber_kem_keypair_derand, (V,)),
                  ("ML-KEM-$lv encaps", C.kyber_kem_enc_derand, (V, V)),
                  ("ML-KEM-$lv decaps", C.kyber_kem_dec, (V, V)))
    end
    for (lv, C) in ((44, MLDSA.Category2), (65, MLDSA.Category3), (87, MLDSA.Category5))
        push!(ps, ("ML-DSA-$lv keygen", C.dilithium_keygen_derand, (V,)),
                  ("ML-DSA-$lv sign", C.dilithium_sign_derand, (V, V, V)),
                  ("ML-DSA-$lv sign (hedged)", C.dilithium_sign, (V, V)))
    end
    for F in (FNDSA.Falcon512, FNDSA.Falcon1024)
        EK = typeof(F.falcon_expand_sk(F.falcon_keygen()[2]))
        push!(ps, ("$(F.IDENTIFIER) expand_sk", F.falcon_expand_sk, (V,)),
                  ("$(F.IDENTIFIER) sign", F.falcon_sign, (V, EK)),
                  ("$(F.IDENTIFIER) sign (byte key)", F.falcon_sign, (V, V)),
                  ("$(F.IDENTIFIER) MRM sign", F.falcon_mrm_sign, (Vector{Int}, V, EK)))
    end
    push!(ps, ("X25519", X25519.x25519, (V, V)),
              ("X-Wing keygen", XWing.xwing_keypair_derand, (V,)),
              ("X-Wing decaps", XWing.xwing_decaps, (V, V)),
              ("X-Wing encaps", XWing.xwing_encaps_derand, (V, V)),
              ("Shamir share bytes", shamir_share_bytes, (V, Int, Int)),
              ("Shamir reconstruct bytes", shamir_reconstruct_bytes, (Vector{Tuple{Int, BigInt}}, Int, Int)))
    ps
end
end # module SpecializationCheck

@testset "No widened Vararg specializations on secret paths" begin
    for (name, f, T) in SpecializationCheck.paths()
        @test (name => SpecializationCheck.widened(f, T)) == (name => String[])
    end
    # catbytes allocates only its result, as vcat does.
    a = rand(UInt8, 32); b = rand(UInt8, 64); c = rand(UInt8, 16)
    cat3(a, b, c) = PQJulia.Bytes.catbytes(a, b, c)
    cat3(a, b, c)
    @test @allocated(cat3(a, b, c)) == @allocated(vcat(a, b, c))
end
