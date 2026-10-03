# Runtime-dispatch check; see README, "Timing test".
using Pkg
Pkg.activate(temp=true; io=devnull)
Pkg.develop(path=dirname(dirname(@__DIR__)); io=devnull)
Pkg.add("JET"; io=devnull)
using PQJulia, JET

V = Vector{UInt8}
paths = Any[]
for (lv, C) in ((512, MLKEM.Category1), (768, MLKEM.Category3), (1024, MLKEM.Category5))
    push!(paths, ("ML-KEM-$lv keygen", C.kyber_kem_keypair_derand, (V,)),
                 ("ML-KEM-$lv encaps", C.kyber_kem_enc_derand, (V, V)),
                 ("ML-KEM-$lv decaps", C.kyber_kem_dec, (V, V)))
end
for (lv, C) in ((44, MLDSA.Category2), (65, MLDSA.Category3), (87, MLDSA.Category5))
    push!(paths, ("ML-DSA-$lv keygen", C.dilithium_keygen_derand, (V,)),
                 ("ML-DSA-$lv sign", C.dilithium_sign_derand, (V, V, V)),
                 ("ML-DSA-$lv unpack_sk", C.unpack_sk, (V,)))
end
for (lv, F) in ((512, FNDSA.Falcon512), (1024, FNDSA.Falcon1024))
    EK = typeof(F.falcon_expand_sk(F.falcon_keygen()[2]))
    push!(paths, ("Falcon-$lv expand_sk", F.falcon_expand_sk, (V,)),
                 ("Falcon-$lv sign", F.falcon_sign, (V, EK)),
                 ("Falcon-$lv MRM sign", F.falcon_mrm_sign, (Vector{Int}, V, EK)))
end
push!(paths, ("X25519", X25519.x25519, (V, V)),
             ("X-Wing decaps", XWing.xwing_decaps, (V, V)),
             ("X-Wing encaps", XWing.xwing_encaps, (V,)))

bad = 0
for (name, f, T) in paths
    n = length(JET.get_reports(report_opt(f, T)))
    global bad += n
    println(rpad(name, 22), n, " runtime-dispatch reports")
end
println(bad == 0 ? "\nNO RUNTIME DISPATCH" : "\nRUNTIME DISPATCH FOUND")
exit(bad == 0 ? 0 : 1)
