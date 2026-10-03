# Runtime-dispatch check; see README, "Timing test".
using Pkg
Pkg.activate(temp=true; io=devnull)
Pkg.develop(path=dirname(dirname(@__DIR__)); io=devnull)
Pkg.add("JET"; io=devnull)
using PQJulia, JET

C = MLKEM.Category3
paths = [
    ("ML-KEM-768 decaps", C.kyber_kem_dec, (Vector{UInt8}, Vector{UInt8})),
    ("ML-KEM-768 encaps", C.kyber_kem_enc, (Vector{UInt8},)),
    ("X25519", X25519.x25519, (Vector{UInt8}, Vector{UInt8})),
    ("X-Wing decaps", XWing.xwing_decaps, (Vector{UInt8}, Vector{UInt8})),
    ("X-Wing encaps", XWing.xwing_encaps, (Vector{UInt8},)),
]
bad = 0
for (name, f, T) in paths
    n = length(JET.get_reports(report_opt(f, T)))
    global bad += n
    println(rpad(name, 22), n, " runtime-dispatch reports")
end
println(bad == 0 ? "\nNO RUNTIME DISPATCH" : "\nRUNTIME DISPATCH FOUND")
exit(bad == 0 ? 0 : 1)
