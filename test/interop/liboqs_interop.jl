# test/interop/liboqs_interop.jl — two-way interop with liboqs (independent C implementations of
# ML-KEM, ML-DSA and Falcon): keys, signatures and ciphertexts cross in both directions. Opt-in:
#   clang -O2 -dynamiclib -o test/interop/liboqs_shim.dylib test/interop/oqs_shim.c \
#     -I$(brew --prefix)/include -Wl,-force_load,$(brew --prefix)/lib/liboqs.a -L$(brew --prefix openssl@3)/lib -lcrypto
#   julia --project=. test/interop/liboqs_interop.jl
using PQJulia, Random
using Libdl
const LIB = joinpath(@__DIR__, "liboqs_shim." * Libdl.dlext)
rb(n) = rand(RandomDevice(), UInt8, n)

struct CSig; p::Ptr{Cvoid}; pk::Int; sk::Int; sig::Int; end
function CSig(name)
    p = ccall((:shim_sig_new, LIB), Ptr{Cvoid}, (Cstring,), name); p == C_NULL && error("no $name")
    CSig(p, [Int(ccall((:shim_sig_len, LIB), Csize_t, (Ptr{Cvoid}, Cint), p, w)) for w in 0:2]...)
end
function ckeypair(s::CSig)
    pk = zeros(UInt8, s.pk); sk = zeros(UInt8, s.sk)
    ccall((:shim_sig_keypair, LIB), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}), s.p, pk, sk) == 0 || error("keypair")
    pk, sk
end
function csign(s::CSig, m, sk; ctx=UInt8[])
    sig = zeros(UInt8, s.sig); len = Ref{Csize_t}(0)
    ccall((:shim_sig_sign, LIB), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Ref{Csize_t}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Ptr{UInt8}),
          s.p, sig, len, m, length(m), ctx, length(ctx), sk) == 0 || error("sign")
    sig[1:len[]]
end
cverify(s::CSig, m, sig, pk; ctx=UInt8[]) =
    ccall((:shim_sig_verify, LIB), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Ptr{UInt8}, Csize_t, Ptr{UInt8}),
          s.p, m, length(m), sig, length(sig), ctx, length(ctx), pk) == 0

flip(v) = (w = copy(v); i = rand(1:length(w)); w[i] ⊻= UInt8(1) << rand(0:7); w)

function sig_interop(name, jkeygen, jsign, jverify, trials; ctx=UInt8[])
    s = CSig(name); c = Dict{String,Int}()
    tick(k, ok) = (c[k] = get(c, k, 0) + ok)
    for t in 1:trials
        m = rb(rand(0:200))
        # C → Julia
        pk, sk = ckeypair(s); sig = csign(s, m, sk; ctx)
        tick("C sig → Julia verify", jverify(m, sig, pk, ctx))
        tick("C sig tampered → Julia rejects", !jverify(m, flip(sig), pk, ctx) && !jverify(flip(m == UInt8[] ? UInt8[0] : m), sig, pk, ctx))
        tick("C sk → Julia sign → C verify", cverify(s, m, jsign(m, sk, ctx), pk; ctx))
        # Julia → C
        jpk, jsk = jkeygen()
        jsig = jsign(m, jsk, ctx)
        tick("Julia sig → C verify", cverify(s, m, jsig, jpk; ctx))
        tick("Julia sig tampered → C rejects", !cverify(s, m, flip(jsig), jpk; ctx))
        tick("Julia sk → C sign → Julia verify", jverify(m, csign(s, m, jsk; ctx), jpk, ctx))
    end
    println(rpad(name * (isempty(ctx) ? "" : " (ctx)"), 26), join(["$k $v/$trials" for (k, v) in sort(collect(c))], " · "))
    all(==(trials), values(c))
end

function kem_interop(name, M, trials)
    p = ccall((:shim_kem_new, LIB), Ptr{Cvoid}, (Cstring,), name)
    L = [Int(ccall((:shim_kem_len, LIB), Csize_t, (Ptr{Cvoid}, Cint), p, w)) for w in 0:3]
    ok = zeros(Int, 4)
    for t in 1:trials
        pk = zeros(UInt8, L[1]); sk = zeros(UInt8, L[2])
        ccall((:shim_kem_keypair, LIB), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}), p, pk, sk)
        ct, ss = M.kyber_kem_enc(pk)                                   # Julia encaps to C key
        ss2 = zeros(UInt8, L[4])
        ccall((:shim_kem_decaps, LIB), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}), p, ss2, ct, sk)
        ok[1] += ss == ss2
        ok[2] += M.kyber_kem_dec(ct, sk) == ss                         # Julia decaps with C sk
        jpk, jsk = M.kyber_kem_keypair()                              # C encaps to Julia key
        ct2 = zeros(UInt8, L[3]); ss3 = zeros(UInt8, L[4])
        ccall((:shim_kem_encaps, LIB), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}), p, ct2, ss3, jpk)
        ok[3] += M.kyber_kem_dec(ct2, jsk) == ss3
        bad = flip(ct2); ss4 = zeros(UInt8, L[4])                      # implicit rejection agrees
        ccall((:shim_kem_decaps, LIB), Cint, (Ptr{Cvoid}, Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8}), p, ss4, bad, jsk)
        ok[4] += M.kyber_kem_dec(bad, jsk) == ss4 != ss3
    end
    println(rpad(name, 26), "Julia enc → C dec $(ok[1])/$trials · C sk → Julia dec $(ok[2])/$trials · C enc → Julia dec $(ok[3])/$trials · tampered ct: same rejection key $(ok[4])/$trials")
    all(==(trials), ok)
end

allok = true
for (name, M, n) in [("Falcon-padded-512", FNDSA.Falcon512, 25), ("Falcon-padded-1024", FNDSA.Falcon1024, 10)]
    global allok &= sig_interop(name, M.falcon_keygen, (m, sk, _) -> M.falcon_sign(m, sk),
                                (m, s, pk, _) -> M.falcon_verify(m, s, pk), n)
end
for (name, C) in [("ML-DSA-44", MLDSA.Category2), ("ML-DSA-65", MLDSA.Category3), ("ML-DSA-87", MLDSA.Category5)], ctx in (UInt8[], Vector{UInt8}("interop-ctx"))
    global allok &= sig_interop(name, C.dilithium_keygen, (m, sk, x) -> C.dilithium_sign(m, sk; context=x),
                                (m, s, pk, x) -> C.dilithium_verify(m, s, pk; context=x), 25; ctx)
end
for (name, M) in [("ML-KEM-512", MLKEM.Category1), ("ML-KEM-768", MLKEM.Category3), ("ML-KEM-1024", MLKEM.Category5)]
    global allok &= kem_interop(name, M, 50)
end
println(allok ? "\nALL INTEROP CHECKS PASS (liboqs $(unsafe_string(ccall((:OQS_version, LIB), Cstring, ()))))" : "\nINTEROP FAILURES")
