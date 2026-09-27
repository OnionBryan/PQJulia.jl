#!/usr/bin/env julia
#= test/falcon_ffsampling.jl — ffSampling signer and fast NTRU solver vs their references.
   Karatsuba negamul = schoolbook; ffSampling signatures verify, satisfy s0 = c − s1·h,
   reject forgeries, and match the Klein/GPV norm distribution; n=512 timing. Run:
       julia test/falcon_ffsampling.jl =#
using Statistics, Printf
include(joinpath(@__DIR__, "..", "src", "falcon", "falcon.jl")); using .Falcon
const F = Falcon
import .FalconNTRUGen as NG

pass = 0; fail = 0
chk(c, m) = (global pass, fail; c ? (pass += 1; println("  ok   ", m)) : (fail += 1; println("  FAIL ", m)))

println("── Karatsuba negamul = schoolbook (BigInt, 4..3000-bit coefficients) ──")
let bad = 0
    for n in (1, 2, 4, 16, 32, 64, 128, 512), _ in 1:6
        bits = rand((4, 60, 300, 3000))
        a = [rand(-big(2)^bits:big(2)^bits) for _ in 1:n]; b = [rand(-big(2)^bits:big(2)^bits) for _ in 1:n]
        NG.negamul(a, b) == NG.negamul_school(a, b) || (bad += 1)
    end
    chk(bad == 0, "negamul == negamul_school on 48 random cases")
end

println("── ffSampling vs Klein/GPV reference ──")
for n in (16, 64)
    sk = F.keygen(n); pk = sk.h
    gs = F.sign_setup(sk); gk = F.sign_setup_klein(sk)
    nf = Float64[]; nk = Float64[]; bad = 0
    for i in 1:400
        msg = Vector{UInt8}("m$i")
        s = F.sign_poly(sk, gs, msg)
        F.verify_poly(pk, n, msg, s) || (bad += 1)
        c = F.hash_to_point(msg, s.salt, n)
        s.s0 == F.centermod.(c .- F.poly_mul_modq(s.s1, pk)) || (bad += 1)
        F.verify_poly(pk, n, Vector{UInt8}("x$i"), s) && (bad += 1)
        push!(nf, s.norm2)
        push!(nk, F.sign_poly_klein(sk, gk, msg).norm2)
    end
    z = abs(mean(nf) - mean(nk)) / sqrt(var(nf) / 400 + var(nk) / 400)
    chk(bad == 0, @sprintf("n=%d: 400 signatures verify, s0 = c − s1·h, forgeries rejected", n))
    chk(z < 4, @sprintf("n=%d: mean‖s‖² ffSampling %.0f vs Klein %.0f (%.1f SE)", n, mean(nf), mean(nk), z))
end

println("── Falcon-512 speed ──")
let sk = F.keygen(512)
    F.sign_setup(sk)
    gs = F.sign_setup(sk)
    F.sign_poly(sk, gs, UInt8[1])
    t = median([(@elapsed F.sign_poly(sk, gs, Vector{UInt8}("p$i"))) for i in 1:10])
    ok = all(F.verify_poly(sk.h, 512, Vector{UInt8}("q$i"), F.sign_poly(sk, gs, Vector{UInt8}("q$i"))) for i in 1:20)
    chk(ok, "n=512: 20/20 signatures verify")
    chk(t < 0.25, @sprintf("n=512: sign median %.4fs", t))
    tk = median([(@elapsed F.keygen(512)) for _ in 1:3])
    chk(tk < 15, @sprintf("n=512: keygen median %.2fs (under the 15s client timeout)", tk))
end

println()
println(fail == 0 ? "ALL $pass FFSAMPLING CHECKS PASS" : "$fail FAILED, $pass passed")
exit(fail == 0 ? 0 : 1)
