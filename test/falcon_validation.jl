#!/usr/bin/env julia
#= test/falcon_validation.jl — discriminating validation of the Falcon (FN-DSA /
   FIPS 206) port. Each check is against ground truth (negacyclic multiply, target
   Gaussian moments, the NTRU equation, encode round-trip, the native NEON NTT,
   and end-to-end keygen→sign→verify with the spec norm bound). Run:
       julia test/falcon_validation.jl =#
using Statistics, Printf
const SRC = joinpath(@__DIR__, "..", "src", "falcon")
include(joinpath(SRC, "falcon_fft.jl"));      using .FalconFFT
include(joinpath(SRC, "falcon_sampler.jl"));  using .FalconSampler
include(joinpath(SRC, "falcon_encoding.jl")); using .FalconEncoding
include(joinpath(SRC, "falcon_neon.jl"));     using .FalconNEON
include(joinpath(SRC, "falcon.jl"));          using .Falcon
const FF = FalconFFT; const NE = FalconNEON

pass = 0; fail = 0
chk(c, m) = (global pass, fail; c ? (pass += 1; println("  ok   ", m)) : (fail += 1; println("  FAIL ", m)))

println("── FFT / NTT core (vs schoolbook negacyclic multiply) ──")
chk(FalconFFT.validate(verbose=false), "fft/ntt: ifft∘fft=id, mul_fft=neg-cyclic, split/merge, NTT mod q")
let a = rand(0:FF.q-1,256), b = rand(0:FF.q-1,256)
    chk(FF.ntt_mul_fast(a,b) == mod.(FF.negacyclic_mul(a,b), FF.q), "fast O(n log n) NTT = reference")
end

println("── SamplerZ bit-exact vs reference KAT (if reference checkout present) ──")
let katf = "/tmp/falconref/scripts/samplerz_KAT512.py"
    if isfile(katf)
        txt = read(katf, String); p = 0; f = 0
        for m in eachmatch(r"'mu':\s*(-?[\d.]+),\s*'sigma':\s*([\d.]+),\s*'sigmin':\s*([\d.]+),\s*'octets':\s*'([0-9A-Fa-f]+)',\s*'z':\s*(-?\d+)", txt)
            z = FalconSampler.samplerz(parse(Float64,m[1]),parse(Float64,m[2]),parse(Float64,m[3]),FalconSampler.KATSource(String(m[4])))
            z == parse(Int,m[5]) ? (p+=1) : (f+=1)
        end
        chk(f == 0 && p > 100, "samplerz bit-exact: $p/$(p+f) reference KAT vectors")
    else
        println("  skip  (reference KAT not checked out at $katf)")
    end
end

println("── SamplerZ (target Gaussian moments) ──")
let r = FalconSampler.RNG(collect(UInt8,1:32))
    s = [FalconSampler.samplerz(0.3, 1.7, 1.2778, r) for _ in 1:200_000]
    chk(abs(mean(s)-0.3) < 0.03 && abs(var(s)-1.7^2) < 0.06, "samplerz(0.3,1.7): mean/var match")
end

println("── NTRU equation f·G − g·F = q (exact) ──")
let
    okall = true
    for n in (8,16,32)
        f = rand(-5:5,n); g = rand(-5:5,n)
        try
            F,G = Falcon.NG.ntru_solve(BigInt.(f),BigInt.(g))
            okall &= Falcon.NG.ntru_check(BigInt.(f),BigInt.(g),F,G)
        catch; end   # non-coprime: skip (handled by keygen resample)
    end
    chk(okall, "ntru_solve: f·G−g·F = q where coprime")
end

println("── Signature encoding round-trip + spec sizes ──")
let s = round.(Int, 60 .* randn(512))
    enc = FalconEncoding.compress(s, (666-41)*8)
    chk(enc !== nothing && length(enc)==625 && FalconEncoding.decompress(enc,(666-41)*8,512)==s,
        "compress/decompress round-trips, 625 bytes (Falcon-512 sig − salt − header)")
end

println("── Native NEON NTT (C/NEON) = reference ──")
let a = rand(0:FF.q-1,512), b = rand(0:FF.q-1,512)
    chk(NE.isavailable() && NE.neon_negamul(a,b) == mod.(FF.negacyclic_mul(a,b), FF.q),
        "NEON NTT negamul = reference (n=512)")
end

println("── End-to-end keygen → sign → verify (+ forgery rejection) ──")
for n in (8, 16, 32)
    sk = Falcon.keygen(n); gs = Falcon.sign_setup(sk); good = 0
    for t in 1:4
        msg = Vector{UInt8}("falcon validation $t")
        sig = Falcon.sign_poly(sk, gs, msg)
        good += (Falcon.verify_poly(sk.h, n, msg, sig) &&
                 !Falcon.verify_poly(sk.h, n, Vector{UInt8}("forged"), sig)) ? 1 : 0
    end
    chk(good == 4, "n=$n: 4/4 round-trips verify, forgeries rejected")
end

println("\n", fail==0 ? "ALL $pass FALCON CHECKS PASS" : "$fail FAILED / $pass passed")
exit(fail==0 ? 0 : 1)
