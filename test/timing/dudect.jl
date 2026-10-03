# dudect timing test; see README, "Timing test".
using PQJulia, Random, Statistics

const N = parse(Int, get(ENV, "DUDECT_N", "100000"))
const CROPS = (1.0, 0.99, 0.95, 0.9, 0.75, 0.5)
const rng = Xoshiro(1)

welch(a, b) = (mean(a) - mean(b)) / sqrt(var(a) / length(a) + var(b) / length(b))

function dudect(name, f, inA, inB; n=N)
    for _ in 1:200; f(inA()); f(inB()); end
    cls = rand(rng, Bool, n)
    ins = [c ? inB() : inA() for c in cls]
    x = Vector{Float64}(undef, n)
    for i in 1:n
        v = ins[i]
        t0 = time_ns(); f(v); x[i] = time_ns() - t0
    end
    ts = map(CROPS) do p
        th = quantile(x, p)
        welch([x[i] for i in 1:n if !cls[i] && x[i] <= th], [x[i] for i in 1:n if cls[i] && x[i] <= th])
    end
    m = maximum(t -> isnan(t) ? 0.0 : abs(t), ts)    # NaN: crop below timer resolution
    println(rpad(name, 44), "max|t| = ", rpad(round(m, digits=2), 6), m > 4.5 ? " LEAK" : " ok",
            "   crops: ", join(round.(ts, digits=1), " "))
    m
end

# Fixed inputs are copied: a reused array stays in cache and times faster.
const K = MLKEM.Category3
const D = MLDSA.Category3
const Q = Int32(8380417)
const G2 = D.GAMMA2
const pk, sk = K.kyber_kem_keypair()
const ct = K.kyber_kem_enc(pk)[1]
const a = rand(rng, UInt8, K.KYBER_CIPHERTEXTBYTES)
const r = rand(rng, UInt8, 32)
const x = rand(rng, UInt8, 32)
const u = X25519.x25519_base(rand(rng, UInt8, 32))
const kfix = rand(RandomDevice(), UInt8, 32)
const M = (UInt64(1) << 51) - 1
const p = (M - 18, M, M, M, M)
const pm1 = (M - 19, M, M, M, M)
const m = rand(RandomDevice(), UInt8, 32)
const sigma = rand(RandomDevice(), UInt8, 32)
const s1 = rand(rng, Int32(-4):Int32(4), 256)
const dsk = D.dilithium_keygen()[2]
const dpool = [D.dilithium_keygen()[2] for _ in 1:32]

decomp(v) = (s = 0; for y in v; a1, a0 = D.decompose(y); s ⊻= a1 ⊻ a0; end; s)
hints(v) = (s = 0; for y in v; s += D.make_hint(y, Int32(1)); end; s)
noise(s) = (t = zeros(Int16, 256); MLKEM.kyber_poly_getnoise_eta1!(t, s, 0x00, 2); MLKEM.kyber_ntt!(t))

# Falcon: integer floating point by rounding path, SamplerZ by center, signing by key.
const FP = PQJulia.FNDSA.FalconFpr
const F5 = FNDSA.Falcon512
rf(lo, hi) = FP.fpr(rand(rng, (-1.0, 1.0)) * ldexp(1.0 + rand(rng), rand(rng, lo:hi)))
mant(lo, hi) = FP.fpr(ldexp(lo + (hi - lo) * rand(rng), rand(rng, -9:9)))
ops(f) = v -> (s = UInt64(0); for (y, z) in v; s ⊻= f(y, z).b; end; s)
pairs64(g) = () -> [g() for _ in 1:64]
const fska, fskb = F5.falcon_keygen()[2], F5.falcon_keygen()[2]
const eka, ekb = F5.falcon_expand_sk(fska), F5.falcon_expand_sk(fskb)
const fmsg = rand(rng, UInt8, 32)
fsign(ek) = FNDSA.Falcon.sign_poly(ek.sk, ek.gs, fmsg)
fsamp(μ) = FP.samplerz(μ, FP.fpr(1.5), FP.fpr(1.2778336969128337),
                       PQJulia.FNDSA.FalconChaCha.ChaCha20(rand(rng, UInt8, 56)))

worst = maximum([
    dudect("ML-KEM-768 decaps: valid vs tampered ct", c -> K.kyber_kem_dec(c, sk),
           () -> copy(ct), () -> (c = copy(ct); c[rand(rng, 1:length(c))] ⊻= 0x01; c)),
    dudect("ML-KEM-768 encaps: fixed vs random m", c -> K.kyber_kem_enc_derand(pk, c), () -> copy(m), () -> rand(rng, UInt8, 32)),
    dudect("ML-KEM noise + NTT: fixed vs random σ", noise, () -> copy(sigma), () -> rand(rng, UInt8, 32)),
    dudect("kyber_verify: equal vs first byte differs", b -> MLKEM.kyber_verify(a, b),
           () -> copy(a), () -> (b = copy(a); b[1] ⊻= 0xff; b)),
    dudect("kyber_cmov!: b = 0 vs b = 1", b -> MLKEM.kyber_cmov!(copy(r), x, b), () -> 0x00, () -> 0x01),
    dudect("ML-DSA-65 unpack_sk: fixed vs random key", D.unpack_sk, () -> copy(dsk), () -> copy(rand(rng, dpool))),
    dudect("ML-DSA ntt!: fixed vs random s1", MLDSA.ntt!, () -> copy(s1), () -> rand(rng, Int32(-4):Int32(4), 256)),
    dudect("ML-DSA decompose: wrap vs middle region", decomp,
           () -> rand(rng, Q - G2:Q - Int32(1), 256), () -> rand(rng, G2 + Int32(1):Q - G2 - Int32(1), 256)),
    dudect("ML-DSA make_hint: hint 0 vs hint 1", hints,
           () -> rand(rng, -G2:G2, 256), () -> rand(rng, G2 + Int32(1):Int32(4) * G2, 256)),
    dudect("X25519: scalar 0x00… vs 0xff…", k -> X25519.x25519(k, u), () -> zeros(UInt8, 32), () -> fill(0xff, 32)),
    dudect("X25519: fixed vs random scalar", k -> X25519.x25519(k, u), () -> copy(kfix), () -> rand(rng, UInt8, 32)),
    dudect("X25519 encoding: p vs p − 1", X25519.fe_tobytes, () -> p, () -> pm1),
    dudect("fpr add: aligned vs far exponents", ops(FP.add),
           pairs64(() -> (y = rf(-20, 20); (y, FP.fpr(Float64(y) * (1 + rand(rng)))))),
           pairs64(() -> (rf(10, 20), rf(-60, -50)))),
    dudect("fpr add: same vs opposite signs", ops(FP.add),
           pairs64(() -> (y = rf(-9, 9); (y, FP.fpr(copysign(abs(Float64(rf(-9, 9))), Float64(y)))))),
           pairs64(() -> (y = rf(-9, 9); (y, FP.fpr(-copysign(abs(Float64(rf(-9, 9))), Float64(y))))))),
    dudect("fpr mul: product mantissa < 2 vs ≥ 2", ops(FP.mul),
           pairs64(() -> (mant(1.0, 1.4), mant(1.0, 1.4))), pairs64(() -> (mant(1.5, 2.0), mant(1.5, 2.0)))),
    dudect("fpr mul: nonzero vs zero operand", ops(FP.mul),
           pairs64(() -> (rf(-20, 20), rf(-20, 20))), pairs64(() -> (rf(-20, 20), FP.ZERO))),
    dudect("fpr div: quotient mantissa < 1 vs ≥ 1", ops(FP.div),
           pairs64(() -> (y = rf(-5, 5); (y, FP.fpr(Float64(y) / (1 + rand(rng)) * 1.0000001)))),
           pairs64(() -> (y = rf(-5, 5); (FP.fpr(Float64(y) * (1 + rand(rng))), y)))),
    dudect("fpr sqrt: even vs odd exponent", v -> (s = UInt64(0); for y in v; s ⊻= sqrt(y).b; end; s),
           pairs64(() -> FP.fpr(ldexp(1.0 + rand(rng), 2rand(rng, -10:10)))),
           pairs64(() -> FP.fpr(ldexp(1.0 + rand(rng), 2rand(rng, -10:10) + 1)))),
    dudect("Falcon SamplerZ: center frac < 0.1 vs ≈ 0.5", fsamp,
           () -> FP.fpr(rand(rng, -50:50) + 0.1rand(rng)), () -> FP.fpr(rand(rng, -50:50) + 0.45 + 0.1rand(rng))),
    dudect("Falcon-512 expand_sk: key A vs key B", F5.falcon_expand_sk, () -> copy(fska), () -> copy(fskb); n = N ÷ 10),
    dudect("Falcon-512 sign: key A vs key B", fsign, () -> eka, () -> ekb; n = N ÷ 10),
])
println(worst > 4.5 ? "\nTIMING DIFFERENCE DETECTED" : "\nNO TIMING DIFFERENCE DETECTED", " (N = $N, $(strip(Sys.cpu_info()[1].model)))")
