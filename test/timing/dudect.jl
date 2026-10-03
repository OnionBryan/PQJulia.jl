# dudect timing test; see README, "Timing test".
using PQJulia, Random, Statistics

const N = parse(Int, get(ENV, "DUDECT_N", "100000"))
const CROPS = (1.0, 0.99, 0.95, 0.9, 0.75, 0.5)
const rng = Xoshiro(1)

welch(a, b) = (mean(a) - mean(b)) / sqrt(var(a) / length(a) + var(b) / length(b))

function dudect(name, f, inA, inB)
    for _ in 1:200; f(inA()); f(inB()); end
    cls = rand(rng, Bool, N)
    ins = [c ? inB() : inA() for c in cls]
    x = Vector{Float64}(undef, N)
    for i in 1:N
        v = ins[i]
        t0 = time_ns(); f(v); x[i] = time_ns() - t0
    end
    ts = map(CROPS) do p
        th = quantile(x, p)
        welch([x[i] for i in 1:N if !cls[i] && x[i] <= th], [x[i] for i in 1:N if cls[i] && x[i] <= th])
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
])
println(worst > 4.5 ? "\nTIMING DIFFERENCE DETECTED" : "\nNO TIMING DIFFERENCE DETECTED", " (N = $N, $(strip(Sys.cpu_info()[1].model)))")
