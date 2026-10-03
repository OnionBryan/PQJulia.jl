# dudect timing test; see README, "Timing test".
using PQJulia, Random, Statistics

const N = parse(Int, get(ENV, "DUDECT_N", "100000"))
const CROPS = (1.0, 0.99, 0.95, 0.9, 0.75, 0.5)
rng = Xoshiro(1)

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

K = MLKEM.Category3
pk, sk = K.kyber_kem_keypair()
ct, _ = K.kyber_kem_enc(pk)
a = rand(rng, UInt8, K.KYBER_CIPHERTEXTBYTES)
r = rand(rng, UInt8, 32); x = rand(rng, UInt8, 32)
u = X25519.x25519_base(rand(rng, UInt8, 32))
kfix = rand(RandomDevice(), UInt8, 32)
M = (UInt64(1) << 51) - 1
p = (M - 18, M, M, M, M); pm1 = (M - 19, M, M, M, M)

worst = maximum([
    dudect("ML-KEM-768 decaps: valid vs tampered ct", c -> K.kyber_kem_dec(c, sk),
           () -> copy(ct), () -> (c = copy(ct); c[rand(rng, 1:length(c))] ⊻= 0x01; c)),
    dudect("kyber_verify: equal vs first byte differs", b -> MLKEM.kyber_verify(a, b),
           () -> copy(a), () -> (b = copy(a); b[1] ⊻= 0xff; b)),
    dudect("kyber_cmov!: b = 0 vs b = 1", b -> MLKEM.kyber_cmov!(copy(r), x, b), () -> 0x00, () -> 0x01),
    dudect("X25519: scalar 0x00… vs 0xff…", k -> X25519.x25519(k, u), () -> zeros(UInt8, 32), () -> fill(0xff, 32)),
    dudect("X25519: fixed vs random scalar", k -> X25519.x25519(k, u), () -> copy(kfix), () -> rand(rng, UInt8, 32)),
    dudect("X25519 encoding: p vs p − 1", X25519.fe_tobytes, () -> p, () -> pm1),
])
println(worst > 4.5 ? "\nTIMING DIFFERENCE DETECTED" : "\nNO TIMING DIFFERENCE DETECTED", " (N = $N, $(strip(Sys.cpu_info()[1].model)))")
