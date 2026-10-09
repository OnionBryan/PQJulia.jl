# Regression tests for three classes of implementation fault reported in 2026, none of which an
# honest round trip detects:
#   - an FO re-encryption check that leaves ciphertext bytes unverified (ePrint 2026/2239);
#   - a missing reduction before the ML-DSA inverse NTT, which wraps Int32 only on sign-aligned
#     inputs (ePrint 2026/1032, §5.2);
#   - signing-path bugs that only derandomised KATs catch (Bernstein 2026,
#     https://cr.yp.to/papers/mldsa-20260601.pdf).
# Included from runtests.jl, which defines h, groups and KAT_DIR.
using Random
import SHA

# ==================== ML-KEM: FO comparison covers every ciphertext byte ====================

# Byte `off` (0-based) of a d-bit-packed section: the bit of least weight within its coefficient.
# Flipping it perturbs one coefficient by the smallest step available in that byte, so decryption
# still returns m (asserted below) and only the re-encryption comparison can reject.
lowbit(off, d) = argmin(b -> (8off + b) % d, 0:7)

function tamper_positions(ct, nu, du, dv)
    [(p, 0x01 << (p <= nu ? lowbit(p - 1, du) : lowbit(p - 1 - nu, dv))) for p in eachindex(ct)]
end

@testset "ML-KEM implicit rejection, every ciphertext byte (ePrint 2026/2239)" begin
    for (name, Cat) in [("512", MLKEM.Category1), ("768", MLKEM.Category3), ("1024", MLKEM.Category5)]
        @testset "ML-KEM-$name" begin
            pk, sk = Cat.kyber_kem_keypair_derand(collect(UInt8, 1:64))
            m = fill(0x5a, 32)
            ct, ss = Cat.kyber_kem_enc_derand(pk, m)
            @test Cat.kyber_kem_dec(ct, sk) == ss
            z = sk[end-31:end]
            skc = sk[1:Cat.KYBER_INDCPA_SECRETKEYBYTES]
            wrong_key = Int[]; m_changed = Int[]
            for (p, bit) in tamper_positions(ct, Cat.KYBER_POLYVECCOMPRESSEDBYTES, Cat.KYBER_DU, Cat.KYBER_DV)
                c = copy(ct); c[p] ⊻= bit
                Cat.kyber_kem_dec(c, sk) == SHA.shake256(vcat(z, c), UInt64(32)) || push!(wrong_key, p)
                Cat.kyber_indcpa_dec(c, skc) == m || push!(m_changed, p)
            end
            @test isempty(wrong_key)                 # K̄ = SHAKE256(z ‖ c′, 32) at all positions
            @test isempty(m_changed)                 # each flip reaches the comparison with m′ = m
            # The single-coordinate model of 2026/2239: sweep the last d_v-bit coefficient of c₂.
            dv = Cat.KYBER_DV
            wrong_key = Int[]
            for v in 0:(1 << dv) - 1
                c = copy(ct); c[end] = (c[end] & (0xff >> dv)) | UInt8(v << (8 - dv))
                c == ct || Cat.kyber_kem_dec(c, sk) == SHA.shake256(vcat(z, c), UInt64(32)) || push!(wrong_key, v)
            end
            @test isempty(wrong_key)
        end
    end

    @testset "X-Wing, ML-KEM-768 part of the ciphertext" begin
        pk, sk = XWing.xwing_keypair_derand(collect(UInt8, 1:32))
        ct, ss = XWing.xwing_encaps_derand(pk, collect(UInt8, 65:128))
        @test XWing.xwing_decaps(ct, sk) == ss
        k = XWing.expand(sk)
        ctX = ct[1089:1120]
        ssX = PQJulia.X25519.x25519(k.skX, ctX)
        C = MLKEM.Category3
        wrong_key = Int[]
        for (p, bit) in tamper_positions(ct[1:1088], C.KYBER_POLYVECCOMPRESSEDBYTES, C.KYBER_DU, C.KYBER_DV)
            c = copy(ct); c[p] ⊻= bit
            ssM = SHA.shake256(vcat(k.skM[end-31:end], c[1:1088]), UInt64(32))
            XWing.xwing_decaps(c, sk) == XWing.combiner(ssM, ssX, ctX, k.pkX) || push!(wrong_key, p)
        end
        @test isempty(wrong_key)
    end
end

# ==================== ML-DSA: reductions before the inverse NTT ====================

const DQ = Int64(MLDSA.Q)
const RMOD = mod(Int64(2)^32, DQ)                         # Montgomery factor R mod q
brv8(i) = Int(bitreverse(UInt8(i)))
# FIPS 204 §7.5 with ζ = 1753, independent of the ZETAS table: NTT(w)ᵢ = Σⱼ wⱼ ζ^((2brv(i)+1)j).
const INTT_REF = [mod(invmod(256, DQ) * powermod(1753, mod(-(2brv8(i) + 1) * j, 512), DQ), DQ)
                  for j in 0:255, i in 0:255]
ref_intt(a) = [mod(sum(Int128(INTT_REF[j, i]) * a[i] for i in 1:256), DQ) for j in 1:256]

# Exact product in Z[X]/(X²⁵⁶ + 1).
function negacyclic(a, b)
    c = zeros(Int128, 256)
    @inbounds for i in 1:256, j in 1:256
        k = i + j - 1
        k <= 256 ? (c[k] += Int128(a[i]) * b[j]) : (c[k-256] -= Int128(a[i]) * b[j])
    end
    c
end

# |invntt! output| bound. Before the final scaling |x| ≤ 256(q−1) when every input satisfies
# |x̂ᵢ| < q, and montgomery_reduce(a) = (a − tq)/2³² with |t| ≤ 2³¹, so
# |out| ≤ 41978·256(q−1)/2³² + q/2 < 4,211,178. Lemma 2(ii) of 2026/1032 gives 4,201,525 under the
# same contract; the extreme input tested below exceeds it.
const INVNTT_BOUND = fld(41978 * 256 * (DQ - 1), Int64(2)^32) + cld(DQ, 2)

# Â whose Montgomery products with v̂ all have sign s and magnitude in [3.4·10⁶, 3.9·10⁶], the
# sign-aligned worst case of 2026/1032 §5.2. montgomery_reduce(a) lies within q/2 + |a|/2³² of 0
# and |Â·v̂| < 9q², so each such residue is returned unchanged. Every entry lies in [0, q), so Â
# is a possible ExpandA output.
function aligned_A(rng, vhat, K, s)
    inv_or_0(x) = iszero(mod(x, DQ)) ? 0 : invmod(mod(x, DQ), DQ)
    [[Int32(mod(mod(s * rand(rng, 3_400_000:3_900_000) * RMOD, DQ) * inv_or_0(x), DQ)) for x in vhat[j]]
     for _ in 1:K, j in eachindex(vhat)]
end

function uniform_A(rng, K, L)
    A = [rand(rng, Int32(0):Int32(DQ - 1), 256) for _ in 1:K, _ in 1:L]
    for a in A; a[rand(rng, 1:256, 32)] .= Int32(DQ - 1) .- rand(rng, Int32(0):Int32(3), 32); end
    A
end

# keygen/sign/verify order: Σⱼ Â[i,j]∘v̂ⱼ (− ĉ∘t̂₁), poly_reduce!, invntt!, (+ s₂), poly_caddq!.
# These rows replay that sequence on chosen inputs; the library's own sites are covered by
# compute_w! below (signing) and by the probe testset (keygen, signing, verification).
function row_pipeline(Ahat, vhat, i; sub=nothing, add=nothing)
    acc = zeros(Int32, 256)
    for j in eachindex(vhat); MLDSA.poly_pointwise_acc!(acc, Ahat[i, j], vhat[j]); end
    sub === nothing || MLDSA.poly_sub!(acc, acc, sub)
    unreduced = copy(acc)
    MLDSA.poly_reduce!(acc)
    nttin = copy(acc)
    MLDSA.invntt!(acc)
    nttout = copy(acc)
    add === nothing || MLDSA.poly_add!(acc, acc, add)
    MLDSA.poly_caddq!(acc)
    (; unreduced, nttin, nttout, out=acc)
end

@testset "ML-DSA reductions before the inverse NTT (ePrint 2026/1032)" begin
    @testset "invntt! extreme input" begin
        # Index 0 of the Gentleman–Sande addition path sums all 256 inputs.
        x = fill(Int32(-(DQ - 1)), 256); x[256] = Int32(-2_141_441_839 + 255 * (DQ - 1))
        out = MLDSA.invntt!(copy(x))
        @test maximum(abs, x) < DQ
        @test mod.(out, DQ) == mod.(RMOD .* ref_intt(x), DQ)
        @test 4_201_525 < maximum(abs, out) <= INVNTT_BOUND
    end

    for (ps, Cat) in [("ML-DSA-44", MLDSA.Category2), ("ML-DSA-65", MLDSA.Category3), ("ML-DSA-87", MLDSA.Category5)]
        @testset "$ps" begin
            K, L, η, γ1, β, τ = Cat.K, Cat.L, Cat.ETA, Cat.GAMMA1, Cat.BETA, Cat.TAU
            rng = MersenneTwister(1032 + K)
            pm(v, n) = rand(rng, Int32[-v, v], n)
            s1 = [pm(η, 256) for _ in 1:L]; s2 = [pm(η, 256) for _ in 1:K]
            y = [pm(γ1 - 1, 256) for _ in 1:L]
            z = [pm(γ1 - β - 1, 256) for _ in 1:L]
            t1 = [rand(rng, Int32[0, 1023], 256) for _ in 1:K]
            c = zeros(Int32, 256); c[randperm(rng, 256)[1:τ]] .= pm(1, τ)
            chat = MLDSA.ntt!(copy(c))
            ct1 = [MLDSA.poly_pointwise!(zeros(Int32, 256), chat, MLDSA.ntt!(MLDSA.poly_shiftl!(copy(t)))) for t in t1]
            hat(v) = [MLDSA.ntt!(copy(p)) for p in v]

            # (site, v, subtracted before reduction, added after invntt!, exact extra term)
            sites = (("keygen t = As₁ + s₂", s1, i -> nothing, i -> s2[i], i -> s2[i]),
                     ("sign w = Ay", y, i -> nothing, i -> nothing, i -> zeros(Int32, 256)),
                     ("verify w′ = Az − ct₁2ᵈ", z, i -> ct1[i], i -> nothing,
                      i -> -negacyclic(c, Int64(2)^13 .* t1[i])))
            for (site, v, sub, add, extra) in sites
                vhat = hat(v)
                for (variant, aligned, Ahat) in (("uniform Â", false, uniform_A(rng, K, L)),
                                                 ("aligned Â, +", true, aligned_A(rng, vhat, K, 1)),
                                                 ("aligned Â, −", true, aligned_A(rng, vhat, K, -1)))
                    @testset "$site, $variant" begin
                        for i in 1:K
                            expect = sum(negacyclic(ref_intt(Ahat[i, j]), v[j]) for j in 1:L) .+ extra(i)
                            r = row_pipeline(Ahat, vhat, i; sub=sub(i), add=add(i))
                            # Without poly_reduce! these rows wrap Int32 at index 0 of invntt!.
                            aligned && @test abs(sum(Int64, r.unreduced)) > typemax(Int32)
                            @test maximum(abs, r.nttin) < DQ
                            @test maximum(abs, r.nttout) <= INVNTT_BOUND
                            @test r.out == mod.(expect, DQ)
                        end
                    end
                end
            end

            @testset "signer compute_w!, aligned Â" begin
                Ahat = aligned_A(rng, hat(y), K, 1)
                w1 = [zeros(Int32, 256) for _ in 1:K]; w0 = [zeros(Int32, 256) for _ in 1:K]
                Cat.compute_w!(w1, w0, Ahat, y, zeros(Int32, 256))
                for i in 1:K
                    expect = mod.(sum(negacyclic(ref_intt(Ahat[i, j]), y[j]) for j in 1:L), DQ)
                    @test mod.(Int64.(w1[i]) .* 2Cat.GAMMA2 .+ w0[i], DQ) == expect
                end
            end

            # Theorem 1 of 2026/1032: |c·u| ≤ τβ and τβ + INVNTT_BOUND < q, so ĉ∘û through invntt!
            # returns the integer product c·u, not only its residue.
            @testset "sparse products c·s, c·t₀ are exact" begin
                cw = zeros(Int32, 256); cw[1:τ] .= 1             # attains τβ against a constant u
                for (u, bound) in ((s1[1], τ * η), (fill(Int32(η), 256), τ * η),
                                   (rand(rng, Int32[-4095, 4096], 256), τ * 4096), (fill(Int32(4096), 256), τ * 4096))
                    uhat = MLDSA.ntt!(copy(u))
                    for cc in (c, cw)
                        out = MLDSA.invntt!(MLDSA.poly_pointwise!(zeros(Int32, 256), MLDSA.ntt!(copy(cc)), uhat))
                        @test out == negacyclic(cc, u)
                        @test maximum(abs, out) <= bound
                    end
                end
            end
        end
    end
end

# The library's keygen and verification sites inline the sequence on ExpandA output, so chosen
# inputs cannot reach them. These modules compile the shipped dilithium_level.jl with invntt!
# wrapped to record its largest input; with every reduction in place that stays below q, while a
# missing one accumulates up to L·q. Outputs must equal the shipped modules' byte for byte.
module DSAProbe
import PQJulia, SHA
using Random
const MAXIN = Ref(0)
for (category, params) in PQJulia.MLDSA.CATEGORY_PARAMS
    @eval module $category
    import PQJulia.MLDSA: Q, N, D, ZETAS, montgomery_reduce, reduce32, caddq, freeze, ntt!,
                          poly_pointwise!, poly_pointwise_acc!, poly_add!, poly_sub!, poly_reduce!,
                          poly_caddq!, poly_shiftl!, poly_chknorm, poly_uniform!, power2round,
                          derived_sizes, wipe!
    import PQJulia, PQJulia.Keccak, SHA
    using Random
    import ..MAXIN
    invntt!(a::Vector{Int32}) = (MAXIN[] = max(MAXIN[], Int(maximum(abs, a))); PQJulia.MLDSA.invntt!(a))
    const (CAT_NUM, K, L, TAU, ETA, OMEGA, _GAMMA2_DIV, _LG_GAMMA1, CTILDEBYTES) = $params
    const BETA = Int32(TAU * ETA)
    const GAMMA1 = Int32(1 << _LG_GAMMA1)
    const GAMMA2 = Int32(div(Q - 1, _GAMMA2_DIV))
    const SEEDBYTES = 32
    const CRHBYTES = 64
    const TRBYTES = 64
    const SIZES = derived_sizes(K, L, ETA, OMEGA, _GAMMA2_DIV, _LG_GAMMA1, CTILDEBYTES)
    const POLYETA_PACKED = SIZES.polyeta_packed
    const POLYZ_PACKED = SIZES.polyz_packed
    const POLYW1_PACKED = SIZES.polyw1_packed
    const POLYT1_PACKED = SIZES.polyt1_packed
    const POLYT0_PACKED = SIZES.polyt0_packed
    const PK_BYTES = SIZES.pk_bytes
    const SK_BYTES = SIZES.sk_bytes
    const SIG_BYTES = SIZES.sig_bytes
    const IDENTIFIER = "ML-DSA-$K$L"
    include(joinpath(pkgdir(PQJulia), "src", "dilithium_level.jl"))
    end
end
end

@testset "ML-DSA library sites: every invntt! input below q (ePrint 2026/1032)" begin
    for (name, Cat) in [(:Category2, MLDSA.Category2), (:Category3, MLDSA.Category3), (:Category5, MLDSA.Category5)]
        P = getfield(DSAProbe, name)
        @testset "$(Cat.IDENTIFIER)" begin
            DSAProbe.MAXIN[] = 0
            same = true
            for seed in 1:4
                xi = fill(UInt8(seed), 32); msg = Vector{UInt8}("probe $seed"); rnd = fill(UInt8(7seed), 32)
                pk, sk = P.dilithium_keygen_derand(xi)
                sig = P.dilithium_sign_derand(msg, sk, rnd)
                same &= (pk, sk) == Cat.dilithium_keygen_derand(xi) && sig == Cat.dilithium_sign_derand(msg, sk, rnd)
                same &= P.dilithium_verify(msg, sig, pk) && !P.dilithium_verify(msg, sig .⊻ [0x01; zeros(UInt8, length(sig) - 1)], pk)
            end
            @test same                                           # the probe is the shipped code
            @test 0 < DSAProbe.MAXIN[] < MLDSA.Q
        end
    end
end

# ==================== ML-DSA: public signing wrappers against ACVP ====================

# The ACVP block in runtests.jl calls the *_derand and internal entry points. Here the public
# dilithium_sign / dilithium_sign_prehash with hedged=false must reproduce the deterministic vectors.
let det = groups("mldsa_siggen_prompt.json",
                 g -> g["signatureInterface"] == "external" && g["deterministic"] && !g["externalMu"])
    for (ps, Cat) in [("ML-DSA-44", MLDSA.Category2), ("ML-DSA-65", MLDSA.Category3), ("ML-DSA-87", MLDSA.Category5)]
        @testset "$ps public sign wrappers, ACVP deterministic" begin
            for g in filter(g -> g["parameterSet"] == ps, det)
                @testset "$(g["preHash"]) ($(length(g["tests"])))" begin
                    for tc in g["tests"]
                        msg, sk, ctx = h(tc["message"]), h(tc["sk"]), h(tc["context"])
                        sig = g["preHash"] == "preHash" ?
                              Cat.dilithium_sign_prehash(msg, sk, tc["hashAlg"]; hedged=false, context=ctx) :
                              Cat.dilithium_sign(msg, sk; hedged=false, context=ctx)
                        @test sig == h(tc["signature"])
                    end
                end
            end
        end
    end
end
