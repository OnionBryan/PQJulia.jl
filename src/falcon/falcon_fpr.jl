# src/falcon/falcon_fpr.jl
# ============================================================================
# Falcon signing without hardware floating point: Pornin's integer emulation of
# IEEE-754 binary64 (fpr.c of the reference, FALCON_FPEMU), branch-free, and the
# FFT, ffLDL tree, ffSampling and SamplerZ written on it. Every operation follows
# the Float64 signer (falcon.jl, falcon_fft.jl, falcon_sampler.jl) in the same
# order, so the signatures are bit-identical; see README, "Falcon signing".
# ============================================================================
module FalconFpr

import ..FalconFFT as FF
import ..FalconSampler as FS
import ...Wipe: wipe!

# ── binary64 as raw bits ────────────────────────────────────────────────────
struct Fpr
    b::UInt64
end
fpr(x::Float64) = Fpr(reinterpret(UInt64, x))
Base.Float64(x::Fpr) = reinterpret(Float64, x.b)

const M52 = (UInt64(1) << 52) - 1
const ZERO = fpr(0.0); const ONE = fpr(1.0); const TWO = fpr(2.0)
const PTWO63 = fpr(2.0^63)

# Shifts by 0..63 without a data-dependent shift count above 31.
@inline ursh(x::UInt64, n::Int) = (x ⊻= (x ⊻ (x >> 32)) & (UInt64(0) - ((n >> 5) % UInt64)); x >> (n & 31))
@inline ulsh(x::UInt64, n::Int) = (x ⊻= (x ⊻ (x << 32)) & (UInt64(0) - ((n >> 5) % UInt64)); x << (n & 31))
@inline irsh(x::Int64, n::Int) = (x ⊻= (x ⊻ (x >> 32)) & -((n >> 5) % Int64); x >> (n & 31))

# Value m·2^e, m in [2^54, 2^55) with two rounding bits; rounds to nearest even.
@inline function mk(s::Int, e::Int, m::UInt64)
    e += 1076
    t = (e % UInt32) >> 31
    m &= (t % UInt64) - UInt64(1)
    t = (m >> 54) % UInt32
    e &= -(t % Int)
    x = (((s % UInt64) << 63) | (m >> 2)) + (((e % UInt32) % UInt64) << 52)
    f = (m % UInt32) & 0x7
    Fpr(x + ((0xC8 >> f) & 1))
end

# Left-normalize m to [2^63, 2^64), adjusting e.
@inline function norm64(m::UInt64, e::Int)
    e -= 63
    for (sh, k) in ((32, 5), (16, 4), (8, 3), (4, 2), (2, 1))
        nt = (m >> (64 - sh)) % UInt32
        nt = (nt | (UInt32(0) - nt)) >> 31
        m ⊻= (m ⊻ (m << sh)) & ((nt % UInt64) - UInt64(1))
        e += (nt << k) % Int
    end
    nt = (m >> 63) % UInt32
    m ⊻= (m ⊻ (m << 1)) & ((nt % UInt64) - UInt64(1))
    m, e + (nt % Int)
end

neg(x::Fpr) = Fpr(x.b ⊻ (UInt64(1) << 63))
# x/2, exact for normal x; ±0 stays ±0.
half(x::Fpr) = Fpr(x.b - ((((((x.b >> 52) & 0x7FF) + 0x7FF) >> 11)) << 52))

function add(x::Fpr, y::Fpr)
    xb, yb = x.b, y.b
    m = (UInt64(1) << 63) - UInt64(1)
    za = (xb & m) - (yb & m)
    cs = ((za >> 63) % UInt32) | ((UInt32(1) - (((UInt64(0) - za) >> 63) % UInt32)) & ((xb >> 63) % UInt32))
    m = (xb ⊻ yb) & (UInt64(0) - (cs % UInt64))
    xb ⊻= m; yb ⊻= m
    ex = (xb >> 52) % Int; sx = ex >> 11; ex &= 0x7FF
    xu = ((xb & M52) | ((((ex + 0x7FF) >> 11) % UInt64) << 52)) << 3
    ex -= 1078
    ey = (yb >> 52) % Int; sy = ey >> 11; ey &= 0x7FF
    yu = ((yb & M52) | ((((ey + 0x7FF) >> 11) % UInt64) << 52)) << 3
    ey -= 1078
    cc = ex - ey
    yu &= UInt64(0) - ((((cc - 60) % UInt32) >> 31) % UInt64)
    cc &= 63
    m = ulsh(UInt64(1), cc) - UInt64(1)
    yu |= (yu & m) + m
    yu = ursh(yu, cc)
    xu += yu - ((yu << 1) & (UInt64(0) - ((sx ⊻ sy) % UInt64)))
    xu, ex = norm64(xu, ex)
    xu |= ((((xu % UInt32) & 0x1FF) + 0x1FF) % UInt64)
    mk(sx, ex + 9, xu >> 9)
end
sub(x::Fpr, y::Fpr) = add(x, neg(y))

function mul(x::Fpr, y::Fpr)
    xu = (x.b & M52) | (UInt64(1) << 52)
    yu = (y.b & M52) | (UInt64(1) << 52)
    p = UInt128(xu) * UInt128(yu)
    lo = (p % UInt64) & ((UInt64(1) << 50) - UInt64(1))
    zu = (p >> 50) % UInt64
    zu |= (lo + ((UInt64(1) << 50) - UInt64(1))) >> 50              # sticky bit
    zv = (zu >> 1) | (zu & UInt64(1))
    w = zu >> 55
    zu ⊻= (zu ⊻ zv) & (UInt64(0) - w)
    ex = ((x.b >> 52) & 0x7FF) % Int; ey = ((y.b >> 52) & 0x7FF) % Int
    e = ex + ey - 2100 + (w % Int)
    s = ((x.b ⊻ y.b) >> 63) % Int
    d = ((ex + 0x7FF) & (ey + 0x7FF)) >> 11
    zu &= UInt64(0) - (d % UInt64)
    mk(s, e, zu)
end
sqr(x::Fpr) = mul(x, x)

function div(x::Fpr, y::Fpr)
    xu = (x.b & M52) | (UInt64(1) << 52)
    yu = (y.b & M52) | (UInt64(1) << 52)
    q = UInt64(0)
    for _ in 1:55
        b = ((xu - yu) >> 63) - UInt64(1)
        xu -= b & yu
        q |= b & UInt64(1)
        xu <<= 1; q <<= 1
    end
    q |= (xu | (UInt64(0) - xu)) >> 63
    q2 = (q >> 1) | (q & UInt64(1))
    w = q >> 55
    q ⊻= (q ⊻ q2) & (UInt64(0) - w)
    ex = ((x.b >> 52) & 0x7FF) % Int; ey = ((y.b >> 52) & 0x7FF) % Int
    e = ex - ey - 55 + (w % Int)
    s = ((x.b ⊻ y.b) >> 63) % Int
    d = (ex + 0x7FF) >> 11
    s &= d; e &= -d
    q &= UInt64(0) - (d % UInt64)
    mk(s, e, q)
end

function Base.sqrt(x::Fpr)
    xu = (x.b & M52) | (UInt64(1) << 52)
    ex = ((x.b >> 52) & 0x7FF) % Int
    e = ex - 1023
    xu += xu & (UInt64(0) - ((e & 1) % UInt64))
    e >>= 1
    xu <<= 1
    q = UInt64(0); s = UInt64(0); r = UInt64(1) << 53
    for _ in 1:54
        t = s + r
        b = ((xu - t) >> 63) - UInt64(1)
        s += (r << 1) & b
        xu -= t & b
        q += r & b
        xu <<= 1; r >>= 1
    end
    q <<= 1
    q |= (xu | (UInt64(0) - xu)) >> 63
    q &= UInt64(0) - ((((ex + 0x7FF) >> 11)) % UInt64)
    mk(0, e - 54, q)
end

# Int64 → binary64, rounded.
function fpr_of(i::Int64)
    s = ((i % UInt64) >> 63) % Int
    i ⊻= -s; i += s
    m = i % UInt64
    m, e = norm64(m, 9)
    m |= ((((m % UInt32) & 0x1FF) + 0x1FF) % UInt64)
    m >>= 9
    t = (((i | -i) % UInt64) >> 63)
    m &= UInt64(0) - t
    e &= -(t % Int)
    mk(s, e, m)
end
fpr_of(i::Integer) = fpr_of(Int64(i))

# Round half to even, as Julia's round(Int, x).
function rint(x::Fpr)
    m = ((x.b << 10) | (UInt64(1) << 62)) & ((UInt64(1) << 63) - UInt64(1))
    e = 1085 - (((x.b >> 52) % Int) & 0x7FF)
    m &= UInt64(0) - ((((e - 64) % UInt32) >> 31) % UInt64)
    e &= 63
    d = ulsh(m, 63 - e)
    dd = (d % UInt32) | (((d >> 32) % UInt32) & 0x1FFFFFFF)
    f = ((d >> 61) % UInt32) | ((dd | (UInt32(0) - dd)) >> 31)
    m = ursh(m, e) + (((0xC8 % UInt32) >> f) & UInt32(1)) % UInt64
    s = (x.b >> 63) % Int64
    ((m % Int64) ⊻ -s) + s
end

# floor, for |x| < 2^63; floor(−0) = 0 as in Julia.
function Base.floor(x::Fpr)
    e = ((x.b >> 52) % Int) & 0x7FF
    t = ((x.b >> 63) % Int64) & ((e + 0x7FF) >> 11)
    xi = (((x.b << 10) | (UInt64(1) << 62)) & ((UInt64(1) << 63) - UInt64(1))) % Int64
    xi = (xi ⊻ -t) + t
    cc = 1085 - e
    xi = irsh(xi, cc & 63)
    xi ⊻= (xi ⊻ -t) & -((((63 - cc) % UInt32) >> 31) % Int64)
    xi
end

function lt(x::Fpr, y::Fpr)
    sx = x.b % Int64; sy = y.b % Int64
    sy &= ~((sx ⊻ sy) >> 63)
    cc0 = ((sx - sy) >> 63) & 1
    cc1 = ((sy - sx) >> 63) & 1
    (cc0 ⊻ ((cc0 ⊻ cc1) & (((x.b & y.b) >> 63) % Int64))) == 1
end

# ── complex numbers on Fpr, with Julia's Complex{Float64} operation order ──
struct CF
    re::Fpr
    im::Fpr
end
cf(z::ComplexF64) = CF(fpr(real(z)), fpr(imag(z)))
Base.ComplexF64(z::CF) = ComplexF64(Float64(z.re), Float64(z.im))
Base.iszero(z::CF) = (z.re.b | z.im.b) == 0
cadd(a::CF, b::CF) = CF(add(a.re, b.re), add(a.im, b.im))
csub(a::CF, b::CF) = CF(sub(a.re, b.re), sub(a.im, b.im))
cmul(a::CF, b::CF) = CF(sub(mul(a.re, b.re), mul(a.im, b.im)), add(mul(a.re, b.im), mul(a.im, b.re)))
cconj(a::CF) = CF(a.re, neg(a.im))
cneg(a::CF) = CF(neg(a.re), neg(a.im))
cdivr(a::CF, x::Fpr) = CF(div(a.re, x), div(a.im, x))
chalf(a::CF) = CF(half(a.re), half(a.im))
# a / b for b with zero imaginary part: Julia's robust division takes the a·(1/c) path.
cdiv_real(a::CF, b::CF) = (t = div(ONE, b.re); CF(mul(a.re, t), mul(a.im, t)))

# ── FFT over ℝ[x]/(xⁿ+1), as FalconFFT ──────────────────────────────────────
const ROOTS = Dict{Int,Vector{CF}}()
const ROOTS_LOCK = ReentrantLock()
roots(n) = lock(() -> get!(() -> cf.(FF.roots(n)), ROOTS, n), ROOTS_LOCK)

function split_fft(F::Vector{CF})
    n = length(F); h = n ÷ 2; ζ = roots(n)
    F0 = Vector{CF}(undef, h); F1 = Vector{CF}(undef, h)
    @inbounds for k in 1:h
        F0[k] = chalf(cadd(F[k], F[k+h]))
        F1[k] = chalf(cmul(csub(F[k], F[k+h]), cconj(ζ[k])))
    end
    F0, F1
end
function merge_fft(F0::Vector{CF}, F1::Vector{CF})
    h = length(F0); n = 2h; ζ = roots(n)
    F = Vector{CF}(undef, n)
    @inbounds for k in 1:h
        F[k]   = cadd(F0[k], cmul(ζ[k], F1[k]))
        F[k+h] = csub(F0[k], cmul(ζ[k], F1[k]))
    end
    F
end
function fft(f::Vector{Fpr})
    length(f) == 1 && return [CF(f[1], ZERO)]
    f0, f1 = f[1:2:end], f[2:2:end]
    a, b = fft(f0), fft(f1)
    out = merge_fft(a, b)
    wipe!(f0, f1, a, b)
    out
end
function ifft(F::Vector{CF})
    length(F) == 1 && return [F[1]]
    F0, F1 = split_fft(F)
    a, b = ifft(F0), ifft(F1)
    out = Vector{CF}(undef, length(F))
    @inbounds for i in eachindex(a); out[2i-1] = a[i]; out[2i] = b[i]; end
    wipe!(F0, F1, a, b)
    out
end

# ── ffLDL tree (falcon.jl ffldl) ────────────────────────────────────────────
mutable struct Leaf
    σ::Fpr
end
struct Node
    l10::Vector{CF}
    t0::Union{Node,Leaf}
    t1::Union{Node,Leaf}
end

function ffldl(g00::Vector{CF}, g01::Vector{CF}, g11::Vector{CF}, σ::Fpr)
    n = length(g00)
    l10 = [cdiv_real(cconj(g01[i]), g00[i]) for i in 1:n]
    d11 = [csub(g11[i], cmul(cmul(l10[i], cconj(l10[i])), g00[i])) for i in 1:n]
    if n > 2
        a0, a1 = split_fft(g00); b0, b1 = split_fft(d11)
        node = Node(l10, ffldl(a0, a1, a0, σ), ffldl(b0, b1, b0, σ))
        wipe!(d11, a0, a1, b0, b1)
        return node
    end
    node = Node(l10, Leaf(div(σ, sqrt(g00[1].re))), Leaf(div(σ, sqrt(d11[1].re))))
    wipe!(d11)
    node
end
leaves(T, out=Float64[]) = T isa Leaf ? push!(out, Float64(T.σ)) : (leaves(T.t0, out); leaves(T.t1, out))
wipe!(T::Leaf) = (T.σ = ZERO; T)
wipe!(T::Node) = (wipe!(T.l10); wipe!(T.t0); wipe!(T.t1); T)

# ── SamplerZ (falcon_sampler.jl) ────────────────────────────────────────────
const LN2 = fpr(FS.LN2); const ILN2 = fpr(FS.ILN2); const INV_2SIGMA2 = fpr(FS.INV_2SIGMA2)

function approxexp(x::Fpr, ccs::Fpr)
    y = Int128(FS.APPROXEXP_C[1])
    z = Int128(floor(mul(x, PTWO63)))
    @inbounds for i in 2:13; y = Int128(FS.APPROXEXP_C[i]) - ((z * y) >> 63); end
    c = Int128(floor(mul(ccs, PTWO63)))
    c = ifelse(ccs.b == ONE.b, Int128(1) << 63, c)                 # floor(2^63) overflows Int64
    ((c << 1) * y) >> 63
end

function berexp(x::Fpr, ccs::Fpr, r)
    s = floor(mul(x, ILN2)); rr = sub(x, mul(fpr_of(s), LN2)); s = ifelse(s > 63, 63, s)
    z = (approxexp(rr, ccs) - 1) >> s; w = 0
    for i in 56:-8:0
        p = Int(byte(r)); w = p - Int((z >> i) & 0xFF)
        w != 0 && break
    end
    w < 0
end

# FalconSampler.basesampler and single bytes, with the randomness wiped after use.
byte(r) = (b = FS.randbytes(r, 1); v = b[1]; wipe!(b); v)
function basesampler(r)
    bytes = FS.randbytes(r, 9)
    u = UInt128(0); for i in 1:9; u |= UInt128(bytes[i]) << (8 * (i - 1)); end
    wipe!(bytes)
    z0 = 0
    @inbounds for c in FS.RCDT; z0 += Int(u < c); end
    z0
end

function samplerz(μ::Fpr, σ::Fpr, σmin::Fpr, r)
    s = floor(μ); rr = sub(μ, fpr_of(s))
    dss = div(ONE, mul(TWO, sqr(σ))); ccs = div(σmin, σ)
    while true
        z0 = basesampler(r)
        b = Int(byte(r)) & 1
        z = b + (2b - 1) * z0
        x = sub(mul(sqr(sub(fpr_of(z), rr)), dss), mul(fpr_of(z0^2), INV_2SIGMA2))
        berexp(x, ccs, r) && return z + s
    end
end

# ── Per-key setup, target and one ffSampling attempt (falcon.jl) ────────────
struct Setup
    b00::Vector{CF}; b01::Vector{CF}; b10::Vector{CF}; b11::Vector{CF}
    T::Node
end

gram(a, b, c, d) = [cadd(cmul(a[i], cconj(c[i])), cmul(b[i], cconj(d[i]))) for i in eachindex(a)]

"Signing data for `sk`; refuses keys with an ffLDL leaf outside [σmin, σmax]."
function setup(sk, σ::Float64, σmin::Float64)
    pg, pf, pG, pF = fpr_of.(sk.g), neg.(fpr_of.(sk.f)), fpr_of.(sk.G), neg.(fpr_of.(sk.F))
    b00, b01, b10, b11 = fft(pg), fft(pf), fft(pG), fft(pF)
    wipe!(pg, pf, pG, pF)
    g00 = gram(b00, b01, b00, b01); g01 = gram(b00, b01, b10, b11); g11 = gram(b10, b11, b10, b11)
    T = ffldl(g00, g01, g11, fpr(σ))
    wipe!(g00, g01, g11)
    lo, hi = fpr(σmin), fpr(FS.MAX_SIGMA)
    ok = all(v -> !lt(v, lo) && !lt(hi, v), leafvals(T))
    ok || (wipe!(T); wipe!(b00, b01, b10, b11);
           throw(ArgumentError("Falcon key has an ffLDL leaf outside [σmin, σmax]; refusing to sign")))
    Setup(b00, b01, b10, b11, T)
end
leafvals(T, out=Fpr[]) = T isa Leaf ? push!(out, T.σ) : (leafvals(T.t0, out); leafvals(T.t1, out))
wipe!(gs::Setup) = (wipe!(gs.b00, gs.b01, gs.b10, gs.b11); wipe!(gs.T); gs)

const Q = fpr(12289.0)
function target(gs::Setup, c::AbstractVector{<:Integer})
    pc = fpr_of.(c); F = fft(pc)
    t0 = [cdivr(cmul(F[i], gs.b11[i]), Q) for i in eachindex(F)]
    t1 = [cdivr(cmul(cneg(F[i]), gs.b01[i]), Q) for i in eachindex(F)]
    wipe!(pc, F)
    t0, t1
end

function ffsample(t0::Vector{CF}, t1::Vector{CF}, T::Union{Node,Leaf}, σmin::Fpr, r)
    if T isa Leaf
        return [CF(fpr_of(samplerz(t0[1].re, T.σ, σmin, r)), ZERO)],
               [CF(fpr_of(samplerz(t1[1].re, T.σ, σmin, r)), ZERO)]
    end
    a, b = split_fft(t1)
    za, zb = ffsample(a, b, T.t1, σmin, r)
    z1 = merge_fft(za, zb)
    t0b = [cadd(t0[i], cmul(csub(t1[i], z1[i]), T.l10[i])) for i in eachindex(t0)]
    c, d = split_fft(t0b)
    zc, zd = ffsample(c, d, T.t0, σmin, r)
    z0 = merge_fft(zc, zd)
    wipe!(a, b, za, zb, t0b, c, d, zc, zd)
    z0, z1
end

function preimage(gs::Setup, c::AbstractVector{<:Integer}, t0, t1, σmin::Float64, r)
    z0, z1 = ffsample(t0, t1, gs.T, fpr(σmin), r)
    u0 = [cadd(cmul(z0[i], gs.b00[i]), cmul(z1[i], gs.b10[i])) for i in eachindex(z0)]
    u1 = [cadd(cmul(z0[i], gs.b01[i]), cmul(z1[i], gs.b11[i])) for i in eachindex(z0)]
    w0, w1 = ifft(u0), ifft(u1)
    v0 = [rint(w0[i].re) for i in eachindex(w0)]
    v1 = [rint(w1[i].re) for i in eachindex(w1)]
    s0, s1 = c .- v0, .-v1
    wipe!(z0, z1, u0, u1, w0, w1, v0, v1)
    s0, s1
end

end # module
