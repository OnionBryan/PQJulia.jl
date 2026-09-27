# Fixed-point Falcon keygen: ePrint 2023/290, pornin/ntrugen.
module FalconFxp

import ..FalconSampler as FS
import ..FalconNTRUGen as NG

const q = 12289

# ── fxr (ng_inner.h) ─────────────────────────────────────────────────────────
fxr_of(j::Integer) = Int64(j) << 32
fxr_mul(x::Int64, y::Int64) = ((Int128(x) * y) >> 32) % Int64
fxr_sqr(x::Int64) = fxr_mul(x, x)
fxr_round(x::Int64) = Int32((x + Int64(0x80000000)) >> 32)
fxr_div2e(x::Int64, n) = (x + ((Int64(1) << n) >> 1)) >> n

# inner_fxr_div (ng_fxp.c): |x|·2³²/|y| bit by bit, rounded, then signed.
function fxr_div(x::Int64, y::Int64)
    ux = reinterpret(UInt64, x); uy = reinterpret(UInt64, y)
    sx = ux >> 63; ux = (ux ⊻ -sx) + sx
    sy = uy >> 63; uy = (uy ⊻ -sy) + sy
    qq = UInt64(0); num = ux >> 31
    for i in 63:-1:33
        b = 1 - ((num - uy) >> 63)
        qq |= b << i; num -= uy & -b; num <<= 1; num |= (ux >> (i - 33)) & 1
    end
    for i in 32:-1:0
        b = 1 - ((num - uy) >> 63)
        qq |= b << i; num -= uy & -b; num <<= 1
    end
    qq += 1 - ((num - uy) >> 63)
    sx ⊻= sy
    reinterpret(Int64, (qq ⊻ -sx) + sx)
end

# ── fxc and the FFT (ng_fxp.c vect_*): re in [1, n/2], im in [n/2+1, n] ────
fxc_mul(ar, ai, br, bi) = (z0 = fxr_mul(ar, br); z1 = fxr_mul(ai, bi); z2 = fxr_mul(ar + ai, br + bi);
                           (z0 - z1, z2 - (z0 + z1)))

brev10(k) = foldl((a, i) -> (a << 1) | ((k >> i) & 1), 0:9; init = 0)
# GM_TAB[k] = e^{iπ·bitrev₁₀(k)/1024} in fxr, rounded to nearest (equal to ntrugen's table).
const GM_TAB = setprecision(BigFloat, 256) do
    [(round(Int64, cospi(big(brev10(k)) / 1024) * big(2)^32),
      round(Int64, sinpi(big(brev10(k)) / 1024) * big(2)^32)) for k in 0:1023]
end

function vect_FFT!(f::Vector{Int64})
    n = length(f); logn = trailing_zeros(n); hn = n >> 1; t = hn
    for lm in 1:logn-1
        m = 1 << lm; ht = t >> 1; j0 = 0
        for i in 0:(m >> 1)-1
            sr, si = GM_TAB[m + i + 1]
            for j in j0:j0+ht-1
                xr, xi = f[j+1], f[j+hn+1]
                yr, yi = fxc_mul(sr, si, f[j+ht+1], f[j+ht+hn+1])
                f[j+1] = xr + yr; f[j+hn+1] = xi + yi
                f[j+ht+1] = xr - yr; f[j+ht+hn+1] = xi - yi
            end
            j0 += t
        end
        t = ht
    end
    f
end

function vect_iFFT!(f::Vector{Int64})
    n = length(f); logn = trailing_zeros(n); hn = n >> 1; ht = 1
    for lm in logn-1:-1:1
        m = 1 << lm; t = ht << 1; j0 = 0
        for i in 0:(m >> 1)-1
            sr, si = GM_TAB[m + i + 1]; si = -si
            for j in j0:j0+ht-1
                xr, xi = f[j+1], f[j+hn+1]
                yr, yi = f[j+ht+1], f[j+ht+hn+1]
                f[j+1] = fxr_div2e(xr + yr, 1); f[j+hn+1] = fxr_div2e(xi + yi, 1)
                f[j+ht+1], f[j+ht+hn+1] = fxc_mul(sr, si, fxr_div2e(xr - yr, 1), fxr_div2e(xi - yi, 1))
            end
            j0 += t
        end
        ht = t
    end
    f
end

function vect_mul_fft!(a, b)
    hn = length(a) >> 1
    for u in 1:hn
        a[u], a[u+hn] = fxc_mul(a[u], a[u+hn], b[u], b[u+hn])
    end
    a
end
vect_adj_fft!(a) = (hn = length(a) >> 1; a[hn+1:end] .= .-a[hn+1:end]; a)
function vect_mul_autoadj_fft!(a, b)
    hn = length(a) >> 1
    for u in 1:hn; a[u] = fxr_mul(a[u], b[u]); a[u+hn] = fxr_mul(a[u+hn], b[u]); end
    a
end
# d = 2^e / (|a|² + |b|²), self-adjoint (real) in FFT form.
function vect_invnorm_fft(a, b, e)
    hn = length(a) >> 1
    [fxr_div(fxr_of(1 << e), fxr_sqr(a[u]) + fxr_sqr(a[u+hn]) + fxr_sqr(b[u]) + fxr_sqr(b[u+hn])) for u in 1:hn]
end

# v·2³²/2^sh rounded to nearest, for BigInt or small integer coefficients.
to_fxr(v, sh) = Int64(div(big(v) << 32 + (big(1) << sh >> 1), big(1) << sh, RoundDown))

# ── (f, g) sampling (ng_gauss.c, ng_falcon.c tables) ─────────────────────────
const GAUSS_256 = UInt16[24,
    1, 3, 6, 11, 22, 40, 73, 129, 222, 371, 602, 950, 1460, 2183, 3179, 4509,
    6231, 8395, 11032, 14150, 17726, 21703, 25995, 30487, 35048, 39540, 43832, 47809, 51385, 54503, 57140, 59304,
    61026, 62356, 63352, 64075, 64585, 64933, 65164, 65313, 65406, 65462, 65495, 65513, 65524, 65529, 65532, 65534]
const GAUSS_512 = UInt16[17,
    1, 4, 11, 28, 65, 146, 308, 615, 1164, 2083, 3535, 5692, 8706, 12669, 17574, 23285,
    29542, 35993, 42250, 47961, 52866, 56829, 59843, 62000, 63452, 64371, 64920, 65227, 65389, 65470, 65507, 65524,
    65531, 65534]
const GAUSS_1024 = UInt16[12,
    2, 8, 28, 94, 280, 742, 1761, 3753, 7197, 12472, 19623, 28206, 37329, 45912, 53063, 58338,
    61782, 63774, 64793, 65255, 65441, 65507, 65527, 65533]

# prng_buffer: 16-bit little-endian draws from 512-byte chunks of the RNG.
mutable struct U16Source
    rng::Any
    buf::Vector{UInt8}
    ptr::Int
end
U16Source(rng) = U16Source(rng, UInt8[], 512)
function next_u16(s::U16Source)
    if s.ptr > 510
        s.buf = FS.randbytes(s.rng, 512); s.ptr = 0
    end
    x = UInt32(s.buf[s.ptr+1]) | UInt32(s.buf[s.ptr+2]) << 8; s.ptr += 2
    x
end

# One table draw: −kmax + #{k : tab[k] < x}.
tabdraw(tab, x) = -Int(tab[1]) + count(k -> tab[k] < x, 2:2*tab[1]+1)

"gauss_sample_poly / gauss_sample_poly_reduced: n coefficients, resampled until their sum is odd."
function gauss_sample_poly(n::Int, rng)
    logn = trailing_zeros(n)
    2 <= logn <= 10 || throw(ArgumentError("fixed-point keygen supports n = 4..1024"))
    src = U16Source(rng)
    tab = logn == 10 ? GAUSS_1024 : logn == 9 ? GAUSS_512 : GAUSS_256
    while true
        f = Int[]
        if logn >= 8
            for _ in 1:n; push!(f, tabdraw(tab, next_u16(src))); end
        else
            while length(f) < n
                y = sum(tabdraw(tab, next_u16(src)) for _ in 1:1 << (8 - logn))
                -127 <= y <= 127 && push!(f, y)
            end
        end
        isodd(sum(f)) && return f
    end
end

# ── Checks of Falcon_keygen (ng_falcon.c) ────────────────────────────────────
# ‖(f,g)‖² < 16823, then q²‖(f*, g*)/(ff* + gg*)‖² < 1.17²·q (fxr 72251709809335 = ⌊1.17²q·2³²⌉).
const GS_LIMIT = Int64(72251709809335)
function gs_ok(f, g; limit::Int64 = GS_LIMIT)
    sum(abs2, f) + sum(abs2, g) < 16823 || return false
    rt1 = vect_FFT!(fxr_of.(f)); rt2 = vect_FFT!(fxr_of.(g))
    rt3 = vect_invnorm_fft(rt1, rt2, 0)
    vect_adj_fft!(rt1); vect_adj_fft!(rt2)
    rt1 .= fxr_mul.(rt1, fxr_of(q)); rt2 .= fxr_mul.(rt2, fxr_of(q))
    vect_mul_autoadj_fft!(rt1, rt3); vect_mul_autoadj_fft!(rt2, rt3)
    vect_iFFT!(rt1); vect_iFFT!(rt2)
    sn = sum(fxr_sqr(rt1[u]) + fxr_sqr(rt2[u]) for u in eachindex(rt1))
    sn < limit
end

# ── Babai reduction in fixed point (ePrint 2023/290 §2) ──────────────────────
"""
    reduce!(f, g, F, G) -> (F, G)

(F, G) −= k·(f, g) with k = round((F f* + G g*)/(f f* + g g*)), repeated while it shrinks (F, G). (f, g) are scaled so their largest coefficient is below 2^(15 − log₂n) (keeps
|FFT|² within the 31 integer bits of fxr); (F, G) are scaled by a further 2^u so that k carries
about R = reduce_bits bits (ntrugen profiles: 10, 9, 8 for n = 256, 512, 1024), and k is applied as k·2^u. f* and g* carry an extra 2^e for precision, removed before rounding.
"""
function reduce!(f, g, F, G)
    n = length(f); logn = trailing_zeros(n); e = max(15 - logn, 0)
    R = logn <= 8 ? 10 : 18 - logn                                    # reduce_bits, SOLVE_Falcon_* profiles
    fb = max(NG.maxbits(f), NG.maxbits(g))
    sf = max(fb - e, 0)
    rf = vect_FFT!(to_fxr.(f, sf)); rg = vect_FFT!(to_fxr.(g, sf))
    d = vect_invnorm_fft(rf, rg, e)
    vect_adj_fft!(rf); vect_adj_fft!(rg)
    vect_mul_autoadj_fft!(rf, d); vect_mul_autoadj_fft!(rg, d)          # 2^e·f*/(ff*+gg*), 2^e·g*/(…)
    # Scaled steps while they shrink (F, G); then unscaled steps until k = 0 or no progress.
    last = typemax(Int); last0 = typemax(Int); unscaled = false
    while true
        sz = max(NG.maxbits(F), NG.maxbits(G))
        unscaled |= sz >= last
        u = unscaled ? 0 : max(sz - fb - R, 0)
        if u == 0
            sz >= last0 && break
            last0 = sz
        end
        last = sz
        rF = vect_FFT!(to_fxr.(F, sf + u)); rG = vect_FFT!(to_fxr.(G, sf + u))
        vect_mul_fft!(rF, rf); vect_mul_fft!(rG, rg)
        rk = vect_iFFT!(rF .+ rG)
        k = BigInt[fxr_round(fxr_div2e(x, e)) for x in rk]
        if all(iszero, k)
            u == 0 && break
            unscaled = true; continue
        end
        F .-= NG.negamul(k, f) .<< u; G .-= NG.negamul(k, g) .<< u
    end
    F, G
end

# Recursive solver of FalconNTRUGen with the fixed-point reduction at every level.
function ntru_solve(f, g)
    n = length(f)
    n == 1 && return NG.ntru_solve(f, g)
    Fp, Gp = ntru_solve(NG.field_norm(f), NG.field_norm(g))
    F = NG.negamul(NG.lift(Fp), NG.galois_conjugate(g))
    G = NG.negamul(NG.lift(Gp), NG.galois_conjugate(f))
    reduce!(f, g, F, G)
end

end # module
