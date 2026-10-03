"""
X25519 (RFC 7748 §5): Diffie–Hellman on Curve25519, y² = x³ + 486662x² + x over GF(2²⁵⁵ − 19),
by the Montgomery ladder on u-coordinates. Two field engines: radix-2⁵¹ limbs (`x25519`) and a
BigInt transcription of RFC 7748 §5 (`x25519_ref`), cross-checked in the tests.
"""
module X25519

using Random

export x25519, x25519_base, x25519_keypair

const P = big(2)^255 - 19
const A24 = 121665
const BASE = vcat(UInt8(9), zeros(UInt8, 31))
const M51 = (UInt64(1) << 51) - 1
const LOW63 = 0x7fffffffffffffff

# ── GF(2²⁵⁵ − 19) in five 51-bit limbs; every operation returns carried limbs ──
const Fe = NTuple{5,UInt64}
const TWO_P = (UInt64(2) * (M51 - 18), UInt64(2) * M51, UInt64(2) * M51, UInt64(2) * M51, UInt64(2) * M51)

function carry(a::NTuple{5,UInt128})
    r0 = a[1]; r1 = a[2] + (r0 >> 51); r0 &= M51
    r2 = a[3] + (r1 >> 51); r1 &= M51
    r3 = a[4] + (r2 >> 51); r2 &= M51
    r4 = a[5] + (r3 >> 51); r3 &= M51
    r0 += 19 * (r4 >> 51); r4 &= M51
    r1 += r0 >> 51; r0 &= M51
    (r0 % UInt64, r1 % UInt64, r2 % UInt64, r3 % UInt64, r4 % UInt64)
end
carry(a::Fe) = carry(UInt128.(a))

fe_add(a::Fe, b::Fe) = carry(ntuple(i -> a[i] + b[i], 5))
fe_sub(a::Fe, b::Fe) = carry(ntuple(i -> a[i] + TWO_P[i] - b[i], 5))
fe_small(a::Fe, k::UInt64) = carry(ntuple(i -> UInt128(a[i]) * k, 5))

function fe_mul(a::Fe, b::Fe)
    a0, a1, a2, a3, a4 = UInt128.(a)
    b0, b1, b2, b3, b4 = UInt128.(b)
    s1, s2, s3, s4 = 19b1, 19b2, 19b3, 19b4
    carry((a0*b0 + a1*s4 + a2*s3 + a3*s2 + a4*s1,
           a0*b1 + a1*b0 + a2*s4 + a3*s3 + a4*s2,
           a0*b2 + a1*b1 + a2*b0 + a3*s4 + a4*s3,
           a0*b3 + a1*b2 + a2*b1 + a3*b0 + a4*s4,
           a0*b4 + a1*b3 + a2*b2 + a3*b1 + a4*b0))
end
fe_sq(a::Fe) = fe_mul(a, a)

const PM2_BITS = reverse(digits(Bool, P - 2, base=2))   # public exponent, MSB first

function fe_inv(a::Fe)
    r = (UInt64(1), UInt64(0), UInt64(0), UInt64(0), UInt64(0))
    for b in PM2_BITS
        r = fe_sq(r)
        b && (r = fe_mul(r, a))
    end
    r
end

# Four little-endian 64-bit words of the limb value (limbs < 2⁵², so it fits).
function fe_words(a::Fe)
    acc = UInt128(a[1]) + (UInt128(a[2]) << 51)
    w0 = acc % UInt64; acc >>= 64
    acc += UInt128(a[3]) << 38
    w1 = acc % UInt64; acc >>= 64
    acc += UInt128(a[4]) << 25
    w2 = acc % UInt64; acc >>= 64
    acc += UInt128(a[5]) << 12
    w3 = acc % UInt64
    (w0, w1, w2, w3), (acc >> 64) % UInt64
end

# w + 19k over four words; the caller keeps the sum below 2²⁵⁶.
function add19(w::NTuple{4,UInt64}, k::UInt64)
    c = UInt128(w[1]) + 19 * UInt128(k); r0 = c % UInt64
    c = UInt128(w[2]) + (c >> 64);       r1 = c % UInt64
    c = UInt128(w[3]) + (c >> 64);       r2 = c % UInt64
    c = UInt128(w[4]) + (c >> 64);       r3 = c % UInt64
    (r0, r1, r2, r3)
end

# Canonical 32-byte encoding, branch-free: fold 2²⁵⁵ ≡ 19 twice, then subtract p if v ≥ p.
# Single assignments: a reassigned captured variable gets boxed.
function fe_tobytes(a::Fe)
    w0, over = fe_words(a)
    w1 = add19((w0[1], w0[2], w0[3], w0[4] & LOW63), (w0[4] >> 63) | (over << 1))
    w2 = add19((w1[1], w1[2], w1[3], w1[4] & LOW63), w1[4] >> 63)
    c = add19(w2, UInt64(1))                 # v + 19 ≥ 2²⁵⁵ ⇔ v ≥ p
    m = -(c[4] >> 63)
    r = map((ci, wi) -> (m & ci) | (~m & wi), c, w2)
    out = Vector{UInt8}(undef, 32)
    for i in 0:31
        out[i+1] = (r[(i >> 3) + 1] >> (8 * (i & 7))) % UInt8
    end
    out[32] &= 0x7f
    out
end

# Little-endian 32 bytes; the top bit is masked (RFC 7748 §5); values in [p, 2²⁵⁵) are accepted.
function fe_frombytes(b::AbstractVector{UInt8})
    w = [reduce(|, UInt64(b[8i+j+1]) << (8j) for j in 0:7) for i in 0:3]
    w[4] &= 0x7fffffffffffffff
    (w[1] & M51, ((w[1] >> 51) | (w[2] << 13)) & M51, ((w[2] >> 38) | (w[3] << 26)) & M51,
     ((w[3] >> 25) | (w[4] << 39)) & M51, w[4] >> 12)
end
tobytes(x::BigInt) = UInt8[(x >> (8i)) & 0xff for i in 0:31]

cswap(s::UInt64, a::Fe, b::Fe) = (m = -s; d = ntuple(i -> m & (a[i] ⊻ b[i]), 5);
                                   (ntuple(i -> a[i] ⊻ d[i], 5), ntuple(i -> b[i] ⊻ d[i], 5)))

# decodeScalar25519 (RFC 7748 §5) as four little-endian 64-bit words.
function clamp_words(k::AbstractVector{UInt8})
    w = ntuple(i -> reduce(|, UInt64(k[8(i-1)+j+1]) << (8j) for j in 0:7), 4)
    (w[1] & ~UInt64(7), w[2], w[3], (w[4] & 0x7fffffffffffffff) | 0x4000000000000000)
end

# decodeScalar25519 (RFC 7748 §5).
function clamp_scalar(k::AbstractVector{UInt8})
    kb = collect(k); kb[1] &= 0xf8; kb[32] &= 0x7f; kb[32] |= 0x40
    sum(big(kb[i+1]) << (8i) for i in 0:31)
end

check32(v, what) = length(v) == 32 || throw(ArgumentError("X25519 $what must be 32 bytes"))

"X25519(k, u) (RFC 7748 §5): the u-coordinate of [k]·u, as 32 little-endian bytes."
function x25519(k::AbstractVector{UInt8}, u::AbstractVector{UInt8})
    check32(k, "scalar"); check32(u, "u-coordinate")
    s = clamp_words(k)
    x1 = fe_frombytes(u)
    one = (UInt64(1), UInt64(0), UInt64(0), UInt64(0), UInt64(0)); zero = (UInt64(0), UInt64(0), UInt64(0), UInt64(0), UInt64(0))
    x2, z2, x3, z3 = one, zero, x1, one
    swap = UInt64(0)
    for t in 254:-1:0
        kt = (s[(t >> 6) + 1] >> (t & 63)) & 1
        swap ⊻= kt
        x2, x3 = cswap(swap, x2, x3); z2, z3 = cswap(swap, z2, z3)
        swap = kt
        A = fe_add(x2, z2); AA = fe_sq(A)
        B = fe_sub(x2, z2); BB = fe_sq(B)
        E = fe_sub(AA, BB)
        C = fe_add(x3, z3); D = fe_sub(x3, z3)
        DA = fe_mul(D, A); CB = fe_mul(C, B)
        x3 = fe_sq(fe_add(DA, CB))
        z3 = fe_mul(x1, fe_sq(fe_sub(DA, CB)))
        x2 = fe_mul(AA, BB)
        z2 = fe_mul(E, fe_add(AA, fe_small(E, UInt64(A24))))
    end
    x2, x3 = cswap(swap, x2, x3); z2, z3 = cswap(swap, z2, z3)
    fe_tobytes(fe_mul(x2, fe_inv(z2)))
end

"RFC 7748 §5 pseudocode over BigInt: the reference the limb engine is tested against."
function x25519_ref(k::AbstractVector{UInt8}, u::AbstractVector{UInt8})
    check32(k, "scalar"); check32(u, "u-coordinate")
    s = clamp_scalar(k)
    ub = collect(u); ub[32] &= 0x7f
    x1 = mod(sum(big(ub[i+1]) << (8i) for i in 0:31), P)
    x2, z2, x3, z3, swap = big(1), big(0), x1, big(1), 0
    for t in 254:-1:0
        kt = Int((s >> t) & 1)
        swap ⊻= kt
        swap == 1 && ((x2, x3) = (x3, x2); (z2, z3) = (z3, z2))
        swap = kt
        A = x2 + z2; AA = mod(A^2, P); B = x2 - z2; BB = mod(B^2, P); E = AA - BB
        C = x3 + z3; D = x3 - z3; DA = mod(D * A, P); CB = mod(C * B, P)
        x3 = mod((DA + CB)^2, P); z3 = mod(x1 * (DA - CB)^2, P)
        x2 = mod(AA * BB, P); z2 = mod(E * (AA + A24 * E), P)
    end
    swap == 1 && ((x2, x3) = (x3, x2); (z2, z3) = (z3, z2))
    tobytes(mod(x2 * powermod(z2, P - 2, P), P))
end

"X25519(k, 9): the public key for scalar k."
x25519_base(k::AbstractVector{UInt8}) = x25519(k, BASE)

"Random X25519 key pair `(pk, sk)` from the OS CSPRNG."
function x25519_keypair()
    sk = rand(RandomDevice(), UInt8, 32)
    x25519_base(sk), sk
end

end # module
