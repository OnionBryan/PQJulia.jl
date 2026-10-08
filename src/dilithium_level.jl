function decompose(a::Int32)
    a1 = (a + 127) >> 7
    if GAMMA2 == Int32(div(Q - 1, 32))
        a1 = (a1 * 1025 + (Int32(1) << 21)) >> 22
        a1 &= 15
    elseif GAMMA2 == Int32(div(Q - 1, 88))
        a1 = (a1 * 11275 + (Int32(1) << 23)) >> 24
        a1 = xor(a1, ((43 - a1) >> 31) & a1)
    end
    a0 = a - a1 * 2 * GAMMA2
    a0 -= ((div(Q - 1, 2) - a0) >> 31) & Q
    return a1, a0
end

# Bitwise, not ||: x86 compiles the short-circuit form to branches.
function make_hint(a0::Int32, a1::Int32)::Bool
    return (a0 > GAMMA2) | (a0 < -GAMMA2) | ((a0 == -GAMMA2) & (a1 != 0))
end

function use_hint(a::Int32, hint::Bool)::Int32
    a1, a0 = decompose(a)
    !hint && return a1
    if GAMMA2 == Int32(div(Q - 1, 32))
        return a0 > 0 ? (a1 + 1) & 15 : (a1 - 1) & 15
    elseif GAMMA2 == Int32(div(Q - 1, 88))
        return a0 > 0 ? (a1 == 43 ? Int32(0) : a1 + 1) : (a1 == 0 ? Int32(43) : a1 - 1)
    end
end

# Rejection sampling, C ref rej_eta.
function rej_eta!(a, offset, len, buf, buflen)
    checkbounds(a, offset+1:offset+len); checkbounds(buf, 1:buflen)
    ctr = 0; pos = 1
    @inbounds while ctr < len && pos <= buflen
        t0 = UInt32(buf[pos]) & 0x0F
        t1 = UInt32(buf[pos]) >> 4
        pos += 1
        if ETA == 2
            if t0 < 15
                t0 = t0 - ((205*t0) >> 10)*5
                a[offset + ctr + 1] = Int32(2) - Int32(t0)
                ctr += 1
            end
            if t1 < 15 && ctr < len
                t1 = t1 - ((205*t1) >> 10)*5
                a[offset + ctr + 1] = Int32(2) - Int32(t1)
                ctr += 1
            end
        elseif ETA == 4
            if t0 < 9
                a[offset + ctr + 1] = Int32(4) - Int32(t0)
                ctr += 1
            end
            if t1 < 9 && ctr < len
                a[offset + ctr + 1] = Int32(4) - Int32(t1)
                ctr += 1
            end
        end
    end
    return ctr
end

function poly_uniform_eta!(a::Vector{Int32}, seed::Vector{UInt8}, nonce::UInt16)
    # SHAKE256(seed || nonce_le16) with re-squeeze loop (C ref: poly.c:435-453)
    # SHAKE256_RATE = 136; initial blocks: eta=2 → 1 block (136B), eta=4 → 2 blocks (272B)
    input = vcat(seed, UInt8[nonce & 0xff, (nonce >> 8) & 0xff])
    nblocks_init = ETA == 2 ? 1 : 2
    total_out = nblocks_init * 136
    buf = Keccak.shake256(input, UInt64(total_out))
    ctr = rej_eta!(a, 0, N, buf, total_out)

    # Re-squeeze loop: one SHAKE256 block at a time
    while ctr < N
        new_total = total_out + 136
        full = Keccak.shake256(input, UInt64(new_total))
        blk = full[total_out+1:new_total]
        total_out = new_total
        ctr += rej_eta!(a, ctr, N - ctr, blk, 136)
        wipe!(full, blk)
    end
    wipe!(input, buf)

    return a
end

function poly_uniform_gamma1!(a::Vector{Int32}, seed::Vector{UInt8}, nonce::UInt16)
    input = vcat(seed, UInt8[nonce & 0xff, (nonce >> 8) & 0xff])
    buf = Keccak.shake256(input, UInt64(POLYZ_PACKED))
    polyz_unpack!(a, buf)
    wipe!(input, buf)
    return a
end

function poly_challenge!(c::Vector{Int32}, seed::Vector{UInt8})
    # SampleInBall (C ref: poly.c:487-519). Uses incremental SHAKE256 squeeze.
    fill!(c, Int32(0))
    total_out = 136  # one SHAKE256 block
    buf = Keccak.shake256(seed, UInt64(total_out))
    signs = UInt64(0)
    for i in 0:7
        signs |= UInt64(buf[i+1]) << (8*i)
    end
    pos = 9
    for i in (N - TAU):(N - 1)
        local b
        while true
            # Re-squeeze if buffer exhausted (C ref: poly.c:509-512)
            if pos > length(buf)
                new_total = total_out + 136
                full = Keccak.shake256(seed, UInt64(new_total))
                buf = vcat(buf, full[total_out+1:new_total])
                total_out = new_total
            end
            b = Int(buf[pos]); pos += 1
            b <= i && break
        end
        c[i+1] = c[b+1]
        c[b+1] = Int32(1 - 2 * Int32(signs & 1))
        signs >>= 1
    end
    return c
end

# ==================== PACKING ====================

function polyt1_pack(a::Vector{Int32})::Vector{UInt8}
    checkbounds(a, 1:N)
    r = zeros(UInt8, POLYT1_PACKED)
    @inbounds for i in 0:(N÷4 - 1)
        r[5i+1] = (a[4i+1]) % UInt8
        r[5i+2] = ((a[4i+1] >> 8) | (a[4i+2] << 2)) % UInt8
        r[5i+3] = ((a[4i+2] >> 6) | (a[4i+3] << 4)) % UInt8
        r[5i+4] = ((a[4i+3] >> 4) | (a[4i+4] << 6)) % UInt8
        r[5i+5] = (a[4i+4] >> 2) % UInt8
    end
    return r
end

function polyt1_unpack!(r::Vector{Int32}, a::Vector{UInt8})
    checkbounds(r, 1:N); checkbounds(a, 1:POLYT1_PACKED)
    @inbounds for i in 0:(N÷4 - 1)
        r[4i+1] = Int32((UInt32(a[5i+1]) | (UInt32(a[5i+2]) << 8)) & 0x3FF)
        r[4i+2] = Int32(((UInt32(a[5i+2]) >> 2) | (UInt32(a[5i+3]) << 6)) & 0x3FF)
        r[4i+3] = Int32(((UInt32(a[5i+3]) >> 4) | (UInt32(a[5i+4]) << 4)) & 0x3FF)
        r[4i+4] = Int32(((UInt32(a[5i+4]) >> 6) | (UInt32(a[5i+5]) << 2)) & 0x3FF)
    end
    return r
end

function polyeta_pack(a::Vector{Int32})::Vector{UInt8}
    checkbounds(a, 1:N)
    r = zeros(UInt8, POLYETA_PACKED)
    if ETA == 4
        @inbounds for i in 0:(N÷2 - 1)
            t0 = (ETA - a[2i+1]) % UInt8
            t1 = (ETA - a[2i+2]) % UInt8
            r[i+1] = t0 | (t1 << 4)
        end
    elseif ETA == 2
        @inbounds for i in 0:(N÷8 - 1)
            t = ntuple(j -> (ETA - a[8i+j]) % UInt8, 8)
            r[3i+1] = t[1] | (t[2] << 3) | (t[3] << 6)
            r[3i+2] = (t[3] >> 2) | (t[4] << 1) | (t[5] << 4) | (t[6] << 7)
            r[3i+3] = (t[6] >> 1) | (t[7] << 2) | (t[8] << 5)
        end
    end
    return r
end

function polyeta_unpack!(r::Vector{Int32}, a::AbstractVector{UInt8})
    checkbounds(r, 1:N); checkbounds(a, 1:POLYETA_PACKED)
    if ETA == 4
        @inbounds for i in 0:(N÷2 - 1)
            r[2i+1] = Int32(ETA) - Int32(a[i+1] & 0x0F)
            r[2i+2] = Int32(ETA) - Int32(a[i+1] >> 4)
        end
    elseif ETA == 2
        @inbounds for i in 0:(N÷8 - 1)
            r[8i+1] = Int32(a[3i+1]) & 7
            r[8i+2] = (Int32(a[3i+1]) >> 3) & 7
            r[8i+3] = ((Int32(a[3i+1]) >> 6) | (Int32(a[3i+2]) << 2)) & 7
            r[8i+4] = (Int32(a[3i+2]) >> 1) & 7
            r[8i+5] = (Int32(a[3i+2]) >> 4) & 7
            r[8i+6] = ((Int32(a[3i+2]) >> 7) | (Int32(a[3i+3]) << 1)) & 7
            r[8i+7] = (Int32(a[3i+3]) >> 2) & 7
            r[8i+8] = (Int32(a[3i+3]) >> 5) & 7
            for j in 1:8
                r[8i+j] = Int32(ETA) - r[8i+j]
            end
        end
    end
    return r
end

function polyz_pack(a::Vector{Int32})::Vector{UInt8}
    checkbounds(a, 1:N)
    r = zeros(UInt8, POLYZ_PACKED)
    if GAMMA1 == Int32(1 << 19)
        # 20-bit packing: 2 coefficients -> 5 bytes
        @inbounds for i in 0:(N÷2 - 1)
            t0 = UInt32(GAMMA1 - a[2i+1])
            t1 = UInt32(GAMMA1 - a[2i+2])
            r[5i+1] = (t0) % UInt8
            r[5i+2] = (t0 >> 8) % UInt8
            r[5i+3] = ((t0 >> 16) | (t1 << 4)) % UInt8
            r[5i+4] = (t1 >> 4) % UInt8
            r[5i+5] = (t1 >> 12) % UInt8
        end
    elseif GAMMA1 == Int32(1 << 17)
        # 18-bit packing: 4 coefficients -> 9 bytes
        @inbounds for i in 0:(N÷4 - 1)
            t0 = UInt32(GAMMA1 - a[4i+1])
            t1 = UInt32(GAMMA1 - a[4i+2])
            t2 = UInt32(GAMMA1 - a[4i+3])
            t3 = UInt32(GAMMA1 - a[4i+4])
            r[9i+1] = (t0) % UInt8
            r[9i+2] = (t0 >> 8) % UInt8
            r[9i+3] = ((t0 >> 16) | (t1 << 2)) % UInt8
            r[9i+4] = (t1 >> 6) % UInt8
            r[9i+5] = ((t1 >> 14) | (t2 << 4)) % UInt8
            r[9i+6] = (t2 >> 4) % UInt8
            r[9i+7] = ((t2 >> 12) | (t3 << 6)) % UInt8
            r[9i+8] = (t3 >> 2) % UInt8
            r[9i+9] = (t3 >> 10) % UInt8
        end
    end
    return r
end
function polyz_unpack!(r::Vector{Int32}, a::Vector{UInt8})
    checkbounds(r, 1:N); checkbounds(a, 1:POLYZ_PACKED)
    if GAMMA1 == Int32(1 << 19)
        @inbounds for i in 0:(N÷2 - 1)
            r[2i+1] = Int32(UInt32(a[5i+1]) | (UInt32(a[5i+2]) << 8) | (UInt32(a[5i+3]) << 16)) & Int32(0xFFFFF)
            r[2i+2] = Int32((UInt32(a[5i+3]) >> 4) | (UInt32(a[5i+4]) << 4) | (UInt32(a[5i+5]) << 12)) & Int32(0xFFFFF)
            r[2i+1] = GAMMA1 - r[2i+1]
            r[2i+2] = GAMMA1 - r[2i+2]
        end
    elseif GAMMA1 == Int32(1 << 17)
        @inbounds for i in 0:(N÷4 - 1)
            r[4i+1] = Int32(UInt32(a[9i+1]) | (UInt32(a[9i+2]) << 8) | (UInt32(a[9i+3]) << 16)) & Int32(0x3FFFF)
            r[4i+2] = Int32((UInt32(a[9i+3]) >> 2) | (UInt32(a[9i+4]) << 6) | (UInt32(a[9i+5]) << 14)) & Int32(0x3FFFF)
            r[4i+3] = Int32((UInt32(a[9i+5]) >> 4) | (UInt32(a[9i+6]) << 4) | (UInt32(a[9i+7]) << 12)) & Int32(0x3FFFF)
            r[4i+4] = Int32((UInt32(a[9i+7]) >> 6) | (UInt32(a[9i+8]) << 2) | (UInt32(a[9i+9]) << 10)) & Int32(0x3FFFF)
            for j in 1:4; r[4i+j] = GAMMA1 - r[4i+j]; end
        end
    end
    return r
end
function polyw1_pack(a::Vector{Int32})::Vector{UInt8}
    checkbounds(a, 1:N)
    r = zeros(UInt8, POLYW1_PACKED)
    if GAMMA2 == Int32(div(Q - 1, 32))
        # 4-bit: 2 coefficients per byte (w1 range 0..15)
        @inbounds for i in 0:(N÷2 - 1)
            r[i+1] = (a[2i+1] | (a[2i+2] << 4)) % UInt8
        end
    elseif GAMMA2 == Int32(div(Q - 1, 88))
        # 6-bit: 4 coefficients per 3 bytes (w1 range 0..43)
        @inbounds for i in 0:(N÷4 - 1)
            r[3i+1] = (a[4i+1] | (a[4i+2] << 6)) % UInt8
            r[3i+2] = ((a[4i+2] >> 2) | (a[4i+3] << 4)) % UInt8
            r[3i+3] = ((a[4i+3] >> 4) | (a[4i+4] << 2)) % UInt8
        end
    end
    return r
end
function polyt0_pack(a::Vector{Int32})::Vector{UInt8}
    # 13-bit packing: 8 coefficients → 13 bytes. C ref: packing.c:664-702
    checkbounds(a, 1:N)
    r = zeros(UInt8, POLYT0_PACKED)
    @inbounds for i in 0:(N÷8 - 1)
        ts = ntuple(k -> UInt32((1 << (D-1)) - a[8i+k]), 8)
        r[13i+1]  = (ts[1]) % UInt8
        r[13i+2]  = ((ts[1] >> 8) | (ts[2] << 5)) % UInt8
        r[13i+3]  = (ts[2] >> 3) % UInt8
        r[13i+4]  = ((ts[2] >> 11) | (ts[3] << 2)) % UInt8
        r[13i+5]  = ((ts[3] >> 6) | (ts[4] << 7)) % UInt8
        r[13i+6]  = (ts[4] >> 1) % UInt8
        r[13i+7]  = ((ts[4] >> 9) | (ts[5] << 4)) % UInt8
        r[13i+8]  = (ts[5] >> 4) % UInt8
        r[13i+9]  = ((ts[5] >> 12) | (ts[6] << 1)) % UInt8
        r[13i+10] = ((ts[6] >> 7) | (ts[7] << 6)) % UInt8
        r[13i+11] = (ts[7] >> 2) % UInt8
        r[13i+12] = ((ts[7] >> 10) | (ts[8] << 3)) % UInt8
        r[13i+13] = (ts[8] >> 5) % UInt8
    end
    return r
end
function polyt0_unpack!(r::Vector{Int32}, a::AbstractVector{UInt8})
    # 13-bit unpacking: 13 bytes → 8 coefficients. C ref: packing.c:712-763
    checkbounds(r, 1:N); checkbounds(a, 1:POLYT0_PACKED)
    @inbounds for i in 0:(N÷8 - 1)
        r[8i+1] = Int32(UInt32(a[13i+1]) | (UInt32(a[13i+2]) << 8)) & Int32(0x1FFF)
        r[8i+2] = Int32((UInt32(a[13i+2]) >> 5) | (UInt32(a[13i+3]) << 3) | (UInt32(a[13i+4]) << 11)) & Int32(0x1FFF)
        r[8i+3] = Int32((UInt32(a[13i+4]) >> 2) | (UInt32(a[13i+5]) << 6)) & Int32(0x1FFF)
        r[8i+4] = Int32((UInt32(a[13i+5]) >> 7) | (UInt32(a[13i+6]) << 1) | (UInt32(a[13i+7]) << 9)) & Int32(0x1FFF)
        r[8i+5] = Int32((UInt32(a[13i+7]) >> 4) | (UInt32(a[13i+8]) << 4) | (UInt32(a[13i+9]) << 12)) & Int32(0x1FFF)
        r[8i+6] = Int32((UInt32(a[13i+9]) >> 1) | (UInt32(a[13i+10]) << 7)) & Int32(0x1FFF)
        r[8i+7] = Int32((UInt32(a[13i+10]) >> 6) | (UInt32(a[13i+11]) << 2) | (UInt32(a[13i+12]) << 10)) & Int32(0x1FFF)
        r[8i+8] = Int32((UInt32(a[13i+12]) >> 3) | (UInt32(a[13i+13]) << 5)) & Int32(0x1FFF)
        for k in 1:8; r[8i+k] = Int32((1 << (D-1))) - r[8i+k]; end
    end
    return r
end
# ==================== KEY GENERATION ====================

function dilithium_keygen_derand(xi::Vector{UInt8})
    length(xi) == SEEDBYTES || throw(ArgumentError("keygen seed must be $SEEDBYTES bytes"))
    seed = vcat(xi, UInt8[K, L])
    expanded = Keccak.shake256(seed, UInt64(2*SEEDBYTES + CRHBYTES))
    rho = expanded[1:SEEDBYTES]
    rhoprime = expanded[SEEDBYTES+1:SEEDBYTES+CRHBYTES]
    key = expanded[SEEDBYTES+CRHBYTES+1:2*SEEDBYTES+CRHBYTES]

    # Expand A
    A = [zeros(Int32, N) for _ in 1:K, _ in 1:L]
    for i in 1:K, j in 1:L
        poly_uniform!(A[i,j], rho, UInt16((i-1) << 8 | (j-1)))
    end

    # Sample s1, s2
    s1 = [zeros(Int32, N) for _ in 1:L]
    for i in 1:L
        poly_uniform_eta!(s1[i], rhoprime, UInt16(i-1))
    end
    s2 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        poly_uniform_eta!(s2[i], rhoprime, UInt16(L + i - 1))
    end

    # t = As1 + s2
    s1hat = [copy(s) for s in s1]
    for i in 1:L; ntt!(s1hat[i]); end

    # One row of t at a time, split by power2round: t = t1*2^D + t0
    t = zeros(Int32, N)
    t1 = [zeros(Int32, N) for _ in 1:K]
    t0 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        fill!(t, Int32(0))
        for j in 1:L
            poly_pointwise_acc!(t, A[i,j], s1hat[j])
        end
        poly_reduce!(t)
        invntt!(t)
        poly_add!(t, t, s2[i])
        poly_caddq!(t)
        for j in 1:N
            t1[i][j], t0[i][j] = power2round(t[j])
        end
    end

    # Pack pk = rho || t1
    pk = copy(rho)
    for i in 1:K; append!(pk, polyt1_pack(t1[i])); end

    # tr = H(pk)
    tr = Keccak.shake256(pk, UInt64(TRBYTES))

    # Pack sk = rho || key || tr || s1 || s2 || t0
    sk = append!(sizehint!(UInt8[], SK_BYTES), rho, key, tr)   # no regrowth copies
    for i in 1:L
        e = polyeta_pack(s1[i]); append!(sk, e); wipe!(e)
    end
    for i in 1:K
        e = polyeta_pack(s2[i]); append!(sk, e); wipe!(e)
    end
    for i in 1:K
        e = polyt0_pack(t0[i]); append!(sk, e); wipe!(e)
    end
    wipe!(seed, expanded, rhoprime, key, s1, s2, s1hat, t, t0)

    return pk, sk
end

function dilithium_keygen()
    xi = rand(RandomDevice(), UInt8, SEEDBYTES)
    kp = dilithium_keygen_derand(xi)
    wipe!(xi)
    return kp
end

# ==================== SIGN ====================

function unpack_sk(sk::Vector{UInt8})
    pos = 1
    rho = sk[pos:pos+SEEDBYTES-1]; pos += SEEDBYTES
    key = sk[pos:pos+SEEDBYTES-1]; pos += SEEDBYTES
    tr = sk[pos:pos+TRBYTES-1]; pos += TRBYTES

    s1 = [zeros(Int32, N) for _ in 1:L]
    for i in 1:L
        polyeta_unpack!(s1[i], view(sk, pos:pos+POLYETA_PACKED-1)); pos += POLYETA_PACKED
    end
    s2 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        polyeta_unpack!(s2[i], view(sk, pos:pos+POLYETA_PACKED-1)); pos += POLYETA_PACKED
    end
    t0 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        polyt0_unpack!(t0[i], view(sk, pos:pos+POLYT0_PACKED-1)); pos += POLYT0_PACKED
    end
    return rho, key, tr, s1, s2, t0
end

function expand_A(rho::Vector{UInt8})
    A = [zeros(Int32, N) for _ in 1:K, _ in 1:L]
    for i in 1:K, j in 1:L
        poly_uniform!(A[i,j], rho, UInt16((i-1) << 8 | (j-1)))
    end
    return A
end

function sample_y!(y::Vector{Vector{Int32}}, rhoprime::Vector{UInt8}, nonce::Int)
    for i in 1:L
        poly_uniform_gamma1!(y[i], rhoprime, (L * nonce + i - 1) % UInt16)
    end
end

function compute_w!(w1::Vector{Vector{Int32}}, w0::Vector{Vector{Int32}}, A::Matrix{Vector{Int32}}, y::Vector{Vector{Int32}}, tmp::Vector{Int32})
    zy = [copy(y[i]) for i in 1:L]
    for i in 1:L; ntt!(zy[i]); end
    for i in 1:K
        fill!(w1[i], Int32(0))
        for j in 1:L
            poly_pointwise_acc!(w1[i], A[i,j], zy[j])
        end
        poly_reduce!(w1[i])
        invntt!(w1[i])
        poly_caddq!(w1[i])
    end

    # Decompose w
    for i in 1:K
        for j in 1:N
            w1[i][j], w0[i][j] = decompose(w1[i][j])
        end
    end
    wipe!(zy)
end

function compute_challenge(mu::Vector{UInt8}, w1::Vector{Vector{Int32}}, cp::Vector{Int32})
    w1_packed = UInt8[]
    for i in 1:K; append!(w1_packed, polyw1_pack(w1[i])); end
    c_tilde = Keccak.shake256(vcat(mu, w1_packed), UInt64(CTILDEBYTES))
    poly_challenge!(cp, c_tilde)
    cp_hat = copy(cp); ntt!(cp_hat)
    return c_tilde, cp_hat
end

function compute_z_and_check_norm!(z::Vector{Vector{Int32}}, cp_hat::Vector{Int32}, s1::Vector{Vector{Int32}}, y::Vector{Vector{Int32}})
    for i in 1:L
        poly_pointwise!(z[i], cp_hat, s1[i])
        invntt!(z[i])
        poly_add!(z[i], z[i], y[i])
        poly_reduce!(z[i])
    end

    for i in 1:L
        if poly_chknorm(z[i], GAMMA1 - BETA)
            return true
        end
    end
    return false
end

function compute_w0_and_check_norm!(w0::Vector{Vector{Int32}}, cp_hat::Vector{Int32}, s2::Vector{Vector{Int32}}, tmp::Vector{Int32})
    for i in 1:K
        poly_pointwise!(tmp, cp_hat, s2[i])
        invntt!(tmp)
        poly_sub!(w0[i], w0[i], tmp)
        poly_reduce!(w0[i])
    end
    for i in 1:K
        if poly_chknorm(w0[i], GAMMA2 - BETA)
            return true
        end
    end
    return false
end

function make_hints_and_check!(h::Vector{Vector{Int32}}, w0::Vector{Vector{Int32}}, w1::Vector{Vector{Int32}}, cp_hat::Vector{Int32}, t0::Vector{Vector{Int32}})
    # ct0
    for i in 1:K
        poly_pointwise!(h[i], cp_hat, t0[i])
        invntt!(h[i])
        poly_reduce!(h[i])
    end
    for i in 1:K
        if poly_chknorm(h[i], GAMMA2)
            return true
        end
    end

    # Make hints
    for i in 1:K
        poly_add!(w0[i], w0[i], h[i])
    end
    hints_count = 0
    for i in 1:K
        for j in 1:N
            h[i][j] = Int32(make_hint(w0[i][j], w1[i][j]))
            hints_count += h[i][j]
        end
    end
    if hints_count > OMEGA
        return true
    end
    return false
end

function pack_signature(c_tilde::Vector{UInt8}, z::Vector{Vector{Int32}}, h::Vector{Vector{Int32}})
    sig = copy(c_tilde)
    for i in 1:L; append!(sig, polyz_pack(z[i])); end
    h_packed = zeros(UInt8, OMEGA + K)
    k_pos = 0
    for i in 1:K
        for j in 1:N
            if h[i][j] != 0
                h_packed[k_pos + 1] = (j - 1) % UInt8
                k_pos += 1
            end
        end
        h_packed[OMEGA + i] = (k_pos) % UInt8
    end
    append!(sig, h_packed)
    return sig
end

"""ML-DSA.Sign_internal core (FIPS 204 Alg. 7) from μ = H(tr ‖ M′). Every signing entry point lands here."""
function sign_mu(mu::Vector{UInt8}, sk::Vector{UInt8}, rnd::Vector{UInt8})
    length(mu) == CRHBYTES || error("mu must be $CRHBYTES bytes")
    length(rnd) == 32 || error("rnd must be 32 bytes")
    length(sk) == SK_BYTES || error("$IDENTIFIER secret key must be $SK_BYTES bytes")

    rho, key, tr, s1, s2, t0 = unpack_sk(sk)
    # skDecode range check: s1, s2 coefficients must lie in [-η, η] (Wycheproof InvalidPrivateKey)
    all(v -> all(x -> -ETA <= x <= ETA, v), s1) && all(v -> all(x -> -ETA <= x <= ETA, v), s2) ||
        throw(ArgumentError("$IDENTIFIER secret key has s1/s2 coefficients outside [-η, η]"))
    A = expand_A(rho)
    for i in 1:L; ntt!(s1[i]); end
    for i in 1:K; ntt!(s2[i]); end
    for i in 1:K; ntt!(t0[i]); end

    kin = vcat(key, rnd, mu)
    rhoprime = Keccak.shake256(kin, UInt64(CRHBYTES))

    nonce = 0  # Int, not UInt16 — avoids overflow at 9362 iterations for L=7 (pq-crystals/dilithium#110)
    y = [zeros(Int32, N) for _ in 1:L]
    z = [zeros(Int32, N) for _ in 1:L]
    w1 = [zeros(Int32, N) for _ in 1:K]
    w0 = [zeros(Int32, N) for _ in 1:K]
    h = [zeros(Int32, N) for _ in 1:K]
    cp = zeros(Int32, N)
    tmp = zeros(Int32, N)

    while true
        sample_y!(y, rhoprime, nonce)
        compute_w!(w1, w0, A, y, tmp)                  # w1 = HighBits(Ay), w0 = LowBits(Ay)
        c_tilde, cp_hat = compute_challenge(mu, w1, cp)
        if compute_z_and_check_norm!(z, cp_hat, s1, y)
            nonce += 1; continue
        end
        if compute_w0_and_check_norm!(w0, cp_hat, s2, tmp)
            nonce += 1; continue
        end
        if make_hints_and_check!(h, w0, w1, cp_hat, t0)
            nonce += 1; continue
        end
        sig = pack_signature(c_tilde, z, h)
        wipe!(key, s1, s2, t0, kin, rhoprime, y, z, w1, w0, h, tmp)
        return sig
    end
end

sk_tr(sk::Vector{UInt8}) = sk[2*SEEDBYTES+1:2*SEEDBYTES+TRBYTES]
mu_of(tr::Vector{UInt8}, mprime::Vector{UInt8}) = Keccak.shake256(vcat(tr, mprime), UInt64(CRHBYTES))

# M′ for pure ML-DSA (FIPS 204 Alg. 2/3): 0x00 ‖ |ctx| ‖ ctx ‖ M.
pure_mprime(msg, context) = vcat(UInt8[0x00, UInt8(length(context))], context, msg)

"""ML-DSA.Sign_internal (FIPS 204 Alg. 7): signs the formatted message M′ as given, no domain separator."""
function dilithium_sign_internal(mprime::Vector{UInt8}, sk::Vector{UInt8}, rnd::Vector{UInt8})
    length(sk) == SK_BYTES || error("$IDENTIFIER secret key must be $SK_BYTES bytes")
    return sign_mu(mu_of(sk_tr(sk), mprime), sk, rnd)
end

"""Sign with explicit μ (FIPS 204 external-μ interface)."""
dilithium_sign_internal_mu(mu::Vector{UInt8}, sk::Vector{UInt8}, rnd::Vector{UInt8}) = sign_mu(mu, sk, rnd)

"""Alias of `dilithium_sign_internal`."""
dilithium_sign_internal_msg(mprime::Vector{UInt8}, sk::Vector{UInt8}, rnd::Vector{UInt8}) =
    dilithium_sign_internal(mprime, sk, rnd)

"""ML-DSA.Sign (FIPS 204 Alg. 2) with caller-supplied 32-byte `rnd` (all zeros = deterministic variant)."""
function dilithium_sign_derand(msg::Vector{UInt8}, sk::Vector{UInt8}, rnd::Vector{UInt8}; context::Vector{UInt8}=UInt8[])
    length(context) > 255 && error("Context string must be ≤ 255 bytes (FIPS 204 §5.2)")
    return dilithium_sign_internal(pure_mprime(msg, context), sk, rnd)
end

"""ML-DSA.Sign (FIPS 204 Alg. 2). Hedged by default (rnd from the OS CSPRNG); `hedged=false` is the deterministic variant."""
function dilithium_sign(msg::Vector{UInt8}, sk::Vector{UInt8}; hedged::Bool=true, context::Vector{UInt8}=UInt8[])
    rnd = hedged ? rand(RandomDevice(), UInt8, 32) : zeros(UInt8, 32)
    sig = dilithium_sign_derand(msg, sk, rnd; context=context)
    wipe!(rnd)
    return sig
end

# ==================== VERIFY ====================

"""ML-DSA.Verify_internal core (FIPS 204 Alg. 8) from μ."""
function dilithium_verify_mu(mu::Vector{UInt8}, sig::Vector{UInt8}, pk::Vector{UInt8})
    length(mu) != CRHBYTES && return false
    length(sig) != SIG_BYTES && return false
    length(pk) != PK_BYTES && return false

    # Unpack pk
    rho = pk[1:SEEDBYTES]
    t1 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        polyt1_unpack!(t1[i], pk[SEEDBYTES + (i-1)*POLYT1_PACKED + 1 : SEEDBYTES + i*POLYT1_PACKED])
    end

    # Unpack sig
    c_tilde = sig[1:CTILDEBYTES]
    z = [zeros(Int32, N) for _ in 1:L]
    pos = CTILDEBYTES + 1
    for i in 1:L
        polyz_unpack!(z[i], sig[pos:pos+POLYZ_PACKED-1]); pos += POLYZ_PACKED
    end

    # Unpack h
    h = [zeros(Int32, N) for _ in 1:K]
    h_raw = sig[pos:end]
    k_pos = 0
    for i in 1:K
        limit = Int(h_raw[OMEGA + i])
        (limit < k_pos || limit > OMEGA) && return false
        for j in (k_pos+1):limit
            idx = Int(h_raw[j]) + 1
            (j > k_pos + 1 && h_raw[j] <= h_raw[j-1]) && return false
            h[i][idx] = Int32(1)
        end
        k_pos = limit
    end

    # Extra hint indices must be zero (C: for(j=k;j<OMEGA;++j) if(sig[j]) return 1)
    for j in (k_pos+1):OMEGA
        h_raw[j] != 0 && return false
    end

    # Check z norm
    for i in 1:L
        poly_chknorm(z[i], GAMMA1 - BETA) && return false
    end

    # w1' = Az - c*t1*2^D
    cp = zeros(Int32, N)
    poly_challenge!(cp, c_tilde)

    for i in 1:L; ntt!(z[i]); end
    w1p = [zeros(Int32, N) for _ in 1:K]
    tmp = zeros(Int32, N)
    # Each A[i,j] is used once, so expand it into tmp one entry at a time
    for i in 1:K
        fill!(w1p[i], Int32(0))
        for j in 1:L
            poly_uniform!(tmp, rho, UInt16((i-1) << 8 | (j-1)))
            poly_pointwise_acc!(w1p[i], tmp, z[j])
        end
    end

    ntt!(cp)
    for i in 1:K
        poly_shiftl!(t1[i])
        ntt!(t1[i])
        poly_pointwise!(tmp, cp, t1[i])
        poly_sub!(w1p[i], w1p[i], tmp)
        poly_reduce!(w1p[i])
        invntt!(w1p[i])
        poly_caddq!(w1p[i])
    end

    # Use hints to recover w1
    for i in 1:K
        for j in 1:N
            w1p[i][j] = use_hint(w1p[i][j], h[i][j] != 0)
        end
    end

    # Recompute challenge
    w1_packed = UInt8[]
    for i in 1:K; append!(w1_packed, polyw1_pack(w1p[i])); end
    c2 = Keccak.shake256(vcat(mu, w1_packed), UInt64(CTILDEBYTES))

    return c_tilde == c2
end

"""ML-DSA.Verify_internal (FIPS 204 Alg. 8) on the formatted message M′."""
function dilithium_verify_internal(mprime::Vector{UInt8}, sig::Vector{UInt8}, pk::Vector{UInt8})
    length(pk) != PK_BYTES && return false
    return dilithium_verify_mu(mu_of(Keccak.shake256(pk, UInt64(TRBYTES)), mprime), sig, pk)
end

"""ML-DSA.Verify (FIPS 204 Alg. 3). A context over 255 bytes is rejected (returns false)."""
function dilithium_verify(msg::Vector{UInt8}, sig::Vector{UInt8}, pk::Vector{UInt8}; context::Vector{UInt8}=UInt8[])
    length(context) > 255 && return false
    return dilithium_verify_internal(pure_mprime(msg, context), sig, pk)
end

# ==================== PREHASH (HashML-DSA) ====================

# OID table: DER-encoded OIDs for NIST hash algorithms
# All are 11 bytes: [0x06, 0x09, 0x60, 0x86, 0x48, 0x01, 0x65, 0x03, 0x04, 0x02, X]
const HASH_OIDS = Dict{String, Vector{UInt8}}(
    "SHA2-256"    => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x01],
    "SHA2-384"    => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x02],
    "SHA2-512"    => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x03],
    "SHA3-256"    => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x08],
    "SHA3-384"    => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x09],
    "SHA3-512"    => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x0A],
    "SHAKE-128"   => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x0B],
    "SHAKE-256"   => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x0C],
    "SHA2-224"    => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x04],
    "SHA3-224"    => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x07],
    "SHA2-512/224" => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x05],
    "SHA2-512/256" => UInt8[0x06,0x09,0x60,0x86,0x48,0x01,0x65,0x03,0x04,0x02,0x06],
)

function prehash_message(msg::Vector{UInt8}, hash_alg::String)::Vector{UInt8}
    if hash_alg == "SHA2-256"
        return SHA.sha256(msg)
    elseif hash_alg == "SHA2-384"
        return SHA.sha384(msg)
    elseif hash_alg == "SHA2-512"
        return SHA.sha512(msg)
    elseif hash_alg == "SHA3-256"
        return Keccak.sha3_256(msg)
    elseif hash_alg == "SHA3-384"
        return SHA.sha3_384(msg)
    elseif hash_alg == "SHA3-512"
        return Keccak.sha3_512(msg)
    elseif hash_alg == "SHAKE-128"
        return Keccak.shake128(msg, UInt64(32))  # 256 bits
    elseif hash_alg == "SHAKE-256"
        return Keccak.shake256(msg, UInt64(64))  # 512 bits
    elseif hash_alg == "SHA2-224"
        return SHA.sha224(msg)
    elseif hash_alg == "SHA3-224"
        return SHA.sha3_224(msg)
    elseif hash_alg == "SHA2-512/224"
        return SHA.sha2_512_224(msg)
    elseif hash_alg == "SHA2-512/256"
        return SHA.sha2_512_256(msg)
    else
        error("Unsupported hash algorithm: $hash_alg")
    end
end

# M′ for HashML-DSA (FIPS 204 Alg. 4/5): 0x01 ‖ |ctx| ‖ ctx ‖ OID ‖ PH(M).
function prehash_mprime(msg::Vector{UInt8}, hash_alg::String, context::Vector{UInt8})
    haskey(HASH_OIDS, hash_alg) || error("Unknown hash algorithm: $hash_alg")
    return vcat(UInt8[0x01, UInt8(length(context))], context, HASH_OIDS[hash_alg], prehash_message(msg, hash_alg))
end

"""HashML-DSA.Sign (FIPS 204 Alg. 4) with caller-supplied 32-byte `rnd`."""
function dilithium_sign_prehash_derand(msg::Vector{UInt8}, sk::Vector{UInt8}, hash_alg::String, rnd::Vector{UInt8};
                                       context::Vector{UInt8}=UInt8[])
    length(context) > 255 && error("Context string must be ≤ 255 bytes (FIPS 204 §5.2)")
    return dilithium_sign_internal(prehash_mprime(msg, hash_alg, context), sk, rnd)
end

"""HashML-DSA.Sign (FIPS 204 Alg. 4). Hedged by default; `hedged=false` is the deterministic variant."""
function dilithium_sign_prehash(msg::Vector{UInt8}, sk::Vector{UInt8}, hash_alg::String;
                                hedged::Bool=true, context::Vector{UInt8}=UInt8[])
    rnd = hedged ? rand(RandomDevice(), UInt8, 32) : zeros(UInt8, 32)
    sig = dilithium_sign_prehash_derand(msg, sk, hash_alg, rnd; context=context)
    wipe!(rnd)
    return sig
end

"""HashML-DSA.Verify (FIPS 204 Alg. 5). A context over 255 bytes is rejected (returns false)."""
function dilithium_verify_prehash(msg::Vector{UInt8}, sig::Vector{UInt8}, pk::Vector{UInt8},
                                  hash_alg::String; context::Vector{UInt8}=UInt8[])
    length(context) > 255 && return false
    return dilithium_verify_internal(prehash_mprime(msg, hash_alg, context), sig, pk)
end
