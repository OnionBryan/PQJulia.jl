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

function make_hint(a0::Int32, a1::Int32)::Bool
    return a0 > GAMMA2 || a0 < -GAMMA2 || (a0 == -GAMMA2 && a1 != 0)
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

function poly_uniform_eta!(a::Vector{Int32}, seed::Vector{UInt8}, nonce::UInt16)
    buf = SHA.shake256(vcat(seed, UInt8[nonce & 0xff, (nonce >> 8) & 0xff]), UInt64(512))
    ctr = 0; pos = 1
    while ctr < N && pos <= length(buf)
        t0 = UInt32(buf[pos]) & 0x0F
        t1 = UInt32(buf[pos]) >> 4
        pos += 1
        if ETA == 2
            if t0 < 15
                t0 = t0 - ((205*t0) >> 10)*5  # mod 5 trick
                ctr += 1; a[ctr] = Int32(2) - Int32(t0)
            end
            if t1 < 15 && ctr < N
                t1 = t1 - ((205*t1) >> 10)*5
                ctr += 1; a[ctr] = Int32(2) - Int32(t1)
            end
        elseif ETA == 4
            if t0 < 9
                ctr += 1; a[ctr] = Int32(4) - Int32(t0)
            end
            if t1 < 9 && ctr < N
                ctr += 1; a[ctr] = Int32(4) - Int32(t1)
            end
        end
    end
    ctr < N && error("poly_uniform_eta: insufficient samples ($ctr/$N)")
    return a
end

function poly_uniform_gamma1!(a::Vector{Int32}, seed::Vector{UInt8}, nonce::UInt16)
    buf = SHA.shake256(vcat(seed, UInt8[nonce & 0xff, (nonce >> 8) & 0xff]), UInt64(POLYZ_PACKED))
    polyz_unpack!(a, buf)
    return a
end

function poly_challenge!(c::Vector{Int32}, seed::Vector{UInt8})
    fill!(c, Int32(0))
    buf = SHA.shake256(seed, UInt64(256))
    signs = UInt64(0)
    for i in 0:7
        signs |= UInt64(buf[i+1]) << (8*i)
    end
    pos = 9
    for i in (N - TAU):(N - 1)
        local b
        while true
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
    r = zeros(UInt8, POLYT1_PACKED)
    for i in 0:(N÷4 - 1)
        r[5i+1] = (a[4i+1]) % UInt8
        r[5i+2] = ((a[4i+1] >> 8) | (a[4i+2] << 2)) % UInt8
        r[5i+3] = ((a[4i+2] >> 6) | (a[4i+3] << 4)) % UInt8
        r[5i+4] = ((a[4i+3] >> 4) | (a[4i+4] << 6)) % UInt8
        r[5i+5] = (a[4i+4] >> 2) % UInt8
    end
    return r
end

function polyt1_unpack!(r::Vector{Int32}, a::Vector{UInt8})
    for i in 0:(N÷4 - 1)
        r[4i+1] = Int32((UInt32(a[5i+1]) | (UInt32(a[5i+2]) << 8)) & 0x3FF)
        r[4i+2] = Int32(((UInt32(a[5i+2]) >> 2) | (UInt32(a[5i+3]) << 6)) & 0x3FF)
        r[4i+3] = Int32(((UInt32(a[5i+3]) >> 4) | (UInt32(a[5i+4]) << 4)) & 0x3FF)
        r[4i+4] = Int32(((UInt32(a[5i+4]) >> 6) | (UInt32(a[5i+5]) << 2)) & 0x3FF)
    end
    return r
end

function polyeta_pack(a::Vector{Int32})::Vector{UInt8}
    r = zeros(UInt8, POLYETA_PACKED)
    if ETA == 4
        for i in 0:(N÷2 - 1)
            t0 = (ETA - a[2i+1]) % UInt8
            t1 = (ETA - a[2i+2]) % UInt8
            r[i+1] = t0 | (t1 << 4)
        end
    elseif ETA == 2
        for i in 0:(N÷8 - 1)
            t = [(ETA - a[8i+j+1]) % UInt8 for j in 0:7]
            r[3i+1] = t[1] | (t[2] << 3) | (t[3] << 6)
            r[3i+2] = (t[3] >> 2) | (t[4] << 1) | (t[5] << 4) | (t[6] << 7)
            r[3i+3] = (t[6] >> 1) | (t[7] << 2) | (t[8] << 5)
        end
    end
    return r
end

function polyeta_unpack!(r::Vector{Int32}, a::Vector{UInt8})
    if ETA == 4
        for i in 0:(N÷2 - 1)
            r[2i+1] = Int32(ETA) - Int32(a[i+1] & 0x0F)
            r[2i+2] = Int32(ETA) - Int32(a[i+1] >> 4)
        end
    elseif ETA == 2
        for i in 0:(N÷8 - 1)
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
    r = zeros(UInt8, POLYZ_PACKED)
    if GAMMA1 == Int32(1 << 19)
        # 20-bit packing: 2 coefficients -> 5 bytes
        for i in 0:(N÷2 - 1)
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
        for i in 0:(N÷4 - 1)
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
    if GAMMA1 == Int32(1 << 19)
        for i in 0:(N÷2 - 1)
            r[2i+1] = Int32(UInt32(a[5i+1]) | (UInt32(a[5i+2]) << 8) | (UInt32(a[5i+3]) << 16)) & Int32(0xFFFFF)
            r[2i+2] = Int32((UInt32(a[5i+3]) >> 4) | (UInt32(a[5i+4]) << 4) | (UInt32(a[5i+5]) << 12)) & Int32(0xFFFFF)
            r[2i+1] = GAMMA1 - r[2i+1]
            r[2i+2] = GAMMA1 - r[2i+2]
        end
    elseif GAMMA1 == Int32(1 << 17)
        for i in 0:(N÷4 - 1)
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
    r = zeros(UInt8, POLYW1_PACKED)
    if GAMMA2 == Int32(div(Q - 1, 32))
        # 4-bit: 2 coefficients per byte (w1 range 0..15)
        for i in 0:(N÷2 - 1)
            r[i+1] = (a[2i+1] | (a[2i+2] << 4)) % UInt8
        end
    elseif GAMMA2 == Int32(div(Q - 1, 88))
        # 6-bit: 4 coefficients per 3 bytes (w1 range 0..43)
        for i in 0:(N÷4 - 1)
            r[3i+1] = (a[4i+1] | (a[4i+2] << 6)) % UInt8
            r[3i+2] = ((a[4i+2] >> 2) | (a[4i+3] << 4)) % UInt8
            r[3i+3] = ((a[4i+3] >> 4) | (a[4i+4] << 2)) % UInt8
        end
    end
    return r
end
# ==================== KEY GENERATION ====================

function dilithium_keygen_derand(xi::Vector{UInt8})
    seed = xi[1:SEEDBYTES]
    expanded = SHA.shake256(vcat(seed, UInt8[K, L]), UInt64(2*SEEDBYTES + CRHBYTES))
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

    t = [zeros(Int32, N) for _ in 1:K]
    tmp = zeros(Int32, N)
    for i in 1:K
        fill!(t[i], Int32(0))
        for j in 1:L
            poly_pointwise!(tmp, A[i,j], s1hat[j])
            poly_add!(t[i], t[i], tmp)
        end
        poly_reduce!(t[i])
        invntt!(t[i])
        poly_add!(t[i], t[i], s2[i])
        poly_caddq!(t[i])
    end

    # power2round: t = t1*2^D + t0
    t1 = [zeros(Int32, N) for _ in 1:K]
    t0 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        for j in 1:N
            t1[i][j], t0[i][j] = power2round(t[i][j])
        end
    end

    # Pack pk = rho || t1
    pk = copy(rho)
    for i in 1:K; append!(pk, polyt1_pack(t1[i])); end

    # tr = H(pk)
    tr = SHA.shake256(pk, UInt64(TRBYTES))

    # Pack sk = rho || key || tr || s1 || s2 || t0
    sk = vcat(rho, key, tr)
    for i in 1:L; append!(sk, polyeta_pack(s1[i])); end
    for i in 1:K; append!(sk, polyeta_pack(s2[i])); end
    for i in 1:K
        # polyt0_pack: simplified 13-bit packing
        buf = zeros(UInt8, POLYT0_PACKED)
        for idx in 0:(N÷8 - 1)
            ts = [UInt32((1 << (D-1)) - t0[i][8idx+k]) for k in 1:8]
            buf[13idx+1]  = (ts[1]) % UInt8
            buf[13idx+2]  = ((ts[1] >> 8) | (ts[2] << 5)) % UInt8
            buf[13idx+3]  = (ts[2] >> 3) % UInt8
            buf[13idx+4]  = ((ts[2] >> 11) | (ts[3] << 2)) % UInt8
            buf[13idx+5]  = ((ts[3] >> 6) | (ts[4] << 7)) % UInt8
            buf[13idx+6]  = (ts[4] >> 1) % UInt8
            buf[13idx+7]  = ((ts[4] >> 9) | (ts[5] << 4)) % UInt8
            buf[13idx+8]  = (ts[5] >> 4) % UInt8
            buf[13idx+9]  = ((ts[5] >> 12) | (ts[6] << 1)) % UInt8
            buf[13idx+10] = ((ts[6] >> 7) | (ts[7] << 6)) % UInt8
            buf[13idx+11] = (ts[7] >> 2) % UInt8
            buf[13idx+12] = ((ts[7] >> 10) | (ts[8] << 3)) % UInt8
            buf[13idx+13] = (ts[8] >> 5) % UInt8
        end
        append!(sk, buf)
    end

    return pk, sk
end

function dilithium_keygen()
    return dilithium_keygen_derand(rand(UInt8, SEEDBYTES))
end

# ==================== SIGN ====================

function dilithium_sign_derand(msg::Vector{UInt8}, sk::Vector{UInt8}, rnd::Vector{UInt8}; context::Vector{UInt8}=UInt8[])
    # Unpack sk
    pos = 1
    rho = sk[pos:pos+SEEDBYTES-1]; pos += SEEDBYTES
    key = sk[pos:pos+SEEDBYTES-1]; pos += SEEDBYTES
    tr = sk[pos:pos+TRBYTES-1]; pos += TRBYTES

    s1 = [zeros(Int32, N) for _ in 1:L]
    for i in 1:L
        polyeta_unpack!(s1[i], sk[pos:pos+POLYETA_PACKED-1]); pos += POLYETA_PACKED
    end
    s2 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        polyeta_unpack!(s2[i], sk[pos:pos+POLYETA_PACKED-1]); pos += POLYETA_PACKED
    end
    t0 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        # polyt0_unpack (simplified)
        buf = sk[pos:pos+POLYT0_PACKED-1]; pos += POLYT0_PACKED
        for idx in 0:(N÷8 - 1)
            t0[i][8idx+1] = Int32(UInt32(buf[13idx+1]) | (UInt32(buf[13idx+2]) << 8)) & Int32(0x1FFF)
            t0[i][8idx+2] = Int32((UInt32(buf[13idx+2]) >> 5) | (UInt32(buf[13idx+3]) << 3) | (UInt32(buf[13idx+4]) << 11)) & Int32(0x1FFF)
            t0[i][8idx+3] = Int32((UInt32(buf[13idx+4]) >> 2) | (UInt32(buf[13idx+5]) << 6)) & Int32(0x1FFF)
            t0[i][8idx+4] = Int32((UInt32(buf[13idx+5]) >> 7) | (UInt32(buf[13idx+6]) << 1) | (UInt32(buf[13idx+7]) << 9)) & Int32(0x1FFF)
            t0[i][8idx+5] = Int32((UInt32(buf[13idx+7]) >> 4) | (UInt32(buf[13idx+8]) << 4) | (UInt32(buf[13idx+9]) << 12)) & Int32(0x1FFF)
            t0[i][8idx+6] = Int32((UInt32(buf[13idx+9]) >> 1) | (UInt32(buf[13idx+10]) << 7)) & Int32(0x1FFF)
            t0[i][8idx+7] = Int32((UInt32(buf[13idx+10]) >> 6) | (UInt32(buf[13idx+11]) << 2) | (UInt32(buf[13idx+12]) << 10)) & Int32(0x1FFF)
            t0[i][8idx+8] = Int32((UInt32(buf[13idx+12]) >> 3) | (UInt32(buf[13idx+13]) << 5)) & Int32(0x1FFF)
            for k in 1:8
                t0[i][8idx+k] = Int32((1 << (D-1))) - t0[i][8idx+k]
            end
        end
    end

    # Expand A, NTT(s1), NTT(s2), NTT(t0)
    A = [zeros(Int32, N) for _ in 1:K, _ in 1:L]
    for i in 1:K, j in 1:L
        poly_uniform!(A[i,j], rho, UInt16((i-1) << 8 | (j-1)))
    end
    for i in 1:L; ntt!(s1[i]); end
    for i in 1:K; ntt!(s2[i]); end
    for i in 1:K; ntt!(t0[i]); end

    # mu = CRH(tr || msg)
    # FIPS 204 pure mode: pre = [0x00, ctxlen, ctx...]
    pre = vcat(UInt8[0x00, UInt8(length(context))], context)
    mu = SHA.shake256(vcat(tr, pre, msg), UInt64(CRHBYTES))

    # rhoprime = CRH(key || zeros || mu) — deterministic signing
    rhoprime = SHA.shake256(vcat(key, rnd[1:32], mu), UInt64(CRHBYTES))

    nonce = UInt16(0)
    y = [zeros(Int32, N) for _ in 1:L]
    z = [zeros(Int32, N) for _ in 1:L]
    w1 = [zeros(Int32, N) for _ in 1:K]
    w0 = [zeros(Int32, N) for _ in 1:K]
    h = [zeros(Int32, N) for _ in 1:K]
    cp = zeros(Int32, N)
    tmp = zeros(Int32, N)

    while true
        # Sample y
        for i in 1:L
            poly_uniform_gamma1!(y[i], rhoprime, (L * nonce + i - 1) % UInt16)
        end

        # w = Ay (in NTT domain)
        zy = [copy(y[i]) for i in 1:L]
        for i in 1:L; ntt!(zy[i]); end
        for i in 1:K
            fill!(w1[i], Int32(0))
            for j in 1:L
                poly_pointwise!(tmp, A[i,j], zy[j])
                poly_add!(w1[i], w1[i], tmp)
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

        # Challenge
        w1_packed = UInt8[]
        for i in 1:K; append!(w1_packed, polyw1_pack(w1[i])); end
        c_tilde = SHA.shake256(vcat(mu, w1_packed), UInt64(CTILDEBYTES))
        poly_challenge!(cp, c_tilde)
        cp_hat = copy(cp); ntt!(cp_hat)

        # z = cs1 + y
        for i in 1:L
            poly_pointwise!(z[i], cp_hat, s1[i])
            invntt!(z[i])
            poly_add!(z[i], z[i], y[i])
            poly_reduce!(z[i])
        end

        # Check ||z||_inf < GAMMA1 - BETA
        reject = false
        for i in 1:L
            if poly_chknorm(z[i], GAMMA1 - BETA)
                reject = true; break
            end
        end
        if reject; nonce += 1; continue; end

        # w0 = w0 - cs2
        for i in 1:K
            poly_pointwise!(tmp, cp_hat, s2[i])
            invntt!(tmp)
            poly_sub!(w0[i], w0[i], tmp)
            poly_reduce!(w0[i])
        end
        reject = false
        for i in 1:K
            if poly_chknorm(w0[i], GAMMA2 - BETA)
                reject = true; break
            end
        end
        if reject; nonce += 1; continue; end

        # ct0
        for i in 1:K
            poly_pointwise!(h[i], cp_hat, t0[i])
            invntt!(h[i])
            poly_reduce!(h[i])
        end
        reject = false
        for i in 1:K
            if poly_chknorm(h[i], GAMMA2)
                reject = true; break
            end
        end
        if reject; nonce += 1; continue; end

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
        if hints_count > OMEGA; nonce += 1; continue; end

        # Pack signature
        sig = copy(c_tilde)
        for i in 1:L; append!(sig, polyz_pack(z[i])); end
        # Pack h (sparse)
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
end

function dilithium_sign(msg::Vector{UInt8}, sk::Vector{UInt8}; hedged::Bool=false, context::Vector{UInt8}=UInt8[])
    rnd = hedged ? rand(UInt8, 32) : zeros(UInt8, 32)
    return dilithium_sign_derand(msg, sk, rnd; context=context)
end

# ==================== VERIFY ====================

function dilithium_verify(msg::Vector{UInt8}, sig::Vector{UInt8}, pk::Vector{UInt8}; context::Vector{UInt8}=UInt8[])
    length(sig) != SIG_BYTES && return false

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

    # mu = CRH(H(pk) || msg)
    tr = SHA.shake256(pk, UInt64(TRBYTES))
    # FIPS 204 pure mode: pre = [0x00, ctxlen, ctx...]
    pre = vcat(UInt8[0x00, UInt8(length(context))], context)
    mu = SHA.shake256(vcat(tr, pre, msg), UInt64(CRHBYTES))

    # Expand A, compute w1' = Az - c*t1*2^D
    cp = zeros(Int32, N)
    poly_challenge!(cp, c_tilde)

    A = [zeros(Int32, N) for _ in 1:K, _ in 1:L]
    for i in 1:K, j in 1:L
        poly_uniform!(A[i,j], rho, UInt16((i-1) << 8 | (j-1)))
    end

    for i in 1:L; ntt!(z[i]); end
    w1p = [zeros(Int32, N) for _ in 1:K]
    tmp = zeros(Int32, N)
    for i in 1:K
        fill!(w1p[i], Int32(0))
        for j in 1:L
            poly_pointwise!(tmp, A[i,j], z[j])
            poly_add!(w1p[i], w1p[i], tmp)
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

    # Check hint count
    hint_count = sum(sum(h[i][j] for j in 1:N) for i in 1:K)
    hint_count > OMEGA && return false

    # Recompute challenge
    w1_packed = UInt8[]
    for i in 1:K; append!(w1_packed, polyw1_pack(w1p[i])); end
    c2 = SHA.shake256(vcat(mu, w1_packed), UInt64(CTILDEBYTES))

    return c_tilde == c2
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
        return SHA.sha3_256(msg)
    elseif hash_alg == "SHA3-384"
        return SHA.sha3_384(msg)
    elseif hash_alg == "SHA3-512"
        return SHA.sha3_512(msg)
    elseif hash_alg == "SHAKE-128"
        return SHA.shake128(msg, UInt64(32))  # 256 bits
    elseif hash_alg == "SHAKE-256"
        return SHA.shake256(msg, UInt64(64))  # 512 bits
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

"""Sign with pre-hash (HashML-DSA, FIPS 204 Algorithm 4)."""
function dilithium_sign_prehash(msg::Vector{UInt8}, sk::Vector{UInt8}, hash_alg::String;
                                 hedged::Bool=false, context::Vector{UInt8}=UInt8[])
    haskey(HASH_OIDS, hash_alg) || error("Unknown hash: $hash_alg")
    oid = HASH_OIDS[hash_alg]
    ph_m = prehash_message(msg, hash_alg)

    # Construct M\' with domain separator 0x01
    # M\' = 0x01 || ctxlen || ctx || OID || PH(M)
    pre = vcat(UInt8[0x01, UInt8(length(context))], context, oid, ph_m)

    # Call sign_internal with pre as the "message"
    # We need to modify the mu derivation to use pre instead of the standard pure prefix
    # Since sign_derand constructs: pre_pure = [0x00, ctxlen, ctx...]; mu = H(tr || pre_pure || msg)
    # For prehash: mu = H(tr || pre)  (pre already contains everything)
    rnd = hedged ? rand(UInt8, 32) : zeros(UInt8, 32)

    # Extract sk components to compute mu directly
    pos = 1
    rho = sk[pos:pos+SEEDBYTES-1]; pos += SEEDBYTES
    key = sk[pos:pos+SEEDBYTES-1]; pos += SEEDBYTES
    tr = sk[pos:pos+TRBYTES-1]; pos += TRBYTES

    # mu = H(tr || pre) where pre = [0x01, ctxlen, ctx, OID, PH(M)]
    mu = SHA.shake256(vcat(tr, pre), UInt64(CRHBYTES))

    # rhoprime = H(key || rnd || mu)
    rhoprime = SHA.shake256(vcat(key, rnd, mu), UInt64(CRHBYTES))

    # Rest of signing is identical — unpack s1/s2/t0, expand A, rejection loop
    # For simplicity, we construct a fake "message" that when combined with
    # the pure-mode prefix [0x00, 0x00], produces the same mu.
    # Actually easier: just call the existing sign with a custom pre.
    # But sign_derand hardcodes the pure pre. We need a lower-level entry.
    # Let me directly use the sign loop with the computed mu and rhoprime.

    # Unpack s1, s2, t0 from sk
    s1 = [zeros(Int32, N) for _ in 1:L]
    for i in 1:L; polyeta_unpack!(s1[i], sk[pos:pos+POLYETA_PACKED-1]); pos += POLYETA_PACKED; end
    s2 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K; polyeta_unpack!(s2[i], sk[pos:pos+POLYETA_PACKED-1]); pos += POLYETA_PACKED; end
    t0 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        buf = sk[pos:pos+POLYT0_PACKED-1]; pos += POLYT0_PACKED
        for idx in 0:(N÷8 - 1)
            t0[i][8idx+1] = Int32(UInt32(buf[13idx+1]) | (UInt32(buf[13idx+2]) << 8)) & Int32(0x1FFF)
            t0[i][8idx+2] = Int32((UInt32(buf[13idx+2]) >> 5) | (UInt32(buf[13idx+3]) << 3) | (UInt32(buf[13idx+4]) << 11)) & Int32(0x1FFF)
            t0[i][8idx+3] = Int32((UInt32(buf[13idx+4]) >> 2) | (UInt32(buf[13idx+5]) << 6)) & Int32(0x1FFF)
            t0[i][8idx+4] = Int32((UInt32(buf[13idx+5]) >> 7) | (UInt32(buf[13idx+6]) << 1) | (UInt32(buf[13idx+7]) << 9)) & Int32(0x1FFF)
            t0[i][8idx+5] = Int32((UInt32(buf[13idx+7]) >> 4) | (UInt32(buf[13idx+8]) << 4) | (UInt32(buf[13idx+9]) << 12)) & Int32(0x1FFF)
            t0[i][8idx+6] = Int32((UInt32(buf[13idx+9]) >> 1) | (UInt32(buf[13idx+10]) << 7)) & Int32(0x1FFF)
            t0[i][8idx+7] = Int32((UInt32(buf[13idx+10]) >> 6) | (UInt32(buf[13idx+11]) << 2) | (UInt32(buf[13idx+12]) << 10)) & Int32(0x1FFF)
            t0[i][8idx+8] = Int32((UInt32(buf[13idx+12]) >> 3) | (UInt32(buf[13idx+13]) << 5)) & Int32(0x1FFF)
            for k in 1:8; t0[i][8idx+k] = Int32((1 << (D-1))) - t0[i][8idx+k]; end
        end
    end

    # Expand A, NTT(s1), NTT(s2), NTT(t0)
    A = [zeros(Int32, N) for _ in 1:K, _ in 1:L]
    for i in 1:K, j in 1:L; poly_uniform!(A[i,j], rho, UInt16((i-1) << 8 | (j-1))); end
    for i in 1:L; ntt!(s1[i]); end
    for i in 1:K; ntt!(s2[i]); end
    for i in 1:K; ntt!(t0[i]); end

    nonce = UInt16(0)
    y = [zeros(Int32, N) for _ in 1:L]
    z = [zeros(Int32, N) for _ in 1:L]
    w1 = [zeros(Int32, N) for _ in 1:K]
    w0 = [zeros(Int32, N) for _ in 1:K]
    h = [zeros(Int32, N) for _ in 1:K]
    cp = zeros(Int32, N)
    tmp = zeros(Int32, N)

    while true
        for i in 1:L; poly_uniform_gamma1!(y[i], rhoprime, (L * nonce + i - 1) % UInt16); end
        zy = [copy(y[i]) for i in 1:L]
        for i in 1:L; ntt!(zy[i]); end
        for i in 1:K
            fill!(w1[i], Int32(0))
            for j in 1:L; poly_pointwise!(tmp, A[i,j], zy[j]); poly_add!(w1[i], w1[i], tmp); end
            poly_reduce!(w1[i]); invntt!(w1[i]); poly_caddq!(w1[i])
        end
        for i in 1:K; for j in 1:N; w1[i][j], w0[i][j] = decompose(w1[i][j]); end; end
        w1_packed = UInt8[]
        for i in 1:K; append!(w1_packed, polyw1_pack(w1[i])); end
        c_tilde = SHA.shake256(vcat(mu, w1_packed), UInt64(CTILDEBYTES))
        poly_challenge!(cp, c_tilde)
        cp_hat = copy(cp); ntt!(cp_hat)

        for i in 1:L; poly_pointwise!(z[i], cp_hat, s1[i]); invntt!(z[i]); poly_add!(z[i], z[i], y[i]); poly_reduce!(z[i]); end
        reject = any(poly_chknorm(z[i], GAMMA1 - BETA) for i in 1:L)
        if reject; nonce += 1; continue; end

        for i in 1:K; poly_pointwise!(tmp, cp_hat, s2[i]); invntt!(tmp); poly_sub!(w0[i], w0[i], tmp); poly_reduce!(w0[i]); end
        reject = any(poly_chknorm(w0[i], GAMMA2 - BETA) for i in 1:K)
        if reject; nonce += 1; continue; end

        for i in 1:K; poly_pointwise!(h[i], cp_hat, t0[i]); invntt!(h[i]); poly_reduce!(h[i]); end
        reject = any(poly_chknorm(h[i], GAMMA2) for i in 1:K)
        if reject; nonce += 1; continue; end

        for i in 1:K; poly_add!(w0[i], w0[i], h[i]); end
        hints_count = 0
        for i in 1:K; for j in 1:N; h[i][j] = Int32(make_hint(w0[i][j], w1[i][j])); hints_count += h[i][j]; end; end
        if hints_count > OMEGA; nonce += 1; continue; end

        sig = copy(c_tilde)
        for i in 1:L; append!(sig, polyz_pack(z[i])); end
        h_packed = zeros(UInt8, OMEGA + K)
        local k_pos = 0
        for i in 1:K
            for j in 1:N; if h[i][j] != 0; h_packed[k_pos + 1] = UInt8(j - 1); k_pos += 1; end; end
            h_packed[OMEGA + i] = UInt8(k_pos)
        end
        append!(sig, h_packed)
        return sig
    end
end

"""Verify with pre-hash (HashML-DSA, FIPS 204 Algorithm 5)."""
function dilithium_verify_prehash(msg::Vector{UInt8}, sig::Vector{UInt8}, pk::Vector{UInt8},
                                    hash_alg::String; context::Vector{UInt8}=UInt8[])
    haskey(HASH_OIDS, hash_alg) || return false
    oid = HASH_OIDS[hash_alg]
    ph_m = prehash_message(msg, hash_alg)

    # Construct pre with 0x01 domain separator
    pre = vcat(UInt8[0x01, UInt8(length(context))], context, oid, ph_m)

    # Verify uses mu = H(tr || pre) — same as sign_prehash
    # We can call the existing verify but with a custom pre
    # Actually, our verify already constructs pre = [0x00, ctxlen, ctx...] internally
    # We need to override that. Simplest: inline the verify with the prehash pre.

    length(sig) != SIG_BYTES && return false

    rho = pk[1:SEEDBYTES]
    t1 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K; polyt1_unpack!(t1[i], pk[SEEDBYTES+(i-1)*POLYT1_PACKED+1:SEEDBYTES+i*POLYT1_PACKED]); end

    c_tilde = sig[1:CTILDEBYTES]
    z = [zeros(Int32, N) for _ in 1:L]
    local vpos = CTILDEBYTES + 1
    for i in 1:L; polyz_unpack!(z[i], sig[vpos:vpos+POLYZ_PACKED-1]); vpos += POLYZ_PACKED; end

    h = [zeros(Int32, N) for _ in 1:K]
    h_raw = sig[vpos:end]
    local k_pos = 0
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
    for j in (k_pos+1):OMEGA; h_raw[j] != 0 && return false; end

    for i in 1:L; poly_chknorm(z[i], GAMMA1 - BETA) && return false; end

    hint_count = sum(sum(h[i][j] for j in 1:N) for i in 1:K)
    hint_count > OMEGA && return false

    # mu with prehash pre
    tr = SHA.shake256(pk, UInt64(TRBYTES))
    mu = SHA.shake256(vcat(tr, pre), UInt64(CRHBYTES))

    cp = zeros(Int32, N); poly_challenge!(cp, c_tilde)
    A = [zeros(Int32, N) for _ in 1:K, _ in 1:L]
    for i in 1:K, j in 1:L; poly_uniform!(A[i,j], rho, UInt16((i-1) << 8 | (j-1))); end

    for i in 1:L; ntt!(z[i]); end
    w1p = [zeros(Int32, N) for _ in 1:K]
    tmp = zeros(Int32, N)
    for i in 1:K
        fill!(w1p[i], Int32(0))
        for j in 1:L; poly_pointwise!(tmp, A[i,j], z[j]); poly_add!(w1p[i], w1p[i], tmp); end
    end

    ntt!(cp)
    for i in 1:K
        poly_shiftl!(t1[i]); ntt!(t1[i])
        poly_pointwise!(tmp, cp, t1[i])
        poly_sub!(w1p[i], w1p[i], tmp)
        poly_reduce!(w1p[i]); invntt!(w1p[i]); poly_caddq!(w1p[i])
    end

    for i in 1:K; for j in 1:N; w1p[i][j] = use_hint(w1p[i][j], h[i][j] != 0); end; end

    w1_packed = UInt8[]
    for i in 1:K; append!(w1_packed, polyw1_pack(w1p[i])); end
    c2 = SHA.shake256(vcat(mu, w1_packed), UInt64(CTILDEBYTES))

    return c_tilde == c2
end

# ==================== INTERNAL SIGNING INTERFACE ====================

"""Sign with explicit mu (internal interface, externalMu=true)."""
function dilithium_sign_internal_mu(mu::Vector{UInt8}, sk::Vector{UInt8}, rnd::Vector{UInt8})
    length(mu) == CRHBYTES || error("mu must be $CRHBYTES bytes")

    pos = 1
    rho = sk[pos:pos+SEEDBYTES-1]; pos += SEEDBYTES
    key = sk[pos:pos+SEEDBYTES-1]; pos += SEEDBYTES
    tr = sk[pos:pos+TRBYTES-1]; pos += TRBYTES

    s1 = [zeros(Int32, N) for _ in 1:L]
    for i in 1:L; polyeta_unpack!(s1[i], sk[pos:pos+POLYETA_PACKED-1]); pos += POLYETA_PACKED; end
    s2 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K; polyeta_unpack!(s2[i], sk[pos:pos+POLYETA_PACKED-1]); pos += POLYETA_PACKED; end
    t0 = [zeros(Int32, N) for _ in 1:K]
    for i in 1:K
        buf = sk[pos:pos+POLYT0_PACKED-1]; pos += POLYT0_PACKED
        for idx in 0:(N÷8 - 1)
            t0[i][8idx+1] = Int32(UInt32(buf[13idx+1]) | (UInt32(buf[13idx+2]) << 8)) & Int32(0x1FFF)
            t0[i][8idx+2] = Int32((UInt32(buf[13idx+2]) >> 5) | (UInt32(buf[13idx+3]) << 3) | (UInt32(buf[13idx+4]) << 11)) & Int32(0x1FFF)
            t0[i][8idx+3] = Int32((UInt32(buf[13idx+4]) >> 2) | (UInt32(buf[13idx+5]) << 6)) & Int32(0x1FFF)
            t0[i][8idx+4] = Int32((UInt32(buf[13idx+5]) >> 7) | (UInt32(buf[13idx+6]) << 1) | (UInt32(buf[13idx+7]) << 9)) & Int32(0x1FFF)
            t0[i][8idx+5] = Int32((UInt32(buf[13idx+7]) >> 4) | (UInt32(buf[13idx+8]) << 4) | (UInt32(buf[13idx+9]) << 12)) & Int32(0x1FFF)
            t0[i][8idx+6] = Int32((UInt32(buf[13idx+9]) >> 1) | (UInt32(buf[13idx+10]) << 7)) & Int32(0x1FFF)
            t0[i][8idx+7] = Int32((UInt32(buf[13idx+10]) >> 6) | (UInt32(buf[13idx+11]) << 2) | (UInt32(buf[13idx+12]) << 10)) & Int32(0x1FFF)
            t0[i][8idx+8] = Int32((UInt32(buf[13idx+12]) >> 3) | (UInt32(buf[13idx+13]) << 5)) & Int32(0x1FFF)
            for k in 1:8; t0[i][8idx+k] = Int32((1 << (D-1))) - t0[i][8idx+k]; end
        end
    end

    A = [zeros(Int32, N) for _ in 1:K, _ in 1:L]
    for i in 1:K, j in 1:L; poly_uniform!(A[i,j], rho, UInt16((i-1) << 8 | (j-1))); end
    for i in 1:L; ntt!(s1[i]); end
    for i in 1:K; ntt!(s2[i]); end
    for i in 1:K; ntt!(t0[i]); end

    rhoprime = SHA.shake256(vcat(key, rnd, mu), UInt64(CRHBYTES))

    nonce = UInt16(0)
    y = [zeros(Int32, N) for _ in 1:L]
    z = [zeros(Int32, N) for _ in 1:L]
    w1 = [zeros(Int32, N) for _ in 1:K]
    w0 = [zeros(Int32, N) for _ in 1:K]
    h = [zeros(Int32, N) for _ in 1:K]
    cp = zeros(Int32, N)
    tmp = zeros(Int32, N)

    while true
        for i in 1:L; poly_uniform_gamma1!(y[i], rhoprime, (L * nonce + i - 1) % UInt16); end
        zy = [copy(y[i]) for i in 1:L]
        for i in 1:L; ntt!(zy[i]); end
        for i in 1:K
            fill!(w1[i], Int32(0))
            for j in 1:L; poly_pointwise!(tmp, A[i,j], zy[j]); poly_add!(w1[i], w1[i], tmp); end
            poly_reduce!(w1[i]); invntt!(w1[i]); poly_caddq!(w1[i])
        end
        for i in 1:K; for j in 1:N; w1[i][j], w0[i][j] = decompose(w1[i][j]); end; end
        w1_packed = UInt8[]; for i in 1:K; append!(w1_packed, polyw1_pack(w1[i])); end
        c_tilde = SHA.shake256(vcat(mu, w1_packed), UInt64(CTILDEBYTES))
        poly_challenge!(cp, c_tilde); cp_hat = copy(cp); ntt!(cp_hat)

        for i in 1:L; poly_pointwise!(z[i], cp_hat, s1[i]); invntt!(z[i]); poly_add!(z[i], z[i], y[i]); poly_reduce!(z[i]); end
        any(poly_chknorm(z[i], GAMMA1 - BETA) for i in 1:L) && (nonce += 1; continue)

        for i in 1:K; poly_pointwise!(tmp, cp_hat, s2[i]); invntt!(tmp); poly_sub!(w0[i], w0[i], tmp); poly_reduce!(w0[i]); end
        any(poly_chknorm(w0[i], GAMMA2 - BETA) for i in 1:K) && (nonce += 1; continue)

        for i in 1:K; poly_pointwise!(h[i], cp_hat, t0[i]); invntt!(h[i]); poly_reduce!(h[i]); end
        any(poly_chknorm(h[i], GAMMA2) for i in 1:K) && (nonce += 1; continue)

        for i in 1:K; poly_add!(w0[i], w0[i], h[i]); end
        hints_count = 0
        for i in 1:K; for j in 1:N; h[i][j] = Int32(make_hint(w0[i][j], w1[i][j])); hints_count += h[i][j]; end; end
        hints_count > OMEGA && (nonce += 1; continue)

        sig = copy(c_tilde)
        for i in 1:L; append!(sig, polyz_pack(z[i])); end
        h_packed = zeros(UInt8, OMEGA + K); local k_pos = 0
        for i in 1:K; for j in 1:N; if h[i][j] != 0; h_packed[k_pos+1] = UInt8(j-1); k_pos += 1; end; end; h_packed[OMEGA+i] = UInt8(k_pos); end
        append!(sig, h_packed)
        return sig
    end
end

"""Sign with message (internal interface, externalMu=false). Computes mu = H(tr || msg)."""
function dilithium_sign_internal_msg(msg::Vector{UInt8}, sk::Vector{UInt8}, rnd::Vector{UInt8})
    tr = sk[2*SEEDBYTES+1:2*SEEDBYTES+TRBYTES]
    mu = SHA.shake256(vcat(tr, msg), UInt64(CRHBYTES))
    return dilithium_sign_internal_mu(mu, sk, rnd)
end

