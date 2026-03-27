# ── rkprf: SHAKE-256(key || input) ───────────────────────────────────────────
# Mirrors symmetric-shake.c:kyber_shake256_rkprf

function kyber_rkprf(key::Vector{UInt8}, input::Vector{UInt8})
    return SHA.shake256(vcat(key, input), UInt64(KYBER_SSBYTES))
end

# ── Polyvec Serialization ────────────────────────────────────────────────────
# Mirrors polyvec.c: polyvec_tobytes / polyvec_frombytes

function kyber_polyvec_tobytes!(r::AbstractVector{UInt8}, a::Vector{Vector{Int16}}, k::Int)
    for i in 1:k
        offset = (i - 1) * KYBER_POLYBYTES
        kyber_poly_tobytes!(@view(r[offset+1:offset+KYBER_POLYBYTES]), a[i])
    end
    return r
end

function kyber_polyvec_frombytes!(r::Vector{Vector{Int16}}, a::AbstractVector{UInt8}, k::Int)
    for i in 1:k
        offset = (i - 1) * KYBER_POLYBYTES
        kyber_poly_frombytes!(r[i], @view(a[offset+1:offset+KYBER_POLYBYTES]))
    end
    return r
end

# ── Polyvec NTT / INTT ──────────────────────────────────────────────────────
# Mirrors polyvec.c: polyvec_ntt / polyvec_invntt_tomont
# Note: poly_ntt in C does ntt() + poly_reduce(). poly_invntt_tomont does invntt() only.

function kyber_polyvec_ntt!(pv::Vector{Vector{Int16}}, k::Int)
    for i in 1:k
        kyber_ntt!(pv[i])
        kyber_poly_reduce!(pv[i])
    end
    return pv
end

function kyber_polyvec_invntt_tomont!(pv::Vector{Vector{Int16}}, k::Int)
    for i in 1:k
        kyber_invntt!(pv[i])
    end
    return pv
end

# ── Polyvec Basemul Accumulate Montgomery ────────────────────────────────────
# Mirrors polyvec.c:polyvec_basemul_acc_montgomery

function kyber_polyvec_basemul_acc!(r::Vector{Int16},
                                    a::Vector{Vector{Int16}},
                                    b::Vector{Vector{Int16}},
                                    k::Int)
    kyber_poly_basemul_montgomery!(r, a[1], b[1])
    tmp = zeros(Int16, KYBER_N)
    for i in 2:k
        kyber_poly_basemul_montgomery!(tmp, a[i], b[i])
        kyber_poly_add!(r, r, tmp)
    end
    kyber_poly_reduce!(r)
    return r
end

# ── Polyvec Add / Reduce ────────────────────────────────────────────────────

function kyber_polyvec_add!(r::Vector{Vector{Int16}},
                            a::Vector{Vector{Int16}},
                            b::Vector{Vector{Int16}},
                            k::Int)
    for i in 1:k
        kyber_poly_add!(r[i], a[i], b[i])
    end
    return r
end

function kyber_polyvec_reduce!(pv::Vector{Vector{Int16}}, k::Int)
    for i in 1:k
        kyber_poly_reduce!(pv[i])
    end
    return pv
end

# ── Polyvec Compress / Decompress (d=10 for Kyber-768) ──────────────────────
# Mirrors polyvec.c with KYBER_POLYVECCOMPRESSEDBYTES = K*320 (d=10)
# 4 coefficients -> 5 bytes (10 bits each)

function kyber_polyvec_compress!(r::AbstractVector{UInt8},
                                 a::Vector{Vector{Int16}},
                                 k::Int)
    idx = 1
    if KYBER_DU == 10
        # d=10: 4 coefficients → 5 bytes (ML-KEM-512, ML-KEM-768)
        for i in 1:k
            for j in 0:(KYBER_N ÷ 4 - 1)
                t = Vector{UInt16}(undef, 4)
                for m in 0:3
                    u = caddq(a[i][4j + m + 1])
                    d0 = UInt64(u % UInt16)
                    d0 <<= 10; d0 += 1665; d0 *= 1290167; d0 >>= 32
                    t[m + 1] = UInt16(d0 & 0x3ff)
                end
                r[idx]     = (t[1]) % UInt8
                r[idx + 1] = ((t[1] >> 8) | (t[2] << 2)) % UInt8
                r[idx + 2] = ((t[2] >> 6) | (t[3] << 4)) % UInt8
                r[idx + 3] = ((t[3] >> 4) | (t[4] << 6)) % UInt8
                r[idx + 4] = (t[4] >> 2) % UInt8
                idx += 5
            end
        end
    elseif KYBER_DU == 11
        # d=11: 8 coefficients → 11 bytes (ML-KEM-1024)
        for i in 1:k
            for j in 0:(KYBER_N ÷ 8 - 1)
                t = Vector{UInt16}(undef, 8)
                for m in 0:7
                    u = caddq(a[i][8j + m + 1])
                    d0 = UInt64(u % UInt16)
                    d0 <<= 11; d0 += 1664; d0 *= 645084; d0 >>= 31
                    t[m + 1] = UInt16(d0 & 0x7ff)
                end
                r[idx]      = (t[1]) % UInt8
                r[idx + 1]  = ((t[1] >> 8) | (t[2] << 3)) % UInt8
                r[idx + 2]  = ((t[2] >> 5) | (t[3] << 6)) % UInt8
                r[idx + 3]  = (t[3] >> 2) % UInt8
                r[idx + 4]  = ((t[3] >> 10) | (t[4] << 1)) % UInt8
                r[idx + 5]  = ((t[4] >> 7) | (t[5] << 4)) % UInt8
                r[idx + 6]  = ((t[5] >> 4) | (t[6] << 7)) % UInt8
                r[idx + 7]  = (t[6] >> 1) % UInt8
                r[idx + 8]  = ((t[6] >> 9) | (t[7] << 2)) % UInt8
                r[idx + 9]  = ((t[7] >> 6) | (t[8] << 5)) % UInt8
                r[idx + 10] = (t[8] >> 3) % UInt8
                idx += 11
            end
        end
    end
    return r
end

function kyber_polyvec_decompress!(r::Vector{Vector{Int16}},
                                   a::AbstractVector{UInt8},
                                   k::Int)
    idx = 1
    if KYBER_DU == 10
        # d=10: 5 bytes → 4 coefficients
        for i in 1:k
            for j in 0:(KYBER_N ÷ 4 - 1)
                t1 = UInt16(a[idx]) | (UInt16(a[idx+1]) << 8); t1 &= 0x3ff
                t2 = (UInt16(a[idx+1]) >> 2) | (UInt16(a[idx+2]) << 6); t2 &= 0x3ff
                t3 = (UInt16(a[idx+2]) >> 4) | (UInt16(a[idx+3]) << 4); t3 &= 0x3ff
                t4 = (UInt16(a[idx+3]) >> 6) | (UInt16(a[idx+4]) << 2); t4 &= 0x3ff
                r[i][4j+1] = kyber_decompress(t1, 10)
                r[i][4j+2] = kyber_decompress(t2, 10)
                r[i][4j+3] = kyber_decompress(t3, 10)
                r[i][4j+4] = kyber_decompress(t4, 10)
                idx += 5
            end
        end
    elseif KYBER_DU == 11
        # d=11: 11 bytes → 8 coefficients
        for i in 1:k
            for j in 0:(KYBER_N ÷ 8 - 1)
                t = Vector{UInt16}(undef, 8)
                t[1] = (UInt16(a[idx]) | (UInt16(a[idx+1]) << 8)) & 0x7ff
                t[2] = ((UInt16(a[idx+1]) >> 3) | (UInt16(a[idx+2]) << 5)) & 0x7ff
                t[3] = ((UInt16(a[idx+2]) >> 6) | (UInt16(a[idx+3]) << 2) | (UInt16(a[idx+4]) << 10)) & 0x7ff
                t[4] = ((UInt16(a[idx+4]) >> 1) | (UInt16(a[idx+5]) << 7)) & 0x7ff
                t[5] = ((UInt16(a[idx+5]) >> 4) | (UInt16(a[idx+6]) << 4)) & 0x7ff
                t[6] = ((UInt16(a[idx+6]) >> 7) | (UInt16(a[idx+7]) << 1) | (UInt16(a[idx+8]) << 9)) & 0x7ff
                t[7] = ((UInt16(a[idx+8]) >> 2) | (UInt16(a[idx+9]) << 6)) & 0x7ff
                t[8] = ((UInt16(a[idx+9]) >> 5) | (UInt16(a[idx+10]) << 3)) & 0x7ff
                for m in 1:8
                    r[i][8j+m] = kyber_decompress(t[m], 11)
                end
                idx += 11
            end
        end
    end
    return r
end

function kyber_gen_matrix(seed::Vector{UInt8}, transposed::Bool, k::Int)
    a = [[zeros(Int16, KYBER_N) for _ in 1:k] for _ in 1:k]
    for i in 0:(k - 1)
        for j in 0:(k - 1)
            if transposed
                kyber_sample_uniform!(a[i+1][j+1], seed, UInt8(i), UInt8(j))
            else
                kyber_sample_uniform!(a[i+1][j+1], seed, UInt8(j), UInt8(i))
            end
        end
    end
    return a
end

# ── Pack / Unpack Helpers ────────────────────────────────────────────────────

function kyber_pack_pk(t::Vector{Vector{Int16}}, rho::Vector{UInt8}, k::Int)
    pk = Vector{UInt8}(undef, KYBER_INDCPA_PUBLICKEYBYTES)
    kyber_polyvec_tobytes!(@view(pk[1:KYBER_POLYVECBYTES]), t, k)
    pk[KYBER_POLYVECBYTES+1:end] .= rho
    return pk
end

function kyber_unpack_pk(pk::Vector{UInt8}, k::Int)
    t = [zeros(Int16, KYBER_N) for _ in 1:k]
    kyber_polyvec_frombytes!(t, @view(pk[1:KYBER_POLYVECBYTES]), k)
    rho = pk[KYBER_POLYVECBYTES+1:KYBER_POLYVECBYTES+KYBER_SYMBYTES]
    return t, rho
end

function kyber_pack_sk(s::Vector{Vector{Int16}}, k::Int)
    sk = Vector{UInt8}(undef, KYBER_INDCPA_SECRETKEYBYTES)
    kyber_polyvec_tobytes!(sk, s, k)
    return sk
end

function kyber_unpack_sk(sk::Vector{UInt8}, k::Int)
    s = [zeros(Int16, KYBER_N) for _ in 1:k]
    kyber_polyvec_frombytes!(s, sk, k)
    return s
end

function kyber_pack_ciphertext(b::Vector{Vector{Int16}}, v::Vector{Int16}, k::Int)
    ct = Vector{UInt8}(undef, KYBER_CIPHERTEXTBYTES)
    kyber_polyvec_compress!(@view(ct[1:KYBER_POLYVECCOMPRESSEDBYTES]), b, k)
    kyber_poly_compress!(@view(ct[KYBER_POLYVECCOMPRESSEDBYTES+1:end]), v, KYBER_DV)
    return ct
end

function kyber_unpack_ciphertext(ct::Vector{UInt8}, k::Int)
    b = [zeros(Int16, KYBER_N) for _ in 1:k]
    kyber_polyvec_decompress!(b, @view(ct[1:KYBER_POLYVECCOMPRESSEDBYTES]), k)
    v = zeros(Int16, KYBER_N)
    kyber_poly_decompress!(v, @view(ct[KYBER_POLYVECCOMPRESSEDBYTES+1:end]), KYBER_DV)
    return b, v
end

# ══════════════════════════════════════════════════════════════════════════════
# IND-CPA KeyGen (indcpa.c:205 -- indcpa_keypair_derand)
# ══════════════════════════════════════════════════════════════════════════════

function kyber_indcpa_keypair_derand(coins::Vector{UInt8})
    k = KYBER_K

    # buf = coins || KYBER_K byte, then hash_g on 33 bytes
    buf_in = vcat(coins[1:KYBER_SYMBYTES], UInt8[UInt8(k)])
    buf = kyber_hash_g(buf_in)  # 64 bytes
    rho = buf[1:KYBER_SYMBYTES]
    sigma = buf[KYBER_SYMBYTES+1:2*KYBER_SYMBYTES]

    # gen_a: A matrix (not transposed)
    a = kyber_gen_matrix(rho, false, k)

    # Sample secret and error vectors
    skpv = [zeros(Int16, KYBER_N) for _ in 1:k]
    e    = [zeros(Int16, KYBER_N) for _ in 1:k]
    nonce = UInt8(0)
    for i in 1:k
        kyber_poly_getnoise_eta1!(skpv[i], sigma, nonce, KYBER_ETA1)
        nonce += UInt8(1)
    end
    for i in 1:k
        kyber_poly_getnoise_eta1!(e[i], sigma, nonce, KYBER_ETA1)
        nonce += UInt8(1)
    end

    # NTT(s), NTT(e) -- poly_ntt does ntt + reduce
    kyber_polyvec_ntt!(skpv, k)
    kyber_polyvec_ntt!(e, k)

    # t = A * s (in NTT domain), with tomont after basemul
    pkpv = [zeros(Int16, KYBER_N) for _ in 1:k]
    for i in 1:k
        kyber_polyvec_basemul_acc!(pkpv[i], a[i], skpv, k)
        kyber_poly_tomont!(pkpv[i])
    end

    # t = t + e
    kyber_polyvec_add!(pkpv, pkpv, e, k)
    kyber_polyvec_reduce!(pkpv, k)

    # Pack
    sk_bytes = kyber_pack_sk(skpv, k)
    pk_bytes = kyber_pack_pk(pkpv, rho, k)

    return pk_bytes, sk_bytes
end

function kyber_indcpa_keypair()
    coins = rand(Random.RandomDevice(), UInt8, KYBER_SYMBYTES)
    return kyber_indcpa_keypair_derand(coins)
end

# ══════════════════════════════════════════════════════════════════════════════
# IND-CPA Encrypt (indcpa.c:260 -- indcpa_enc)
# ══════════════════════════════════════════════════════════════════════════════

function kyber_indcpa_enc(m::Vector{UInt8}, pk::Vector{UInt8}, coins::Vector{UInt8})
    k = KYBER_K

    # Unpack public key
    t_hat, rho = kyber_unpack_pk(pk, k)

    # Message to polynomial
    msg_poly = zeros(Int16, KYBER_N)
    kyber_poly_frommsg!(msg_poly, m)

    # gen_at: A^T (transposed)
    at = kyber_gen_matrix(rho, true, k)

    # Sample r (sp), e1 (ep), e2 (epp)
    sp  = [zeros(Int16, KYBER_N) for _ in 1:k]
    ep  = [zeros(Int16, KYBER_N) for _ in 1:k]
    epp = zeros(Int16, KYBER_N)
    nonce = UInt8(0)
    for i in 1:k
        kyber_poly_getnoise_eta1!(sp[i], coins, nonce, KYBER_ETA1)
        nonce += UInt8(1)
    end
    for i in 1:k
        kyber_poly_getnoise_eta2!(ep[i], coins, nonce)
        nonce += UInt8(1)
    end
    kyber_poly_getnoise_eta2!(epp, coins, nonce)

    # NTT(r)
    kyber_polyvec_ntt!(sp, k)

    # b = A^T * r  (in NTT domain)
    b = [zeros(Int16, KYBER_N) for _ in 1:k]
    for i in 1:k
        kyber_polyvec_basemul_acc!(b[i], at[i], sp, k)
    end

    # v = t^T * r  (inner product in NTT domain)
    v = zeros(Int16, KYBER_N)
    kyber_polyvec_basemul_acc!(v, t_hat, sp, k)

    # INTT
    kyber_polyvec_invntt_tomont!(b, k)
    kyber_invntt!(v)

    # b = b + e1
    kyber_polyvec_add!(b, b, ep, k)

    # v = v + e2 + msg
    kyber_poly_add!(v, v, epp)
    kyber_poly_add!(v, v, msg_poly)

    # Reduce
    kyber_polyvec_reduce!(b, k)
    kyber_poly_reduce!(v)

    # Pack ciphertext
    ct = kyber_pack_ciphertext(b, v, k)
    return ct
end

# ══════════════════════════════════════════════════════════════════════════════
# IND-CPA Decrypt (indcpa.c:314 -- indcpa_dec)
# ══════════════════════════════════════════════════════════════════════════════

function kyber_indcpa_dec(ct::Vector{UInt8}, sk::Vector{UInt8})
    k = KYBER_K

    # Unpack ciphertext
    b, v = kyber_unpack_ciphertext(ct, k)

    # Unpack secret key
    s_hat = kyber_unpack_sk(sk, k)

    # NTT(u)
    kyber_polyvec_ntt!(b, k)

    # mp = s^T * u (in NTT domain)
    mp = zeros(Int16, KYBER_N)
    kyber_polyvec_basemul_acc!(mp, s_hat, b, k)

    # INTT
    kyber_invntt!(mp)

    # mp = v - mp
    kyber_poly_sub!(mp, v, mp)
    kyber_poly_reduce!(mp)

    # Convert back to message
    msg = Vector{UInt8}(undef, KYBER_INDCPA_MSGBYTES)
    kyber_poly_tomsg!(msg, mp)
    return msg
end

# ══════════════════════════════════════════════════════════════════════════════
# CCA2 KEM KeyGen (kem.c:25 -- crypto_kem_keypair_derand)
# ══════════════════════════════════════════════════════════════════════════════
# sk_kem = sk_cpa || pk || H(pk) || z

function kyber_kem_keypair_derand(coins::Vector{UInt8})
    # coins is 2*KYBER_SYMBYTES = 64 bytes
    # First 32 bytes -> indcpa keygen seed, last 32 bytes -> z
    pk, sk_cpa = kyber_indcpa_keypair_derand(coins[1:KYBER_SYMBYTES])

    # Build KEM secret key: sk_cpa || pk || H(pk) || z
    h_pk = kyber_hash_h(pk)
    z = coins[KYBER_SYMBYTES+1:2*KYBER_SYMBYTES]

    sk = vcat(sk_cpa, pk, h_pk, z)
    return pk, sk
end

function kyber_kem_keypair()
    coins = rand(Random.RandomDevice(), UInt8, 2 * KYBER_SYMBYTES)
    return kyber_kem_keypair_derand(coins)
end

# ══════════════════════════════════════════════════════════════════════════════
# CCA2 KEM Encapsulate (kem.c:76 -- crypto_kem_enc_derand)
# ══════════════════════════════════════════════════════════════════════════════

function kyber_kem_enc_derand(pk::Vector{UInt8}, coins::Vector{UInt8})
    # buf = coins (the message m)
    buf = copy(coins[1:KYBER_SYMBYTES])

    # Multitarget countermeasure: buf = m || H(pk)
    h_pk = kyber_hash_h(pk)
    buf_full = vcat(buf, h_pk)

    # kr = G(m || H(pk))  -> 64 bytes: first 32 = shared secret K, last 32 = encryption coins
    kr = kyber_hash_g(buf_full)

    # Encrypt: ct = indcpa_enc(m, pk, kr[33:64])
    ct = kyber_indcpa_enc(buf, pk, kr[KYBER_SYMBYTES+1:2*KYBER_SYMBYTES])

    # Shared secret = first 32 bytes of kr
    ss = kr[1:KYBER_SYMBYTES]

    return ct, ss
end

function kyber_kem_enc(pk::Vector{UInt8})
    coins = rand(Random.RandomDevice(), UInt8, KYBER_SYMBYTES)
    return kyber_kem_enc_derand(pk, coins)
end

# ══════════════════════════════════════════════════════════════════════════════
# CCA2 KEM Decapsulate (kem.c:140 -- crypto_kem_dec)
# ══════════════════════════════════════════════════════════════════════════════

function kyber_kem_dec(ct::Vector{UInt8}, sk::Vector{UInt8})
    # Parse the KEM secret key
    sk_cpa = sk[1:KYBER_INDCPA_SECRETKEYBYTES]
    pk     = sk[KYBER_INDCPA_SECRETKEYBYTES+1:KYBER_INDCPA_SECRETKEYBYTES+KYBER_PUBLICKEYBYTES]
    h_pk   = sk[KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES+1:KYBER_SECRETKEYBYTES-KYBER_SYMBYTES]
    z      = sk[KYBER_SECRETKEYBYTES-KYBER_SYMBYTES+1:KYBER_SECRETKEYBYTES]

    # Decrypt: m' = indcpa_dec(ct, sk_cpa)
    buf = kyber_indcpa_dec(ct, sk_cpa)

    # buf_full = m' || H(pk)  (H(pk) was stored in sk)
    buf_full = vcat(buf, h_pk)

    # kr = G(m' || H(pk))
    kr = kyber_hash_g(buf_full)

    # Re-encrypt: cmp = indcpa_enc(m', pk, kr[33:64])
    cmp = kyber_indcpa_enc(buf, pk, kr[KYBER_SYMBYTES+1:2*KYBER_SYMBYTES])

    # Constant-time verify
    fail = kyber_verify(ct, cmp)

    # Compute rejection key: rkprf(z, ct) = SHAKE256(z || ct)
    ss = kyber_rkprf(z, ct)

    # If verification passed (fail==0), copy true key kr[1:32] into ss
    # cmov copies x into r when b==1; C does cmov(ss, kr, KYBER_SYMBYTES, !fail)
    # !fail: fail=0 -> b=1 (copy), fail=1 -> b=0 (keep rejection key)
    b = UInt8(1 - fail)
    kyber_cmov!(ss, kr[1:KYBER_SYMBYTES], b)

    return ss
end
