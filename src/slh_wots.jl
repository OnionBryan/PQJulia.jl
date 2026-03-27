# WOTS+ One-Time Signature
# Ported from sphincsplus/ref/wots.c

# gen_chain: iterate thash from position start for steps iterations
function gen_chain!(out::Vector{UInt8}, inp::Vector{UInt8}, start::Int, steps::Int,
                    p::SLHParams, ctx::SLHContext, addr::Vector{UInt8})
    out[1:p.n] .= inp[1:p.n]
    for i in start:(start + steps - 1)
        i >= p.wots_w && break
        set_hash_addr!(addr, i)
        out[1:p.n] .= thash(p, ctx, addr, out[1:p.n])
    end
end

# base_w: interpret bytes as base-W digits
function base_w(out_len::Int, input::Vector{UInt8}, logw::Int)
    output = zeros(UInt32, out_len)
    w_mask = UInt32((1 << logw) - 1)
    inp_idx = 1; bits = 0; total = UInt8(0)
    for i in 1:out_len
        if bits == 0
            total = input[inp_idx]; inp_idx += 1; bits = 8
        end
        bits -= logw
        output[i] = UInt32((total >> bits) & w_mask)
    end
    return output
end

# wots_checksum: compute checksum over message digits
function wots_checksum(msg_base_w::Vector{UInt32}, p::SLHParams)
    csum = UInt32(0)
    for i in 1:p.wots_len1
        csum += UInt32(p.wots_w - 1) - msg_base_w[i]
    end
    # Left-shift for padding alignment
    csum <<= ((8 - ((p.wots_len2 * p.wots_logw) % 8)) % 8)
    # Convert to bytes then base_w
    csum_bytes_len = (p.wots_len2 * p.wots_logw + 7) ÷ 8
    csum_bytes = zeros(UInt8, csum_bytes_len)
    for i in 0:csum_bytes_len-1
        csum_bytes[csum_bytes_len - i] = UInt8((csum >> (8*i)) & 0xff)
    end
    return base_w(p.wots_len2, csum_bytes, p.wots_logw)
end

# chain_lengths: message → chain positions (msg digits + checksum digits)
function chain_lengths(msg::Vector{UInt8}, p::SLHParams)
    msg_bw = base_w(p.wots_len1, msg, p.wots_logw)
    csum_bw = wots_checksum(msg_bw, p)
    return vcat(msg_bw, csum_bw)
end

# wots_pk_from_sig: recover WOTS+ public key from signature + message
function wots_pk_from_sig(sig::Vector{UInt8}, msg::Vector{UInt8},
                           p::SLHParams, ctx::SLHContext, addr::Vector{UInt8})
    lengths = chain_lengths(msg, p)
    pk = zeros(UInt8, p.wots_bytes)
    tmp = zeros(UInt8, p.n)
    for i in 0:(p.wots_len - 1)
        set_chain_addr!(addr, i)
        gen_chain!(tmp, sig[i*p.n+1:(i+1)*p.n], Int(lengths[i+1]),
                   p.wots_w - 1 - Int(lengths[i+1]), p, ctx, addr)
        pk[i*p.n+1:(i+1)*p.n] .= tmp[1:p.n]
    end
    return pk
end

# wots_gen_sk: generate WOTS+ secret key for chain i
function wots_gen_sk(chain_idx::Int, p::SLHParams, ctx::SLHContext, addr::Vector{UInt8})
    set_chain_addr!(addr, chain_idx)
    set_hash_addr!(addr, 0)
    sk_addr = copy(addr)
    set_type!(sk_addr, ADDR_TYPE_WOTSPRF)
    set_keypair_addr!(sk_addr, 0)  # will be overwritten
    # Copy keypair address from original addr
    copy_keypair_addr!(sk_addr, addr)
    return prf_addr(p, ctx, sk_addr)
end

# wots_sign: generate WOTS+ signature for n-byte message
function wots_sign(msg::Vector{UInt8}, p::SLHParams, ctx::SLHContext, addr::Vector{UInt8})
    lengths = chain_lengths(msg, p)
    sig = zeros(UInt8, p.wots_bytes)
    tmp = zeros(UInt8, p.n)
    sk_addr = copy(addr)
    set_type!(sk_addr, ADDR_TYPE_WOTSPRF)
    copy_keypair_addr!(sk_addr, addr)

    for i in 0:(p.wots_len - 1)
        set_chain_addr!(sk_addr, i)
        set_hash_addr!(sk_addr, 0)
        sk = prf_addr(p, ctx, sk_addr)

        set_chain_addr!(addr, i)
        gen_chain!(tmp, sk, 0, Int(lengths[i+1]), p, ctx, addr)
        sig[i*p.n+1:(i+1)*p.n] .= tmp[1:p.n]
    end
    return sig
end
