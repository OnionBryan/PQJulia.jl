# SLH-DSA Hash Functions (SHAKE variant)
# Ported from sphincsplus/ref/hash_shake.c + thash_shake_simple.c

using SHA

# Context: holds pub_seed and sk_seed
struct SLHContext
    pub_seed::Vector{UInt8}
    sk_seed::Vector{UInt8}
end

# thash — tweakable hash: SHAKE256(pub_seed || addr || in) → n bytes
# Mirrors thash_shake_simple.c
function thash(p::SLHParams, ctx::SLHContext, addr::Vector{UInt8},
               input::Vector{UInt8})
    buf = vcat(ctx.pub_seed, addr, input)
    return SHA.shake256(buf, UInt64(p.n))
end

# PRF — pseudo-random function for secret key derivation
# SHAKE256(pub_seed || addr || sk_seed) → n bytes
function prf_addr(p::SLHParams, ctx::SLHContext, addr::Vector{UInt8})
    buf = vcat(ctx.pub_seed, addr, ctx.sk_seed)
    return SHA.shake256(buf, UInt64(p.n))
end

# gen_message_random — generate randomizer R for signing
# SHAKE256(sk_prf || optrand || msg) → n bytes
function gen_message_random(p::SLHParams, sk_prf::Vector{UInt8},
                             optrand::Vector{UInt8}, msg::Vector{UInt8})
    return SHA.shake256(vcat(sk_prf, optrand, msg), UInt64(p.n))
end

# hash_message — hash message to get FORS digest + tree index + leaf index
# SHAKE256(R || PK || msg) → fors_msg_bytes + ceil(full_height - tree_height * 8) / 8 + ...
function hash_message(p::SLHParams, R::Vector{UInt8}, pk::Vector{UInt8}, msg::Vector{UInt8})
    tree_bits = p.full_height - p.tree_height
    leaf_bits = p.tree_height
    tree_bytes = (tree_bits + 7) ÷ 8
    leaf_bytes = (leaf_bits + 7) ÷ 8
    total = p.fors_msg_bytes + tree_bytes + leaf_bytes

    buf = SHA.shake256(vcat(R, pk, msg), UInt64(total))

    md = buf[1:p.fors_msg_bytes]

    # Extract tree index
    tree_val = UInt64(0)
    for i in 1:tree_bytes
        tree_val = (tree_val << 8) | UInt64(buf[p.fors_msg_bytes + i])
    end
    tree_val &= (UInt64(1) << tree_bits) - 1

    # Extract leaf index
    leaf_val = UInt32(0)
    for i in 1:leaf_bytes
        leaf_val = (leaf_val << 8) | UInt32(buf[p.fors_msg_bytes + tree_bytes + i])
    end
    leaf_val &= (UInt32(1) << leaf_bits) - 1

    return md, tree_val, leaf_val
end

# Utility: big-endian uint64 to bytes
function u64_to_bytes(val::UInt64, len::Int)
    out = zeros(UInt8, len)
    for i in 0:len-1
        out[len - i] = UInt8((val >> (8*i)) & 0xff)
    end
    return out
end

function u32_to_bytes(val::UInt32)
    return UInt8[(val >> 24) & 0xff, (val >> 16) & 0xff, (val >> 8) & 0xff, val & 0xff]
end

function bytes_to_u64(data::AbstractVector{UInt8})
    val = UInt64(0)
    for b in data
        val = (val << 8) | UInt64(b)
    end
    return val
end
