# SLH-DSA Parameters and Address Structure
# Ported from sphincsplus/ref/params/*.h + shake_offsets.h

# Parameter struct
struct SLHParams
    n::Int          # Hash output bytes
    full_height::Int
    d::Int          # Number of layers
    fors_height::Int
    fors_trees::Int
    wots_w::Int     # Winternitz parameter (always 16)
    # Derived
    wots_logw::Int
    wots_len1::Int
    wots_len2::Int
    wots_len::Int
    wots_bytes::Int
    tree_height::Int
    fors_msg_bytes::Int
    fors_bytes::Int
    sig_bytes::Int
    pk_bytes::Int
    sk_bytes::Int
end

function SLHParams(n, full_height, d, fors_height, fors_trees; wots_w=16)
    wots_logw = wots_w == 16 ? 4 : 8
    wots_len1 = 8 * n ÷ wots_logw
    wots_len2 = n <= 8 ? 2 : (n <= 136 ? 3 : 4)  # for W=16
    wots_len = wots_len1 + wots_len2
    wots_bytes = wots_len * n
    tree_height = full_height ÷ d
    fors_msg_bytes = (fors_height * fors_trees + 7) ÷ 8
    fors_bytes = (fors_height + 1) * fors_trees * n
    sig_bytes = n + fors_bytes + d * wots_bytes + full_height * n
    pk_bytes = 2 * n
    sk_bytes = 2 * n + pk_bytes
    SLHParams(n, full_height, d, fors_height, fors_trees, wots_w,
              wots_logw, wots_len1, wots_len2, wots_len, wots_bytes,
              tree_height, fors_msg_bytes, fors_bytes, sig_bytes, pk_bytes, sk_bytes)
end

# All 6 SHAKE parameter sets
const SHAKE_128F = SLHParams(16, 66, 22, 6, 33)
const SHAKE_128S = SLHParams(16, 63, 7, 12, 14)
const SHAKE_192F = SLHParams(24, 66, 22, 8, 33)
const SHAKE_192S = SLHParams(24, 63, 7, 14, 17)
const SHAKE_256F = SLHParams(32, 68, 17, 9, 35)
const SHAKE_256S = SLHParams(32, 64, 8, 14, 22)

# ADRS — 32-byte address structure (SHAKE offsets)
const ADDR_BYTES = 32
const OFFSET_LAYER      = 4   # byte 3 (1-indexed = 4)
const OFFSET_TREE       = 9   # bytes 8-15 (1-indexed = 9)
const OFFSET_TYPE       = 20  # byte 19 (1-indexed = 20)
const OFFSET_KP_ADDR    = 21  # bytes 20-23
const OFFSET_CHAIN_ADDR = 28  # byte 27
const OFFSET_HASH_ADDR  = 32  # byte 31
const OFFSET_TREE_HGT   = 28  # byte 27 (overlaps CHAIN_ADDR)
const OFFSET_TREE_INDEX  = 29  # bytes 28-31 (overlaps HASH_ADDR area)

# Address types
const ADDR_TYPE_WOTS     = UInt8(0)
const ADDR_TYPE_WOTSPK   = UInt8(1)
const ADDR_TYPE_HASHTREE = UInt8(2)
const ADDR_TYPE_FORSTREE = UInt8(3)
const ADDR_TYPE_FORSPK   = UInt8(4)
const ADDR_TYPE_WOTSPRF  = UInt8(5)
const ADDR_TYPE_FORSPRF  = UInt8(6)

# Address manipulation (matches address.c)
function set_layer_addr!(addr::Vector{UInt8}, layer::Integer)
    addr[OFFSET_LAYER] = UInt8(layer)
end

function set_tree_addr!(addr::Vector{UInt8}, tree::Integer)
    # 8 bytes big-endian at OFFSET_TREE
    for i in 0:7
        addr[OFFSET_TREE + 7 - i] = UInt8((tree >> (8*i)) & 0xff)
    end
end

function set_type!(addr::Vector{UInt8}, t::UInt8)
    addr[OFFSET_TYPE] = t
    # Setting type clears the 12 trailing bytes (C ref: memset to 0)
    for i in (OFFSET_TYPE+1):ADDR_BYTES
        addr[i] = 0x00
    end
end

function set_keypair_addr!(addr::Vector{UInt8}, kp::Integer)
    # 4 bytes big-endian at OFFSET_KP_ADDR
    addr[OFFSET_KP_ADDR]     = UInt8((kp >> 24) & 0xff)
    addr[OFFSET_KP_ADDR + 1] = UInt8((kp >> 16) & 0xff)
    addr[OFFSET_KP_ADDR + 2] = UInt8((kp >> 8) & 0xff)
    addr[OFFSET_KP_ADDR + 3] = UInt8(kp & 0xff)
end

function copy_keypair_addr!(out::Vector{UInt8}, inp::Vector{UInt8})
    out[1:OFFSET_KP_ADDR+3] .= inp[1:OFFSET_KP_ADDR+3]
    for i in (OFFSET_KP_ADDR+4):ADDR_BYTES; out[i] = 0x00; end
end

function set_chain_addr!(addr::Vector{UInt8}, chain::Integer)
    addr[OFFSET_CHAIN_ADDR] = UInt8(chain)
end

function set_hash_addr!(addr::Vector{UInt8}, hash::Integer)
    addr[OFFSET_HASH_ADDR] = UInt8(hash)
end

function set_tree_height!(addr::Vector{UInt8}, h::Integer)
    addr[OFFSET_TREE_HGT] = UInt8(h)
end

function set_tree_index!(addr::Vector{UInt8}, idx::Integer)
    # 4 bytes big-endian at OFFSET_TREE_INDEX
    addr[OFFSET_TREE_INDEX]     = UInt8((idx >> 24) & 0xff)
    addr[OFFSET_TREE_INDEX + 1] = UInt8((idx >> 16) & 0xff)
    addr[OFFSET_TREE_INDEX + 2] = UInt8((idx >> 8) & 0xff)
    addr[OFFSET_TREE_INDEX + 3] = UInt8(idx & 0xff)
end

function copy_subtree_addr!(out::Vector{UInt8}, inp::Vector{UInt8})
    # Copy layer + tree address (first OFFSET_TYPE bytes)
    out[1:OFFSET_TYPE-1] .= inp[1:OFFSET_TYPE-1]
end
