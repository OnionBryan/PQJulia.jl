# SLH-DSA Complete Implementation
# Ported from sphincsplus/ref/ (SHAKE simple variant)

module SLHDSA

using SHA
using Random

include("slh_params.jl")
include("slh_hash.jl")
include("slh_wots.jl")

# ==================== TREE UTILITIES ====================
# Ported from utils.c: compute_root, treehash

function compute_root!(root::Vector{UInt8}, leaf::Vector{UInt8}, leaf_idx::UInt32,
                       idx_offset::UInt32, auth_path::Vector{UInt8}, tree_height::Int,
                       p::SLHParams, ctx::SLHContext, addr::Vector{UInt8})
    buf = zeros(UInt8, 2 * p.n)
    current = copy(leaf[1:p.n])

    for i in 0:(tree_height - 1)
        set_tree_height!(addr, i)
        if ((leaf_idx >> i) & 1) == 0
            set_tree_index!(addr, (leaf_idx >> (i+1)) + (idx_offset >> (i+1)))
            buf[1:p.n] .= current
            buf[p.n+1:2*p.n] .= auth_path[i*p.n+1:(i+1)*p.n]
        else
            set_tree_index!(addr, ((leaf_idx - 1) >> (i+1)) + (idx_offset >> (i+1)))
            buf[1:p.n] .= auth_path[i*p.n+1:(i+1)*p.n]
            buf[p.n+1:2*p.n] .= current
        end
        current .= thash(p, ctx, addr, buf)
    end
    root[1:p.n] .= current
end

# Simplified treehash: builds full tree, extracts root + auth path for leaf_idx
function treehash!(root::Vector{UInt8}, auth_path::Vector{UInt8},
                   p::SLHParams, ctx::SLHContext, leaf_idx::UInt32, idx_offset::UInt32,
                   tree_height::Int, gen_leaf::Function, tree_addr::Vector{UInt8})
    num_leaves = 1 << tree_height
    # Stack-based tree building
    stack = Vector{Vector{UInt8}}()
    stack_heights = Int[]

    for i in 0:(num_leaves - 1)
        leaf = gen_leaf(UInt32(i + idx_offset), p, ctx, tree_addr)
        node = copy(leaf)
        node_height = 0

        # Save auth path nodes
        if UInt32(i) == leaf_idx
            # nothing special, but we track it
        end

        # Merge with stack
        while !isempty(stack_heights) && stack_heights[end] == node_height
            sibling = pop!(stack)
            pop!(stack_heights)

            tree_idx = UInt32(i >> (node_height + 1))
            set_tree_height!(tree_addr, node_height)
            set_tree_index!(tree_addr, tree_idx + (idx_offset >> (node_height + 1)))

            # Determine order
            if ((i >> node_height) & 1) == 0
                error("Stack merge on even index should not happen")
            end
            buf = vcat(sibling, node)
            node = thash(p, ctx, tree_addr, buf)
            node_height += 1
        end

        # Extract auth path: if i is the sibling of leaf_idx at some height
        for h in 0:(tree_height - 1)
            if i == Int(leaf_idx ⊻ (1 << h)) && ((i >> h) & 1) != ((leaf_idx >> h) & 1)
                # This doesn't work for stack-based approach
            end
        end

        push!(stack, node)
        push!(stack_heights, node_height)
    end

    root[1:p.n] .= stack[1]
end

# Simpler approach: compute all leaves, build tree in array, extract auth path
function treehash_full(p::SLHParams, ctx::SLHContext, leaf_idx::UInt32, idx_offset::UInt32,
                        tree_height::Int, gen_leaf::Function, tree_addr::Vector{UInt8})
    num_leaves = 1 << tree_height
    # Generate all leaves
    leaves = [gen_leaf(UInt32(i + idx_offset), p, ctx, tree_addr) for i in 0:(num_leaves-1)]

    # Build tree bottom-up
    tree = Vector{Vector{UInt8}}[leaves]
    for h in 0:(tree_height - 1)
        prev = tree[end]
        next = Vector{UInt8}[]
        for i in 0:2:(length(prev) - 1)
            set_tree_height!(tree_addr, h)
            set_tree_index!(tree_addr, UInt32(i ÷ 2) + (idx_offset >> (h + 1)))
            buf = vcat(prev[i + 1], prev[i + 2])
            push!(next, thash(p, ctx, tree_addr, buf))
        end
        push!(tree, next)
    end

    root = tree[end][1]

    # Extract auth path
    auth_path = UInt8[]
    idx = Int(leaf_idx)
    for h in 0:(tree_height - 1)
        sibling_idx = idx ⊻ 1
        append!(auth_path, tree[h + 1][sibling_idx + 1])
        idx >>= 1
    end

    return root, auth_path
end

# ==================== FORS ====================

function message_to_indices(msg::Vector{UInt8}, p::SLHParams)
    indices = zeros(UInt32, p.fors_trees)
    offset = 0
    for i in 1:p.fors_trees
        idx = UInt32(0)
        for j in 0:(p.fors_height - 1)
            byte_idx = offset ÷ 8 + 1
            bit_idx = offset % 8
            idx |= UInt32((msg[byte_idx] >> bit_idx) & 1) << j
            offset += 1
        end
        indices[i] = idx
    end
    return indices
end

function fors_gen_leaf(idx::UInt32, p::SLHParams, ctx::SLHContext, addr::Vector{UInt8})
    sk_addr = copy(addr)
    set_type!(sk_addr, ADDR_TYPE_FORSPRF)
    set_tree_index!(sk_addr, idx)
    sk = prf_addr(p, ctx, sk_addr)

    leaf_addr = copy(addr)
    set_type!(leaf_addr, ADDR_TYPE_FORSTREE)
    set_tree_height!(leaf_addr, 0)
    set_tree_index!(leaf_addr, idx)
    return thash(p, ctx, leaf_addr, sk)
end

function fors_sign(msg_digest::Vector{UInt8}, p::SLHParams, ctx::SLHContext, fors_addr::Vector{UInt8})
    indices = message_to_indices(msg_digest, p)
    sig = UInt8[]
    roots = UInt8[]

    for i in 0:(p.fors_trees - 1)
        idx = indices[i + 1]
        idx_offset = UInt32(i * (1 << p.fors_height))

        # Secret key for selected leaf
        sk_addr = copy(fors_addr)
        set_type!(sk_addr, ADDR_TYPE_FORSPRF)
        set_tree_index!(sk_addr, idx + idx_offset)
        sk = prf_addr(p, ctx, sk_addr)
        append!(sig, sk)

        # Auth path via treehash
        tree_addr = copy(fors_addr)
        set_type!(tree_addr, ADDR_TYPE_FORSTREE)
        root, auth = treehash_full(p, ctx, idx, idx_offset, p.fors_height, fors_gen_leaf, tree_addr)
        append!(sig, auth)
        append!(roots, root)
    end

    # Compress roots
    pk_addr = copy(fors_addr)
    set_type!(pk_addr, ADDR_TYPE_FORSPK)
    fors_pk = thash(p, ctx, pk_addr, roots)

    return sig, fors_pk
end

function fors_pk_from_sig(sig::Vector{UInt8}, msg_digest::Vector{UInt8},
                           p::SLHParams, ctx::SLHContext, fors_addr::Vector{UInt8})
    indices = message_to_indices(msg_digest, p)
    roots = UInt8[]
    sig_offset = 1

    for i in 0:(p.fors_trees - 1)
        idx = indices[i + 1]
        idx_offset = UInt32(i * (1 << p.fors_height))

        # Recover leaf from SK
        sk = sig[sig_offset:sig_offset + p.n - 1]; sig_offset += p.n
        leaf_addr = copy(fors_addr)
        set_type!(leaf_addr, ADDR_TYPE_FORSTREE)
        set_tree_height!(leaf_addr, 0)
        set_tree_index!(leaf_addr, idx + idx_offset)
        leaf = thash(p, ctx, leaf_addr, sk)

        # Compute root from auth path
        auth = sig[sig_offset:sig_offset + p.fors_height * p.n - 1]; sig_offset += p.fors_height * p.n
        tree_addr = copy(fors_addr)
        set_type!(tree_addr, ADDR_TYPE_FORSTREE)
        root = zeros(UInt8, p.n)
        compute_root!(root, leaf, idx, idx_offset, auth, p.fors_height, p, ctx, tree_addr)
        append!(roots, root)
    end

    pk_addr = copy(fors_addr)
    set_type!(pk_addr, ADDR_TYPE_FORSPK)
    return thash(p, ctx, pk_addr, roots)
end

# ==================== MERKLE / XMSS ====================

function wots_gen_leaf(kp_idx::UInt32, p::SLHParams, ctx::SLHContext, tree_addr::Vector{UInt8})
    wots_addr = copy(tree_addr)
    set_type!(wots_addr, ADDR_TYPE_WOTS)
    set_keypair_addr!(wots_addr, kp_idx)

    sk_addr = copy(wots_addr)
    set_type!(sk_addr, ADDR_TYPE_WOTSPRF)
    copy_keypair_addr!(sk_addr, wots_addr)

    pk_buf = UInt8[]
    for i in 0:(p.wots_len - 1)
        set_chain_addr!(sk_addr, i)
        sk = prf_addr(p, ctx, sk_addr)
        set_chain_addr!(wots_addr, i)
        chain_val = zeros(UInt8, p.n)
        gen_chain!(chain_val, sk, 0, p.wots_w - 1, p, ctx, wots_addr)
        append!(pk_buf, chain_val)
    end

    pk_addr = copy(wots_addr)
    set_type!(pk_addr, ADDR_TYPE_WOTSPK)
    copy_keypair_addr!(pk_addr, wots_addr)
    return thash(p, ctx, pk_addr, pk_buf)
end

function merkle_sign(root_val::Vector{UInt8}, p::SLHParams, ctx::SLHContext,
                      wots_addr::Vector{UInt8}, tree_addr::Vector{UInt8}, idx_leaf::UInt32)
    # WOTS+ sign the current root
    set_type!(wots_addr, ADDR_TYPE_WOTS)
    set_keypair_addr!(wots_addr, idx_leaf)
    wots_sig = wots_sign(root_val, p, ctx, wots_addr)

    # Build subtree, get auth path
    ht_addr = copy(tree_addr)
    set_type!(ht_addr, ADDR_TYPE_HASHTREE)
    tree_root, auth_path = treehash_full(p, ctx, idx_leaf, UInt32(0), p.tree_height, wots_gen_leaf, ht_addr)

    return wots_sig, auth_path, tree_root
end

# ==================== TOP-LEVEL API ====================

function slh_keygen(p::SLHParams; seed::Vector{UInt8}=rand(UInt8, 3*p.n))
    sk_seed = seed[1:p.n]
    sk_prf = seed[p.n+1:2*p.n]
    pub_seed = seed[2*p.n+1:3*p.n]

    ctx = SLHContext(pub_seed, sk_seed)

    # Build top tree to get root
    top_tree_addr = zeros(UInt8, ADDR_BYTES)
    set_layer_addr!(top_tree_addr, p.d - 1)
    set_tree_addr!(top_tree_addr, 0)

    ht_addr = copy(top_tree_addr)
    set_type!(ht_addr, ADDR_TYPE_HASHTREE)
    root, _ = treehash_full(p, ctx, UInt32(0), UInt32(0), p.tree_height, wots_gen_leaf, ht_addr)

    pk = vcat(pub_seed, root)
    sk = vcat(sk_seed, sk_prf, pk)

    return pk, sk
end

function slh_sign(msg::Vector{UInt8}, sk::Vector{UInt8}, p::SLHParams;
                   randomize::Bool=true)
    sk_seed = sk[1:p.n]
    sk_prf = sk[p.n+1:2*p.n]
    pub_seed = sk[2*p.n+1:3*p.n]
    pk_root = sk[3*p.n+1:4*p.n]

    ctx = SLHContext(pub_seed, sk_seed)
    pk = vcat(pub_seed, pk_root)

    # Generate randomizer R
    optrand = randomize ? rand(UInt8, p.n) : zeros(UInt8, p.n)
    R = gen_message_random(p, sk_prf, optrand, msg)

    # Hash message to get FORS indices + tree + leaf
    md, tree_idx, leaf_idx = hash_message(p, R, pk, msg)

    sig = copy(R)  # sig starts with R

    # FORS signature
    fors_addr = zeros(UInt8, ADDR_BYTES)
    set_tree_addr!(fors_addr, tree_idx)
    set_type!(fors_addr, ADDR_TYPE_FORSTREE)
    set_keypair_addr!(fors_addr, leaf_idx)

    fors_sig, fors_pk = fors_sign(md, p, ctx, fors_addr)
    append!(sig, fors_sig)

    # Hypertree signature
    ht_root = copy(fors_pk)
    current_tree = tree_idx
    current_leaf = leaf_idx

    for layer in 0:(p.d - 1)
        wots_addr = zeros(UInt8, ADDR_BYTES)
        set_layer_addr!(wots_addr, layer)
        set_tree_addr!(wots_addr, current_tree)

        tree_addr = copy(wots_addr)

        wots_sig, auth_path, tree_root = merkle_sign(ht_root, p, ctx, wots_addr, tree_addr, UInt32(current_leaf))
        append!(sig, wots_sig)
        append!(sig, auth_path)

        ht_root .= tree_root
        current_leaf = UInt32(current_tree & ((1 << p.tree_height) - 1))
        current_tree >>= p.tree_height
    end

    return sig
end

function slh_verify(msg::Vector{UInt8}, sig::Vector{UInt8}, pk::Vector{UInt8}, p::SLHParams)
    length(sig) != p.sig_bytes && return false

    pub_seed = pk[1:p.n]
    pk_root = pk[p.n+1:2*p.n]
    ctx = SLHContext(pub_seed, UInt8[])  # no sk_seed needed for verify

    sig_offset = 1
    R = sig[sig_offset:sig_offset + p.n - 1]; sig_offset += p.n

    md, tree_idx, leaf_idx = hash_message(p, R, pk, msg)

    # FORS verify
    fors_addr = zeros(UInt8, ADDR_BYTES)
    set_tree_addr!(fors_addr, tree_idx)
    set_type!(fors_addr, ADDR_TYPE_FORSTREE)
    set_keypair_addr!(fors_addr, leaf_idx)

    fors_sig_len = p.fors_trees * (p.n + p.fors_height * p.n)
    fors_sig = sig[sig_offset:sig_offset + fors_sig_len - 1]; sig_offset += fors_sig_len
    fors_pk = fors_pk_from_sig(fors_sig, md, p, ctx, fors_addr)

    # Hypertree verify
    root = copy(fors_pk)
    current_tree = tree_idx
    current_leaf = leaf_idx

    for layer in 0:(p.d - 1)
        wots_addr = zeros(UInt8, ADDR_BYTES)
        set_layer_addr!(wots_addr, layer)
        set_tree_addr!(wots_addr, current_tree)
        set_type!(wots_addr, ADDR_TYPE_WOTS)
        set_keypair_addr!(wots_addr, current_leaf)

        wots_sig = sig[sig_offset:sig_offset + p.wots_bytes - 1]; sig_offset += p.wots_bytes
        auth_path = sig[sig_offset:sig_offset + p.tree_height * p.n - 1]; sig_offset += p.tree_height * p.n

        wots_pk = wots_pk_from_sig(wots_sig, root, p, ctx, wots_addr)

        pk_addr = copy(wots_addr)
        set_type!(pk_addr, ADDR_TYPE_WOTSPK)
        copy_keypair_addr!(pk_addr, wots_addr)
        leaf = thash(p, ctx, pk_addr, wots_pk)

        tree_addr = copy(wots_addr)
        set_type!(tree_addr, ADDR_TYPE_HASHTREE)
        compute_root!(root, leaf, UInt32(current_leaf), UInt32(0), auth_path, p.tree_height, p, ctx, tree_addr)

        current_leaf = UInt32(current_tree & ((1 << p.tree_height) - 1))
        current_tree >>= p.tree_height
    end

    return root == pk_root
end

export SHAKE_128F, SHAKE_128S, SHAKE_192F, SHAKE_192S, SHAKE_256F, SHAKE_256S
export slh_keygen, slh_sign, slh_verify

end # module SLHDSA
