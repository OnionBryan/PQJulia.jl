"""
Keccak-f[1600] and the FIPS 202 functions the package uses: SHAKE128, SHAKE256 and SHA3-224/256/
384/512. The 25-lane state is an immutable tuple, so it stays in registers and on the stack and
leaves no heap copy; the round function is straight-line code with constant rotation counts.
Byte-identical to SHA.jl (test/keccak.jl).
"""
module Keccak

export shake128, shake256, sha3_224, sha3_256, sha3_384, sha3_512

const State = NTuple{25,UInt64}

const RC = (0x0000000000000001, 0x0000000000008082, 0x800000000000808a, 0x8000000080008000,
            0x000000000000808b, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
            0x000000000000008a, 0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
            0x000000008000808b, 0x800000000000008b, 0x8000000000008089, 0x8000000000008003,
            0x8000000000008002, 0x8000000000000080, 0x000000000000800a, 0x800000008000000a,
            0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008)

# ρ offsets r[x, y] (FIPS 202 Table 2); lane (x, y) is tuple index x + 5y + 1.
const ROT = ((0, 36, 3, 41, 18), (1, 44, 10, 45, 2), (62, 6, 43, 15, 61),
             (28, 55, 25, 21, 56), (27, 20, 39, 8, 14))

# One round as straight-line code on locals a1..a25: θ, ρ and π into b, χ back into a, ι.
function round_body()
    a(x, y) = Symbol(:a, mod(x, 5) + 5mod(y, 5) + 1)
    b(x, y) = Symbol(:b, mod(x, 5) + 5mod(y, 5) + 1)
    c(x) = Symbol(:c, mod(x, 5))
    d(x) = Symbol(:d, x)
    ex = Expr[]
    for x in 0:4
        push!(ex, :($(c(x)) = $(a(x, 0)) ⊻ $(a(x, 1)) ⊻ $(a(x, 2)) ⊻ $(a(x, 3)) ⊻ $(a(x, 4))))
    end
    for x in 0:4
        push!(ex, :($(d(x)) = $(c(x - 1)) ⊻ bitrotate($(c(x + 1)), 1)))
    end
    for x in 0:4, y in 0:4
        push!(ex, :($(b(y, 2x + 3y)) = bitrotate($(a(x, y)) ⊻ $(d(x)), $(ROT[x+1][y+1]))))
    end
    for x in 0:4, y in 0:4
        push!(ex, :($(a(x, y)) = $(b(x, y)) ⊻ (~$(b(x + 1, y)) & $(b(x + 2, y)))))
    end
    push!(ex, :(a1 ⊻= rc))
    Expr(:block, ex...)
end

const LANES = [Symbol(:a, i) for i in 1:25]

@eval function keccak_f1600(s::State)
    ($(LANES...),) = s
    for rc in RC
        $(round_body())
    end
    ($(LANES...),)
end

# Lane j (0-based) of the final block: the message tail of `nrem` bytes at `off`, then the
# domain byte `ds`, zeros, and 0x80 in the last byte of the `rate`-byte block.
@inline function tail_lane(data::AbstractVector{UInt8}, off::Int, nrem::Int, j::Int, ds::UInt8, rate::Int)
    v = UInt64(0)
    @inbounds for k in 0:7
        p = 8j + k
        byte = p < nrem ? data[off+p+1] : (p == nrem ? ds : 0x00)
        byte ⊻= p == rate - 1 ? 0x80 : 0x00
        v |= UInt64(byte) << (8k)
    end
    v
end

@inline function load_lane(data::AbstractVector{UInt8}, o::Int)
    @inbounds UInt64(data[o+1]) | UInt64(data[o+2]) << 8 | UInt64(data[o+3]) << 16 |
              UInt64(data[o+4]) << 24 | UInt64(data[o+5]) << 32 | UInt64(data[o+6]) << 40 |
              UInt64(data[o+7]) << 48 | UInt64(data[o+8]) << 56
end

# State ⊻ one block: lanes 1..W from `data` at `off` (whole block, or the padded tail of `nrem`
# bytes). Functions rather than closures: a closure over the reassigned state would box it.
@inline xor_block(st::State, data, off::Int, ::Val{W}) where {W} =
    ntuple(i -> i <= W ? st[i] ⊻ load_lane(data, off + 8(i-1)) : st[i], Val(25))
@inline xor_tail(st::State, data, off::Int, nrem::Int, ds::UInt8, ::Val{W}) where {W} =
    ntuple(i -> i <= W ? st[i] ⊻ tail_lane(data, off, nrem, i - 1, ds, 8W) : st[i], Val(25))

# Absorb `data` at rate R bytes with domain byte `ds` and pad10*1; returns the permuted state.
function absorb(::Val{R}, ds::UInt8, data::AbstractVector{UInt8}) where {R}
    st = ntuple(_ -> UInt64(0), Val(25))
    off = firstindex(data) - 1; stop = off + length(data)
    while stop - off >= R
        st = keccak_f1600(xor_block(st, data, off, Val(R ÷ 8)))
        off += R
    end
    keccak_f1600(xor_tail(st, data, off, stop - off, ds, Val(R ÷ 8)))
end

# Bytes 1..n of the state, little-endian, into out[pos+1 : pos+n]: whole lanes, then the rest.
@inline function store_block!(out::Vector{UInt8}, pos::Int, st::State, n::Int)
    full = n >> 3
    GC.@preserve out begin
        p = Ptr{UInt64}(pointer(out, pos + 1))
        for j in 1:full
            unsafe_store!(p, htol(st[j]), j)              # unaligned store
        end
    end
    @inbounds for i in 8full:n-1
        out[pos+i+1] = (st[(i >> 3) + 1] >> (8 * (i & 7))) % UInt8
    end
end

# Squeeze `outlen` bytes at rate R bytes.
function squeeze(::Val{R}, st::State, outlen::Int) where {R}
    out = Vector{UInt8}(undef, outlen)
    pos = 0
    while true
        n = min(R, outlen - pos)
        store_block!(out, pos, st, n)
        pos += n
        pos == outlen && return out
        st = keccak_f1600(st)
    end
end

"SHAKE128(data, outlen) (FIPS 202)."
shake128(data::AbstractVector{UInt8}, outlen::Integer) = squeeze(Val(168), absorb(Val(168), 0x1f, data), Int(outlen))
"SHAKE256(data, outlen) (FIPS 202)."
shake256(data::AbstractVector{UInt8}, outlen::Integer) = squeeze(Val(136), absorb(Val(136), 0x1f, data), Int(outlen))
"SHA3-224(data) (FIPS 202)."
sha3_224(data::AbstractVector{UInt8}) = squeeze(Val(144), absorb(Val(144), 0x06, data), 28)
"SHA3-256(data) (FIPS 202)."
sha3_256(data::AbstractVector{UInt8}) = squeeze(Val(136), absorb(Val(136), 0x06, data), 32)
"SHA3-384(data) (FIPS 202)."
sha3_384(data::AbstractVector{UInt8}) = squeeze(Val(104), absorb(Val(104), 0x06, data), 48)
"SHA3-512(data) (FIPS 202)."
sha3_512(data::AbstractVector{UInt8}) = squeeze(Val(72), absorb(Val(72), 0x06, data), 64)

end # module Keccak
