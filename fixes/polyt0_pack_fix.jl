# Fix for Critical #3: polyt0_pack inlined twice — extract to shared function
# Source: dilithium/ref/packing.c lines 664-763
# D = 13, N = 256, POLYT0_PACKEDBYTES = 416

"""
    polyt0_pack!(r::Vector{UInt8}, a::Vector{Int32})

Bit-pack polynomial t0 with coefficients in ]-2^(D-1), 2^(D-1)] into bytes.
8 coefficients → 13 bytes, total: 416 bytes.
Centering: t = (1 << (D-1)) - coeff = 4096 - coeff.

C ref: packing.c lines 664-702.
"""
function polyt0_pack!(r::Vector{UInt8}, a::Vector{Int32})
    for i in 0:31  # 32 groups of 8 coefficients
        t0 = UInt32((1 << 12) - a[8*i+1])
        t1 = UInt32((1 << 12) - a[8*i+2])
        t2 = UInt32((1 << 12) - a[8*i+3])
        t3 = UInt32((1 << 12) - a[8*i+4])
        t4 = UInt32((1 << 12) - a[8*i+5])
        t5 = UInt32((1 << 12) - a[8*i+6])
        t6 = UInt32((1 << 12) - a[8*i+7])
        t7 = UInt32((1 << 12) - a[8*i+8])

        b = 13*i + 1  # byte base, 1-indexed

        r[b+ 0]  = UInt8(t0)
        r[b+ 1]  = UInt8(t0 >> 8)
        r[b+ 1] |= UInt8(t1 << 5)
        r[b+ 2]  = UInt8(t1 >> 3)
        r[b+ 3]  = UInt8(t1 >> 11)
        r[b+ 3] |= UInt8(t2 << 2)
        r[b+ 4]  = UInt8(t2 >> 6)
        r[b+ 4] |= UInt8(t3 << 7)
        r[b+ 5]  = UInt8(t3 >> 1)
        r[b+ 6]  = UInt8(t3 >> 9)
        r[b+ 6] |= UInt8(t4 << 4)
        r[b+ 7]  = UInt8(t4 >> 4)
        r[b+ 8]  = UInt8(t4 >> 12)
        r[b+ 8] |= UInt8(t5 << 1)
        r[b+ 9]  = UInt8(t5 >> 7)
        r[b+ 9] |= UInt8(t6 << 6)
        r[b+10]  = UInt8(t6 >> 2)
        r[b+11]  = UInt8(t6 >> 10)
        r[b+11] |= UInt8(t7 << 3)
        r[b+12]  = UInt8(t7 >> 5)
    end
    return r
end

"""
    polyt0_unpack!(r::Vector{Int32}, a::Vector{UInt8})

Unpack bytes to polynomial t0 with coefficients in ]-2^(D-1), 2^(D-1)].
13 bytes → 8 coefficients, total: 416 bytes → 256 coefficients.

C ref: packing.c lines 712-763.
"""
function polyt0_unpack!(r::Vector{Int32}, a::Vector{UInt8})
    for i in 0:31
        b = 13*i + 1
        c = 8*i

        r[c+1] = Int32((UInt32(a[b+ 0])       | (UInt32(a[b+ 1]) << 8)) & 0x1FFF)
        r[c+2] = Int32((UInt32(a[b+ 1]) >> 5  | (UInt32(a[b+ 2]) << 3) | (UInt32(a[b+ 3]) << 11)) & 0x1FFF)
        r[c+3] = Int32((UInt32(a[b+ 3]) >> 2  | (UInt32(a[b+ 4]) << 6)) & 0x1FFF)
        r[c+4] = Int32((UInt32(a[b+ 4]) >> 7  | (UInt32(a[b+ 5]) << 1) | (UInt32(a[b+ 6]) << 9)) & 0x1FFF)
        r[c+5] = Int32((UInt32(a[b+ 6]) >> 4  | (UInt32(a[b+ 7]) << 4) | (UInt32(a[b+ 8]) << 12)) & 0x1FFF)
        r[c+6] = Int32((UInt32(a[b+ 8]) >> 1  | (UInt32(a[b+ 9]) << 7)) & 0x1FFF)
        r[c+7] = Int32((UInt32(a[b+ 9]) >> 6  | (UInt32(a[b+10]) << 2) | (UInt32(a[b+11]) << 10)) & 0x1FFF)
        r[c+8] = Int32((UInt32(a[b+11]) >> 3  | (UInt32(a[b+12]) << 5)) & 0x1FFF)

        # Undo centering: coeff = 4096 - t
        for j in 1:8
            r[c+j] = Int32(1 << 12) - r[c+j]
        end
    end
    return r
end