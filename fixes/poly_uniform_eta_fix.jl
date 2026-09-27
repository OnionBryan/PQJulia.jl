# Fix for Critical #2: poly_uniform_eta! missing re-squeeze loop
# Source: dilithium/ref/poly.c lines 370-453
# The existing implementation uses a fixed 512-byte buffer and error() on exhaustion.
# The C reference uses a while(ctr < N) re-squeeze loop.

using SHA

const SHAKE256_RATE = 136  # fips202.h; STREAM256_BLOCKBYTES = SHAKE256_RATE

# poly.c lines 384-417
# Sample coefficients in [-eta, eta] via rejection sampling on buf.
# Returns count of accepted coefficients.
function rej_eta(a::Vector{Int32}, offset::Int, len::Int,
                 buf::Vector{UInt8}, buflen::Int, eta::Int)::Int
    ctr = 0
    pos = 1  # 1-indexed

    # C line 394: while(ctr < len && pos < buflen)
    while ctr < len && pos <= buflen
        # C lines 395-396: nibble extraction
        t0 = UInt32(buf[pos] & 0x0F)
        t1 = UInt32(buf[pos] >> 4)
        pos += 1

        if eta == 2
            # C lines 399-406
            if t0 < 15
                t0 = t0 - (205 * t0 >> 10) * 5
                a[offset + ctr + 1] = Int32(2) - Int32(t0)
                ctr += 1
            end
            if t1 < 15 && ctr < len
                t1 = t1 - (205 * t1 >> 10) * 5
                a[offset + ctr + 1] = Int32(2) - Int32(t1)
                ctr += 1
            end
        else  # eta == 4
            # C lines 407-412
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

# poly.c lines 435-453
# Fill a[1:256] with coefficients in [-eta, eta] via SHAKE256(seed || nonce)
# with re-squeeze loop.
function poly_uniform_eta!(a::Vector{Int32}, seed::Vector{UInt8},
                           nonce::UInt16, eta::Int)
    # C lines 430-434: initial buffer size depends on eta
    nblocks_init = eta == 2 ? 1 : 2
    buflen = nblocks_init * SHAKE256_RATE

    input = vcat(seed, UInt8(nonce & 0xFF), UInt8((nonce >> 8) & 0xFF))

    # Track total bytes squeezed for deterministic re-squeeze
    bytes_squeezed = 0

    function squeeze_blocks(nblocks::Int)::Vector{UInt8}
        start  = bytes_squeezed + 1
        nbytes = nblocks * SHAKE256_RATE
        out    = SHA.shake256(input, UInt64(bytes_squeezed + nbytes))
        bytes_squeezed += nbytes
        return out[start:end]
    end

    # C line 445: initial squeeze
    buf = squeeze_blocks(nblocks_init)

    # C line 447: first rejection pass
    ctr = rej_eta(a, 0, 256, buf, buflen, eta)

    # C lines 449-452: re-squeeze loop
    while ctr < 256
        buf  = squeeze_blocks(1)
        ctr += rej_eta(a, ctr, 256 - ctr, buf, SHAKE256_RATE, eta)
    end

    return nothing
end