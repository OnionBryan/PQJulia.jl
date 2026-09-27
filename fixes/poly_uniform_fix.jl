# Fix for Critical #1: poly_uniform! missing re-squeeze loop
# Source: dilithium/ref/poly.c lines 295-368
# The existing implementation uses a fixed 1280-byte buffer and error() on exhaustion.
# The C reference uses a while(ctr < N) re-squeeze loop, squeezing one SHAKE128 block
# at a time until all 256 coefficients are filled.

using SHA

const SHAKE128_RATE = 168          # fips202.h; STREAM128_BLOCKBYTES = SHAKE128_RATE
const DIL_N         = 256          # params.h
const DIL_Q         = Int32(8380417)  # params.h

# poly.c lines 309-331
# Sample uniformly random coefficients in [0, Q-1] by rejection sampling on buf.
# Returns count of accepted coefficients; may be < len if buf is exhausted.
function rej_uniform(a::Vector{Int32}, len::Int, buf::Vector{UInt8}, buflen::Int)::Int
    ctr = 0
    pos = 1  # 1-indexed (C: pos = 0)
    # C line 319: while(ctr < len && pos + 3 <= buflen)
    while ctr < len && pos + 2 <= buflen
        # C lines 320-323: assemble 23-bit value from 3 bytes LE
        t  = UInt32(buf[pos])
        t |= UInt32(buf[pos + 1]) << 8
        t |= UInt32(buf[pos + 2]) << 16
        t &= UInt32(0x7FFFFF)
        pos += 3
        # C lines 325-326: reject if >= Q
        if t < UInt32(DIL_Q)
            ctr += 1
            a[ctr] = Int32(t)
        end
    end
    return ctr
end

# poly.c lines 344-368
# Sample polynomial with uniformly random coefficients in [0, Q-1] via SHAKE128(seed || nonce).
# Uses re-squeeze loop until all N=256 coefficients are filled.
function poly_uniform!(a::Vector{Int32}, seed::Vector{UInt8}, nonce::UInt16)
    # stream128_init absorbs seed || nonce_le16
    input = vcat(seed, UInt8[nonce % UInt8, (nonce >> 8) % UInt8])

    # C line 344: POLY_UNIFORM_NBLOCKS = ceil(768/168) = 5
    nblocks   = (768 + SHAKE128_RATE - 1) ÷ SHAKE128_RATE   # = 5
    buflen    = nblocks * SHAKE128_RATE                        # = 840
    total_out = buflen

    # C line 355: stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state)
    buf = SHA.shake128(input, UInt64(total_out))

    # C line 357: ctr = rej_uniform(a->coeffs, N, buf, buflen)
    ctr = rej_uniform(a, DIL_N, buf, buflen)

    # C lines 359-367: re-squeeze loop
    while ctr < DIL_N
        # C line 360: off = buflen % 3
        off = buflen % 3

        # C line 364: squeeze one more block
        new_total = total_out + SHAKE128_RATE
        full_stream = SHA.shake128(input, UInt64(new_total))

        # C lines 361-362: carry trailing bytes
        new_buf = Vector{UInt8}(undef, SHAKE128_RATE + off)
        for i in 1:off
            new_buf[i] = buf[buflen - off + i]
        end
        for i in 1:SHAKE128_RATE
            new_buf[off + i] = full_stream[total_out + i]
        end

        buflen    = SHAKE128_RATE + off
        total_out = new_total
        buf       = new_buf

        # C line 366: ctr += rej_uniform(a->coeffs + ctr, N - ctr, buf, buflen)
        ctr += rej_uniform(view(a, ctr + 1 : length(a)), DIL_N - ctr, buf, buflen)
    end
    return nothing
end