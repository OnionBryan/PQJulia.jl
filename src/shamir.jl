# Shamir (k,n)-threshold Secret Sharing
# Prime must be > any secret value. 2^521-1 (13th Mersenne prime) handles
# secrets up to 520 bits (covers 256-bit keys, 512-bit hashes, etc.)

import Random

const PRIME_P = big(2)^521 - 1

function mod_inverse(a::Integer, m::Integer)::Integer
    a = mod(a, m)
    gcd(a, m) != 1 && throw(ArgumentError("No modular inverse (a=$a, m=$m)"))
    old_r, r = a, m
    old_s, s = 1, 0
    while r != 0
        q = old_r ÷ r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
    end
    return mod(old_s, m)
end

function shamir_share(secret::Integer, k::Int, n::Int)
    p = PRIME_P
    k < 1 && throw(ArgumentError("Threshold k must be >= 1, got $k"))
    n < 1 && throw(ArgumentError("Number of shares n must be >= 1, got $n"))
    k > n && throw(ArgumentError("Threshold k=$k exceeds number of shares n=$n"))
    (secret < 0 || secret >= p) && throw(ArgumentError("Secret out of range [0, p)"))

    # Build random polynomial: f(x) = secret + a1*x + ... + a_{k-1}*x^{k-1}
    coeffs = Vector{BigInt}(undef, k)
    coeffs[1] = BigInt(secret)
    for j in 2:k
        coeffs[j] = rand(Random.RandomDevice(), big(0):(p-1))
    end
    # Ensure leading coefficient is nonzero to preserve threshold
    # (rwot8 best practice: zero leading coeff reduces threshold by 1)
    if k >= 2 && coeffs[k] == 0
        coeffs[k] = rand(Random.RandomDevice(), big(1):(p-1))
    end

    # Evaluate polynomial at x=1..n using Horner's method
    shares = Vector{Tuple{Int, BigInt}}(undef, n)
    for i in 1:n
        x = big(i)
        fx = big(0)
        for j in k:-1:1
            fx = mod(fx * x + coeffs[j], p)
        end
        shares[i] = (i, fx)
    end
    return shares
end

function shamir_reconstruct(shares::Vector{Tuple{Int, T}}, k::Int) where T <: Integer
    p = PRIME_P
    length(shares) < k && throw(ArgumentError("Need >= $k shares, got $(length(shares))"))

    subset = shares[1:k]

    # Validate: no duplicate x-coordinates
    xs = [s[1] for s in subset]
    length(unique(xs)) == length(xs) || throw(ArgumentError("Duplicate share x-coordinates detected"))

    # Lagrange interpolation at x=0
    result = big(0)
    for i in 1:length(subset)
        xi, yi = subset[i]
        num, den = big(1), big(1)
        for j in 1:length(subset)
            i == j && continue
            xj, _ = subset[j]
            num = mod(num * (0 - xj), p)
            den = mod(den * (xi - xj), p)
        end
        result = mod(result + yi * num * mod_inverse(den, p), p)
    end
    return result
end

# ── Byte-level API ────────────────────────────────────────────────────────

"""
    shamir_share_bytes(secret_bytes, k, n) -> Vector{Tuple{Int, BigInt}}

Split raw bytes into Shamir shares. Convenience wrapper that handles
BigInt encoding/decoding for byte-oriented secret material (keys, hashes).
"""
function shamir_share_bytes(secret_bytes::Vector{UInt8}, k::Int, n::Int)
    secret_int = parse(BigInt, bytes2hex(secret_bytes), base=16)
    return shamir_share(secret_int, k, n)
end

"""
    shamir_reconstruct_bytes(shares, k, nbytes) -> Vector{UInt8}

Reconstruct raw bytes from Shamir shares. `nbytes` is the expected output
length (e.g., 32 for a 256-bit key).
"""
function shamir_reconstruct_bytes(shares::Vector{Tuple{Int, T}}, k::Int, nbytes::Int) where T <: Integer
    secret_int = shamir_reconstruct(shares, k)
    hex = string(secret_int, base=16)
    # Pad to expected length
    hex = lpad(hex, 2 * nbytes, '0')
    return hex2bytes(hex)
end
