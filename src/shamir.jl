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

# GMP's Miller–Rabin (mpz_probab_prime_p): 0 composite, 1 probably prime, 2 prime.
isprobableprime(p::Integer) = p >= 2 &&
    ccall((:__gmpz_probab_prime_p, Base.GMP.libgmp), Cint, (Ref{BigInt}, Cint), BigInt(p), 40) > 0

# The field modulus must be prime, or interpolation has no inverses and sharing is not secret.
check_modulus(p::Integer) = p == PRIME_P || isprobableprime(p) ||
    throw(ArgumentError("modulus p must be prime, got $p"))

function shamir_share(secret::Integer, k::Int, n::Int, p::Integer=PRIME_P)
    k < 1 && throw(ArgumentError("Threshold k must be >= 1, got $k"))
    n < 1 && throw(ArgumentError("Number of shares n must be >= 1, got $n"))
    k > n && throw(ArgumentError("Threshold k=$k exceeds number of shares n=$n"))
    check_modulus(p)
    # Shares sit at x = 1..n; x = p would be f(0), the secret itself.
    n < p || throw(ArgumentError("Number of shares n=$n must be < p"))
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

function shamir_reconstruct(shares::Vector{Tuple{Int, T}}, k::Int, p::Integer=PRIME_P) where T <: Integer
    k < 1 && throw(ArgumentError("Threshold k must be >= 1, got $k"))
    length(shares) < k && throw(ArgumentError("Need >= $k shares, got $(length(shares))"))
    check_modulus(p)

    subset = view(shares, 1:k)          # no copy: Int shares would leave their values in plain data

    # Validate: no duplicate x-coordinates in the field (x and x + p are the same point)
    xs = [mod(big(s[1]), p) for s in subset]
    length(unique(xs)) == length(xs) || throw(ArgumentError("Duplicate share x-coordinates detected"))

    # Lagrange interpolation at x=0, over BigInt (Int differences of x-coordinates can overflow)
    result = big(0)
    for i in 1:length(subset)
        xi, yi = xs[i], subset[i][2]
        num, den = big(1), big(1)
        for j in 1:length(subset)
            i == j && continue
            xj = xs[j]
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
function shamir_share_bytes(secret_bytes::Vector{UInt8}, k::Int, n::Int, p::Integer=PRIME_P)
    isempty(secret_bytes) && throw(ArgumentError("secret must be at least 1 byte"))
    # Big-endian bytes straight into one BigInt (mpz_import), with no hex String copy.
    secret_int = BigInt()
    GC.@preserve secret_bytes ccall((:__gmpz_import, Base.GMP.libgmp), Cvoid,
        (Ref{BigInt}, Csize_t, Cint, Csize_t, Cint, Csize_t, Ptr{UInt8}),
        secret_int, length(secret_bytes), 1, 1, 1, 0, secret_bytes)
    return shamir_share(secret_int, k, n, p)
end

"""
    shamir_reconstruct_bytes(shares, k, nbytes) -> Vector{UInt8}

Reconstruct raw bytes from Shamir shares. `nbytes` is the expected output
length (e.g., 32 for a 256-bit key).
"""
function shamir_reconstruct_bytes(shares::Vector{Tuple{Int, T}}, k::Int, nbytes::Int, p::Integer=PRIME_P) where T <: Integer
    nbytes >= 1 || throw(ArgumentError("nbytes must be >= 1, got $nbytes"))
    secret_int = shamir_reconstruct(shares, k, p)
    # A value of more than nbytes bytes means inconsistent shares or a wrong nbytes.
    used = iszero(secret_int) ? 0 : cld(ndigits(secret_int; base=2), 8)
    used <= nbytes ||
        throw(ArgumentError("reconstructed secret does not fit in $nbytes bytes (wrong shares or nbytes)"))
    # Big-endian, right-aligned, as shamir_share_bytes reads it (mpz_export, no shifted copies).
    out = zeros(UInt8, nbytes)
    count = Ref{Csize_t}(0)
    GC.@preserve out ccall((:__gmpz_export, Base.GMP.libgmp), Ptr{Cvoid},
        (Ptr{UInt8}, Ref{Csize_t}, Cint, Csize_t, Cint, Csize_t, Ref{BigInt}),
        pointer(out, nbytes - used + 1), count, 1, 1, 1, 0, secret_int)
    count[] == used || error("mpz_export wrote $(count[]) bytes, expected $used")
    return out
end
