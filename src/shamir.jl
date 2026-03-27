# Shamir (k,n)-threshold Secret Sharing over GF(2^127-1)

const PRIME_P = big(2)^127 - 1

function mod_inverse(a::Integer, m::Integer)::Integer
    a = mod(a, m)
    gcd(a, m) != 1 && throw(ArgumentError("No modular inverse"))
    old_r, r = a, m
    old_s, s = 1, 0
    while r != 0
        q = old_r ÷ r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
    end
    return mod(old_s, m)
end

function shamir_share(secret::Integer, k::Int, n::Int, p::Integer=PRIME_P)
    k > n && throw(ArgumentError("k > n"))
    (secret < 0 || secret >= p) && throw(ArgumentError("Secret out of range"))
    coeffs = [i == 1 ? BigInt(secret) : rand(big(0):(p-1)) for i in 1:k]
    shares = [(i, mod(sum(coeffs[j] * powermod(big(i), j-1, p) for j in 1:k), p)) for i in 1:n]
    return shares
end

function shamir_reconstruct(shares::Vector{Tuple{Int, T}}, k::Int, p::Integer=PRIME_P) where T <: Integer
    length(shares) < k && throw(ArgumentError("Insufficient shares"))
    subset = shares[1:k]
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
