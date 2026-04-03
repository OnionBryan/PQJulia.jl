#!/usr/bin/env julia
"""
PQJulia.jl Benchmark Suite
===========================
33 benchmarks: 6 ML-KEM ops × 3 levels + 5 ML-DSA ops × 3 levels
Uses BenchmarkTools @benchmark for statistical rigor (min/median/mean/max).
"""

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using PQJulia
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
BenchmarkTools.DEFAULT_PARAMETERS.samples = 50

println("=" ^ 70)
println("PQJulia.jl Benchmark Suite")
println("=" ^ 70)

results = Dict{String, Any}()

# ═══════════════════════════════════════════════════════════════
# ML-KEM Benchmarks
# ═══════════════════════════════════════════════════════════════
for (name, Cat) in [
    ("ML-KEM-512", MLKEM.Category1),
    ("ML-KEM-768", MLKEM.Category3),
    ("ML-KEM-1024", MLKEM.Category5)
]
    println("\n--- $name ---")

    # KeyGen
    b = @benchmark $Cat.kyber_kem_keypair()
    results["$(name)/KeyGen"] = b
    println("  KeyGen:  $(BenchmarkTools.prettytime(median(b).time))  ($(length(b.times)) samples)")

    # Setup for Encaps/Decaps
    pk, sk = Cat.kyber_kem_keypair()

    # Encaps
    b = @benchmark $Cat.kyber_kem_enc($pk)
    results["$(name)/Encaps"] = b
    println("  Encaps:  $(BenchmarkTools.prettytime(median(b).time))")

    # Decaps
    ct, _ = Cat.kyber_kem_enc(pk)
    b = @benchmark $Cat.kyber_kem_dec($ct, $sk)
    results["$(name)/Decaps"] = b
    println("  Decaps:  $(BenchmarkTools.prettytime(median(b).time))")
end

# ═══════════════════════════════════════════════════════════════
# ML-DSA Benchmarks
# ═══════════════════════════════════════════════════════════════
msg = Vector{UInt8}("Benchmark message for ML-DSA signing")

for (name, Cat) in [
    ("ML-DSA-44", MLDSA.Category2),
    ("ML-DSA-65", MLDSA.Category3),
    ("ML-DSA-87", MLDSA.Category5)
]
    println("\n--- $name ---")

    # KeyGen
    b = @benchmark $Cat.dilithium_keygen()
    results["$(name)/KeyGen"] = b
    println("  KeyGen:  $(BenchmarkTools.prettytime(median(b).time))")

    # Setup
    pk, sk = Cat.dilithium_keygen()

    # Sign
    b = @benchmark $Cat.dilithium_sign($msg, $sk)
    results["$(name)/Sign"] = b
    println("  Sign:    $(BenchmarkTools.prettytime(median(b).time))")

    # Verify
    sig = Cat.dilithium_sign(msg, sk)
    b = @benchmark $Cat.dilithium_verify($msg, $sig, $pk)
    results["$(name)/Verify"] = b
    println("  Verify:  $(BenchmarkTools.prettytime(median(b).time))")
end

# ═══════════════════════════════════════════════════════════════
# Shamir Benchmarks
# ═══════════════════════════════════════════════════════════════
println("\n--- Shamir (3,5) ---")
secret = parse(BigInt, bytes2hex(rand(UInt8, 32)), base=16)

b = @benchmark shamir_share($secret, 3, 5)
results["Shamir/Share(3,5)"] = b
println("  Share:       $(BenchmarkTools.prettytime(median(b).time))")

shares = shamir_share(secret, 3, 5)
b = @benchmark shamir_reconstruct($([shares[1], shares[3], shares[5]]), 3)
results["Shamir/Reconstruct(3,5)"] = b
println("  Reconstruct: $(BenchmarkTools.prettytime(median(b).time))")

# ═══════════════════════════════════════════════════════════════
# Summary Table
# ═══════════════════════════════════════════════════════════════
println("\n" * "=" ^ 70)
println("SUMMARY (median times)")
println("=" ^ 70)
println(rpad("Operation", 30) * rpad("Median", 15) * "Allocs")
println("-" ^ 60)
for key in sort(collect(keys(results)))
    b = results[key]
    m = median(b)
    println(rpad(key, 30) * rpad(BenchmarkTools.prettytime(m.time), 15) * string(m.allocs))
end
println("=" ^ 70)
