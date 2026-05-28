using Pkg
Pkg.activate(".")

include("src/PQJulia.jl")
using .PQJulia
using BenchmarkTools

# Setup for kyber_poly_decompress!
# d=5: 160 bytes -> 256 coefficients
const KYBER_N = 256

function run_bench()
    a_bytes = rand(UInt8, 160)
    r = Vector{Int16}(undef, 256)

    wrapper() = PQJulia.MLKEM.kyber_poly_decompress!(r, a_bytes, 5)

    b = @benchmark $wrapper()
    display(b)
    println()
end

run_bench()
