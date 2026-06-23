using Pkg
Pkg.activate(joinpath(@__DIR__, "."))
using PQJulia
using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 5
BenchmarkTools.DEFAULT_PARAMETERS.samples = 100

seed = rand(UInt8, 32)
b = @benchmark PQJulia.MLKEM.Category1.kyber_gen_matrix($seed, false, 2)
println("Original Category 1 kyber_gen_matrix: ", median(b))
b = @benchmark PQJulia.MLKEM.Category3.kyber_gen_matrix($seed, false, 3)
println("Original Category 3 kyber_gen_matrix: ", median(b))
b = @benchmark PQJulia.MLKEM.Category5.kyber_gen_matrix($seed, false, 4)
println("Original Category 5 kyber_gen_matrix: ", median(b))
