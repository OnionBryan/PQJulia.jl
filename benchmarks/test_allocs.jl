using BenchmarkTools
include("../src/PQJulia.jl")
using .PQJulia

println("Setting up benchmark for Dilithium keygen Category 2")
b2 = @benchmark PQJulia.MLDSA.Category2.dilithium_keygen()
display(b2)
println()

println("Setting up benchmark for Dilithium keygen Category 3")
b3 = @benchmark PQJulia.MLDSA.Category3.dilithium_keygen()
display(b3)
println()

println("Setting up benchmark for Dilithium keygen Category 5")
b5 = @benchmark PQJulia.MLDSA.Category5.dilithium_keygen()
display(b5)
println()
