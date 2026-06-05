using BenchmarkTools
include("../src/PQJulia.jl")
using .PQJulia

# Function to run benchmark and return median time in ms
function run_bench()
    b2 = @benchmark PQJulia.MLDSA.Category2.dilithium_keygen() samples=100
    b3 = @benchmark PQJulia.MLDSA.Category3.dilithium_keygen() samples=100
    b5 = @benchmark PQJulia.MLDSA.Category5.dilithium_keygen() samples=100

    return (
        median(b2.times) / 1e6,
        median(b3.times) / 1e6,
        median(b5.times) / 1e6
    )
end

println("Warmup...")
PQJulia.MLDSA.Category2.dilithium_keygen()
PQJulia.MLDSA.Category3.dilithium_keygen()
PQJulia.MLDSA.Category5.dilithium_keygen()

println("Running benchmark...")
t2, t3, t5 = run_bench()
println("Category 2: $t2 ms")
println("Category 3: $t3 ms")
println("Category 5: $t5 ms")
