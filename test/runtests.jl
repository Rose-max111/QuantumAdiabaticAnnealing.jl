using QuantumAdiabaticAnnealing
using Test

@testset "generate graphs" begin
    include("generate_graph.jl")
end

# @testset "energy_calculate.jl" begin
#     include("energy_calculate.jl")
# end

@testset "cuda" begin
    include("cusa.jl")
end