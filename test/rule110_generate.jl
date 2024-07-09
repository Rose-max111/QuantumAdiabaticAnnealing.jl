using QuantumAdiabaticAnnealing, Test
using QuantumAdiabaticAnnealing:rule110_generate
using GenericTensorNetworks
using LinearAlgebra
using Graphs

@testset "rule110_generate" begin
    graph, weights, locations, P, Q, R, Target = rule110_generate()

    weights[P] = -50
    weights[Q] = -50
    weights[R] = -50

    @info "Test 000"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 0

    weights[P] = -50
    weights[Q] = -50
    weights[R] = 50

    @info "Test 001"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = -50
    weights[Q] = 50
    weights[R] = -50

    @info "Test 010"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = -50
    weights[Q] = 50
    weights[R] = 50

    @info "Test 011"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = 50
    weights[Q] = -50
    weights[R] = -50

    @info "Test 100"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 0

    weights[P] = 50
    weights[Q] = -50
    weights[R] = 50

    @info "Test 101"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = 50
    weights[Q] = 50
    weights[R] = -50

    @info "Test 110"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = 50
    weights[Q] = 50
    weights[R] = 50

    @info "Test 111"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 0

    Max_radius_square = 4.2
    for i in 1:nv(graph)
        for j in i+1:nv(graph)
            @info "now testing edge $i $j"
            if has_edge(graph, i, j)
                @test (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 <= Max_radius_square
            else
                @test (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 > Max_radius_square
            end
        end
    end
end