using QuantumAdiabaticAnnealing, Test
using QuantumAdiabaticAnnealing:rule54_generate
using GenericTensorNetworks
using LinearAlgebra

@testset "rule54_generate" begin
    graph, locations, weights, P, Q, R, Target = rule54_generate()
    # locations = map(t->t.*100, locations)
    # locations = map(t->(t[1], -t[2]), locations)

    # test if the graph could give the correct output with pined nodes
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
    @test max_config_weighted.c.data[Target] == 0

    weights[P] = 50
    weights[Q] = -50
    weights[R] = -50

    @info "Test 100"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 1

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
    @test max_config_weighted.c.data[Target] == 0

    weights[P] = 50
    weights[Q] = 50
    weights[R] = 50

    @info "Test 111"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[Target] == 0


    # Test the graph is a unit disk graph with radius sqrt(5)
    Max_radius_square = 5
    for i in 1:Target
        for j in i+1:Target
            @info "now testing edge $i $j"
            if has_edge(graph, i, j)
                @test (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 <= Max_radius_square
            else
                @test (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 > Max_radius_square
            end
        end
    end
end