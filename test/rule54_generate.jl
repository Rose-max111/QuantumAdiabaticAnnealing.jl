using QuantumAdiabaticAnnealing, Test
using QuantumAdiabaticAnnealing:rule54_generate
using GenericTensorNetworks
using LinearAlgebra

@testset "rule54_generate" begin
    graph, locations, weights, P, Q, R = rule54_generate()
    locations = map(t->t.*100, locations)
    locations = map(t->(t[1], -t[2]), locations)

    weights[P] = -50
    weights[Q] = -50
    weights[R] = -50

    @info "Test 000"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 0

    weights[P] = -50
    weights[Q] = -50
    weights[R] = 50

    @info "Test 001"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 1

    weights[P] = -50
    weights[Q] = 50
    weights[R] = -50

    @info "Test 010"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 1

    weights[P] = -50
    weights[Q] = 50
    weights[R] = 50

    @info "Test 011"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 0

    weights[P] = 50
    weights[Q] = -50
    weights[R] = -50

    @info "Test 100"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 1

    weights[P] = 50
    weights[Q] = -50
    weights[R] = 50

    @info "Test 101"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 1

    weights[P] = 50
    weights[Q] = 50
    weights[R] = -50

    @info "Test 110"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 0

    weights[P] = 50
    weights[Q] = 50
    weights[R] = 50

    @info "Test 111"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 0
end