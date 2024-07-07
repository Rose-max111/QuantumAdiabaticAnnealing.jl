using QuantumAdiabaticAnnealing, Test
using QuantumAdiabaticAnnealing:rule54_generate

@test "rule54_generate" begin
    graph, locations, weights, P, Q, R = rule54_generate()
    locations = map(t->t.*100, locations)
    locations = map(t->(t[1], -t[2]), locations)

    weights[P] = -50
    weights[Q] = -50
    weights[R] = -50

    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    @test max_config_weighted.c.data[60] == 0
end