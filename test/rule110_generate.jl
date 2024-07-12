using QuantumAdiabaticAnnealing, Test
using QuantumAdiabaticAnnealing:rule110_generate, transversal_graph, rule110_transverse_generate
using GenericTensorNetworks
using LinearAlgebra
using Graphs
# using LuxorGraphPlot.Luxor, LuxorGraphPlot

@testset "rule110_generate" begin
    graph, weights, locations, P, Q, R, Target = rule110_generate()

    weights[P] = -50
    weights[Q] = -50
    weights[R] = -50

    @info "Test 000"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    max_independent_set_size = solve(problem, ConfigsMax())[]
    @info "max_independent_set_size = $max_independent_set_size"
    @test max_config_weighted.c.data[Target] == 0


    weights[P] = -50
    weights[Q] = -50
    weights[R] = 50

    @info "Test 001"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    max_independent_set_size = solve(problem, ConfigsMax())[]
    @info "max_independent_set_size = $max_independent_set_size"
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = -50
    weights[Q] = 50
    weights[R] = -50

    @info "Test 010"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    max_independent_set_size = solve(problem, ConfigsMax())[]
    @info "max_independent_set_size = $max_independent_set_size"
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = -50
    weights[Q] = 50
    weights[R] = 50

    @info "Test 011"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    max_independent_set_size = solve(problem, ConfigsMax())[]
    @info "max_independent_set_size = $max_independent_set_size"
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = 50
    weights[Q] = -50
    weights[R] = -50

    @info "Test 100"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    max_independent_set_size = solve(problem, ConfigsMax())[]
    @info "max_independent_set_size = $max_independent_set_size"
    @test max_config_weighted.c.data[Target] == 0

    weights[P] = 50
    weights[Q] = -50
    weights[R] = 50

    @info "Test 101"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    max_independent_set_size = solve(problem, ConfigsMax())[]
    @info "max_independent_set_size = $max_independent_set_size"
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = 50
    weights[Q] = 50
    weights[R] = -50

    @info "Test 110"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    max_independent_set_size = solve(problem, ConfigsMax())[]
    @info "max_independent_set_size = $max_independent_set_size"
    @test max_config_weighted.c.data[Target] == 1

    weights[P] = 50
    weights[Q] = 50
    weights[R] = 50

    @info "Test 111"
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_config_weighted = solve(problem, SingleConfigMax())[]
    max_independent_set_size = solve(problem, ConfigsMax())[]
    @info "max_independent_set_size = $max_independent_set_size"
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

@testset "transversal_graph" begin
    locations, weights, input, output = transversal_graph()
    radius_square = 4.2

    el = Edge.([(1,2),
            (2,3),
            (3,4),
            (4,6),
            (4,8),
            (5,6),
            (5,7),
            (6,7),
            (6,8),
            (6,9),
            (7,8),
            (7,9),
            (8,9),
            (7,10),
            (9,10),
            (8,11),
            (9,11),
            (11,18),
            (2,12),
            (12,15),
            (12,16),
            (12,17),
            (13,14),
            (13,15),
            (14,15),
            (14,16),
            (15,16),
            (15,17),
            (16,17),
            (17,18),
            (17,19),
            (16,19),
            (18,19),
            (16,20),
            (19,20),
            (15,21),
            (20,21),
            (21,22),
            (22,23)])
    G = SimpleGraph(el)
    for i in 1:length(locations)
        for j in i+1:length(locations)
            @info "now testing edge $i $j"
            if (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 <= radius_square
                @test has_edge(G, i, j)
            else
                @test !has_edge(G, i, j)
            end
        end
    end
end

# @testset "rule110_transverse_generate" begin
#     locations, weights, input_id, output_id = rule110_transverse_generate(3,2)
#     G = SimpleGraph(length(weights))
#     radius_square = 4.2
#     for i in 1:length(locations)
#         for j in i+1:length(locations)
#             if (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 <= radius_square
#                 add_edge!(G, i, j)
#             end
#         end
#     end
    
#     target_output = [0, 1, 1, 1, 0, 1, 1, 0]
#     for input_mask in 0:7
        
# end

function show_transversal_graph(n, m)
    locations, weights, input_id, output_id, input_layer_id = rule110_transverse_generate(n,m)
    G = SimpleGraph(length(weights))
    radius_square = 4.2
    for i in 1:length(locations)
        for j in i+1:length(locations)
            if (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 <= radius_square
                add_edge!(G, i, j)
            end
        end
    end

    colorspace = ["Blue", "Red", "Green", "Black"]

    show_graph(G, locs = locations; format = :svg, vertex_colors = colorspace[weights], vertexlabels = 1:nv(G))
end

# n=2
# m=2
# locations, weights, input_id, output_id, input_layer_id = rule110_transverse_generate(n,m)
# G = SimpleGraph(length(weights))
# radius_square = 4.2
# for i in 1:length(locations)
#     for j in i+1:length(locations)
#         if (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 <= radius_square
#             add_edge!(G, i, j)
#         end
#     end
# end
# problem = GenericTensorNetwork(IndependentSet(G, weights));
# max_independent_set_size = solve(problem, ConfigsMax())[]
# colorspace = ["Blue", "Red", "Green", "Black"]

# show_graph(G, locs = locations; format = :svg, vertex_colors = colorspace[weights], show_number = true)