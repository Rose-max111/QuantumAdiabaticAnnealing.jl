using QuantumAdiabaticAnnealing, Test
using GenericTensorNetworks
using LinearAlgebra
using Graphs
using LuxorGraphPlot

@testset "mis rules" begin
    for ridx in (54, 110)
        rule = CellularAutomata1D{ridx}()
        g = mis_gadget(rule)
        locations, weights, (P, Q, R), Target = g.coordinates, g.weights, g.inputs, g.outputs[]
        graph = topology(g)
        configs = solve(GenericTensorNetwork(IndependentSet(graph, weights)), ConfigsMax())[].c.data
        @test length(unique!(getindex.(configs, Ref([P, Q, R])))) == 8
        for c in configs
            input = c[[P, Q, R]]
            input_idx = input[3] + 2*input[2] + 4*input[1]
            @test c[Target] == ((ridx >> input_idx) & 1)
        end
    end
end

@testset "transversal_graph" begin
    g = QuantumAdiabaticAnnealing.transversal_mis_gadget(CellularAutomata1D{110}())
    locations = g.coordinates

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
            if (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 <= g.radius^2
                @test has_edge(G, i, j)
            else
                @test !has_edge(G, i, j)
            end
        end
    end
end


@testset "rule110_transverse_generate" begin
    rule = CellularAutomata1D{110}()
    locations, weights, inputs_id, outputs_id, input_layer_id = transversal_mis_gadget(rule, 1, 1)
    P = input_layer_id[1][1]
    Q = inputs_id[1]
    R = input_layer_id[1][2]
    Target = outputs_id[1]
    locations = map(t -> (Float64(t[1]), Float64(t[2])), locations)
    graph = unit_disk_graph(locations, 2.05)

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

# n=1
# m=3
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
# independent_set_size = solve(problem, CountingAll())[]
# colorspace = ["Blue", "Red", "Green", "Black"]

# show_graph(G, locs = locations; format = :svg, vertex_colors = colorspace[weights], show_number = true)


