using QuantumAdiabaticAnnealing
using GenericTensorNetworks
using Graphs

function annealing(graph, a, b, MIS_size, run_annealing_iters, tempscale = 6 .- (1:200 .- 1) .* 0.03, niters = 4000)
    success_time = 0
    SA = SimulatedAnnealingMIS(graph; a = a, b = b)
    for i in 1:run_annealing_iters
        SA = SimulatedAnnealingMIS(graph; a = a, b = b)
        track_equilibration!(SA, -(MIS_size), 1000, tempscale, niters)
        @info "runtime = $i, this_time_best_obj = $(SA.best_obj), success_time = $(success_time)"
        if abs(SA.best_obj + (MIS_size)) <= 1e-2
            success_time += 1
        end
    end
    return success_time
end

function print_graph(G, locations, colorconfig)
    colorspace = ["Black", "Red"]
    vertex_colors = colorspace[colorconfig]
    show_graph(G, locs = locations; format = :svg, vertex_colors = vertex_colors, texts = map(string, 1:nv(G)), vertex_text_colors = fill("White", nv(G)))
end

function control_input!(weights, inputs_id, input_layer_id)
    for i in inputs_id
        weights[i] = 100
    end
    for i in input_layer_id
        weights[i[1]] = weights[i[2]] = 100
    end
end

function control_output!(weights, outputs_id, input_layer_id)
    for i in outputs_id
        weights[i] = 100
    end
    for i in input_layer_id
        weights[i[1]] = weights[i[2]] = 100
    end
end

n = 7
m = 7
energy_gap = 1.4
# energy_gap = ARGS[1]
# energy_gap = parse(Float64, energy_gap)

locations, weights, inputs_id, outputs_id, input_layer_id = rule110_transverse_generate(n, m; gradient = m*energy_gap:-energy_gap:1)
locations = map(t -> (Float64(t[1]), Float64(t[2])), locations)
weights = Float64.(weights)

control_input!(weights, inputs_id, input_layer_id)
# control_output!(weights, outputs_id, input_layer_id)

graph = unit_disk_graph(locations, 2.05)

# First use TensorNetwork method to get the size of weighted MIS problem
problem = GenericTensorNetwork(IndependentSet(graph, weights));
max_independent_set = solve(problem, SizeMax(1))[]
# max_independent_set = solve(problem, ConfigsMax())[]
# colorconfig = [max_independent_set.c.data[1][i] == 0 ? 1 : 2 for i in 1:nv(graph)]
# print_graph(graph, locations, colorconfig)

MIS_size = max_independent_set.orders[1].n

# Now try Simulated annealing
detuning = -weights
success_time = annealing(graph, detuning, Inf, MIS_size, 20000)
println("success_time = ", success_time, "energy_gap = ", energy_gap)

# colorconfig = [SA.best_IS_bitarr[i] + 1 for i in 1:nv(graph)]
# print_graph(graph, locations, colorconfig)
