using QuantumAdiabaticAnnealing
using GenericTensorNetworks
using Graphs
using Profile

function annealing(graph, a, b, MIS_size, run_annealing_iters, tempscale = 15 .- (1:300 .- 1) .* 0.05, niters = 4000)
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

function init(begin_energy_gap, end_energy_gap, step, n, m)
    ret_weights = []
    ret_MIS_size = []
    ret_graph = []
    for energy_gap in begin_energy_gap:step:end_energy_gap
        locations, weights, inputs_id, outputs_id, input_layer_id = rule110_transverse_generate(n, m; gradient = 1:energy_gap:m*energy_gap)
        locations = map(t -> (Float64(t[1]), Float64(t[2])), locations)
        weights = Float64.(weights)

        # control_input!(weights, inputs_id, input_layer_id)
        control_output!(weights, outputs_id, input_layer_id)

        graph = unit_disk_graph(locations, 2.05)
        problem = GenericTensorNetwork(IndependentSet(graph, weights));
        max_independent_set = solve(problem, SizeMax(1))[]
        MIS_size = max_independent_set.orders[1].n  

        push!(ret_graph, graph)
        push!(ret_MIS_size, MIS_size)
        push!(ret_weights, weights)
    end
    return ret_graph, ret_weights, ret_MIS_size
end
n = 7
m = 7

graph_arr, weights_arr, MIS_size_arr = init(1.0, 3.1, 0.3, n, m)



# Now try Simulated annealing
total_success_time = zeros(Int, length(MIS_size_arr))

begin_time = time()
for stage in 1:2000
    for energy_id in 1:1
        graph = graph_arr[energy_id]
        detuning = - weights_arr[energy_id]
        MIS_size = MIS_size_arr[energy_id]

        success_time = annealing(graph, detuning, Inf, MIS_size, 1)
        total_success_time[energy_id] += success_time
        @info "stage = $stage, energy_id= $energy_id, total_success_time = $total_success_time"
    end
    @info "stage = $stage, total_success_time = $total_success_time"
end
end_time = time()
# colorconfig = [SA.best_IS_bitarr[i] + 1 for i in 1:nv(graph)]
# print_graph(graph, locations, colorconfig)
