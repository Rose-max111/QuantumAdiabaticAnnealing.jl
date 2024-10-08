using QuantumAdiabaticAnnealing
using GenericTensorNetworks
using Graphs
using CairoMakie

function annealing(graph, a, b, MIS_size, run_annealing_iters, tempscale = 6 .- (1:30 .- 1) .* 0.2, niters = 500)
    success_time = 0
    SA = SimulatedAnnealingMIS(graph; a = a, b = b)
    for i in 1:run_annealing_iters
        SA = SimulatedAnnealingMIS(graph; a = a, b = b)
        track_equilibration!(SA, -(MIS_size), 1000, tempscale, niters)
        # @info "runtime = $i, this_time_best_obj = $(SA.best_obj), success_time = $(success_time)"
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
        weights[i] = 1000
    end
    for i in input_layer_id
        weights[i[1]] = weights[i[2]] = 1000
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

function evaluate_slice_performance(n, energy_gap)
    m = 1
    locations, weights, inputs_id, outputs_id, input_layer_id = rule110_transverse_generate(n, m; gradient = [energy_gap^i for i in m-1:-1:0])
    locations = map(t -> (Float64(t[1]), Float64(t[2])), locations)
    weights = Float64.(weights)

    control_input!(weights, inputs_id, input_layer_id)
    graph = unit_disk_graph(locations, 2.05)

    # First use TensorNetwork method to get the size of weighted MIS problem   
    problem = GenericTensorNetwork(IndependentSet(graph, weights));
    max_independent_set = solve(problem, SizeMax(1))[]

    MIS_size = max_independent_set.orders[1].n

    # Now try Simulated annealing

    prob = []
    run_time = 1000:500:21000
    for simulate_time in run_time
        Temp_max = 10.0
        tempscale = Temp_max .- (1:simulate_time .- 1) .* (Temp_max / simulate_time)
        Real_tempscale = [tempscale[i] * energy_gap^j for j in m-1:-1:0 for i in 1:length(tempscale)]
        niters = 1

        detuning = -weights
        success_time = annealing(graph, detuning, Inf, MIS_size, 2000, Real_tempscale, niters)
        @info "simulate_time = $simulate_time, success_time = $success_time, success_prob = $(1.0 * success_time / 2000)"
        push!(prob, 1.0 * success_time / 2000)
    end
    return prob
end


prob_3_gadgets = evaluate_slice_performance(3, 1)

filepath = joinpath(@__DIR__, "data/prob_3_gadgets.txt")
open(filepath,"w") do file
    for i in prob_3_gadgets
        println(file, i)
    end
end
