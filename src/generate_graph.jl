using Bloqade
using UnitDiskMapping
using Graphs
using KrylovKit
using CairoMakie

function countbits(u::Int)
    if u == 0
        return 0
    end
    c = 0
    while u != 0 
        u = u >> 1
        c += 1
    end
    return c
end

function Distance(u, v)
    return sqrt((u[1] - v[1])^2 + (u[2] - v[2])^2)
end

function generate_factoring_new_graph(x1, x2)
    mres = UnitDiskMapping.map_factoring(x1, x2)
    vaild = ones(length(mres.grid_graph.nodes))
    radius = mres.grid_graph.radius
    
    # remove all pins_zeros nodes
    vaild[mres.pins_zeros] .= 0

    output_pins_zeros = []
    output_pins_ones = []
    answer = x1 * x2
    bit_of_answer = length(mres.pins_output)
    for i in 1:bit_of_answer
        if answer & (1 << (i-1)) != 0
            push!(output_pins_ones, mres.pins_output[i])
        else
            push!(output_pins_zeros, mres.pins_output[i])
        end
    end

    # remove all output_pins_zeros nodes
    vaild[output_pins_zeros] .= 0

    # remove all output_pins_ones nodes
    vaild[output_pins_ones] .= 0

    # remove all nodes that are too close to output_pins_ones
    for i in 1:length(mres.grid_graph.nodes)
        for pin_one in output_pins_ones
            if Distance(mres.grid_graph.nodes[i].loc, mres.grid_graph.nodes[pin_one].loc) < radius
                vaild[i] = 0
            end
        end
    end
    my = mres.grid_graph.size[2]
    mx = mres.grid_graph.size[1]
    new_graph_nodes = [tuple(float(mres.grid_graph.nodes[t].loc[2]), float(mx - mres.grid_graph.nodes[t].loc[1] + 1)) for t in 1:length(vaild) if vaild[t] == 1.0]
    new_graph_weights = [mres.grid_graph.nodes[t].weight for t in 1:length(vaild) if vaild[t] == 1.0]
    return new_graph_nodes, new_graph_weights, radius
end



# new_graph_nodes, new_graph_weights, radius = generate_factoring_new_graph(1, 1)
# atoms = AtomList(new_graph_nodes)

function generate_some_graph()

    origin_graph = [(4,2),
    (2,6),
    (6,8),
    (10,6),
    (8,2),
    (0,6),
    (14,6),
    (6,12),
    (20,6),
    (24,6),
    (28,8),
    (32,6),
    (36,6),
    (30,2),
    (26,2),
    (28,12),
    (10,14),
    (26,14),
    (14,16),
    (16,12),
    (20,12),
    (22,16),
    (18,18)]

    weights = [1,
    2,
    2,
    2,
    1,
    1,
    1,
    2,
    1,
    2,
    2,
    2,
    1,
    1,
    1,
    2,
    2,
    2,
    2,
    1,
    1,
    2,
    1]


    vaild = ones(length(weights))
    vaild[[6, 2, 4, 7, 9, 13]] .= 0
    # vaild[end] = 0



    new_graph = [tuple(float(origin_graph[i][1]), float(origin_graph[i][2])) for i in 1:length(vaild) if vaild[i] == 1.0]
    new_weight = [weights[i] for i in 1:length(vaild) if vaild[i] == 1.0]
    return new_graph, new_weight
end

function calculate_energy(state, Δ, T_max, nodes)
    energy = 0.0
    for i in 1:length(state)
        for j in (i+1):length(state)
            if state[i] == state[j] && state[i] == 1
                energy += 2*π*862690 / Distance(nodes[i], nodes[j])^6
            end
        end
    end

    for i in 1:length(state)
        if state[i] == 1
            energy -= Δ[i](T_max)
        end
    end
    return energy
end

function generate_energy_plot(hamiltonian, T_max, T_step, howmany, subspace)
    clocks = 0.0:T_step:T_max
    val = []
    for t in clocks
        h = hamiltonian |> attime(t)
        mat_h = mat(h)
        eigvals, eigvecs, info = eigsolve(mat_h, howmany, :SR)
        append!(val, [eigvals])
    end
    return val, clocks
end

new_graph_nodes, new_graph_weights = generate_some_graph()

atoms = AtomList(new_graph_nodes)

radius = 4.8

detuning_max = 400
detuning_min = - detuning_max

T_max = 10.0
empty_time = 2.0
Δ = map(1:length(new_graph_nodes)) do idx
    piecewise_linear(clocks = [0.0, empty_time, T_max - empty_time, T_max], values = [detuning_min * new_graph_weights[idx], detuning_min * new_graph_weights[idx], detuning_max * new_graph_weights[idx], detuning_max * new_graph_weights[idx]])
end

Ω_max = 2π * 4
Ω = piecewise_linear(clocks = [0.0, empty_time, T_max - empty_time, T_max], values = [0, Ω_max, Ω_max, 0])

hamiltonian = rydberg_h(atoms, Ω = Ω, Δ = Δ)
subspace = blockade_subspace(atoms, radius)

eigvals, times = generate_energy_plot(hamiltonian, T_max, 0.01, 3, subspace)

PlotTime = [t for t in times]
fig = Figure(resolution = (1500, 800))
ax = Axis(fig[1, 1], title = "Energy Line", xlabel = "Time", ylabel = "log(Energy minus ground state energy)")

lines!(ax, PlotTime, [log(val[2] - val[1]) for val in eigvals], color = :blue, linewidth = 2, label = "first excitation")
lines!(ax, PlotTime, [log(val[3] - val[1]) for val in eigvals], color = :green, linewidth = 2, label = "second excitation")
# lines!(ax, PlotTime, [log(val[4] - val[1]) for val in eigvals], color = :brown, linewidth = 2, label = "third excitation")
# lines!(ax, PlotTime, [log(val[5] - val[1]) for val in eigvals], color = :red, linewidth = 2, label = "fourth excitation")
axislegend(ax; position = :lb, labelsize = 15)


fig
# nsites = length(new_graph_nodes)
# clocks = 0.0:1e-2:T_max
# prob2 = KrylovEvolution(zero_state(subspace), clocks, hamltonian)
# emulate!(prob2)

# bitstring_hist(prob2.reg; nlargest = 20)

# best_bit_strings = most_probable(prob2.reg, 10)
# desire = [0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0]
# Bloqade.plot(atoms, blockade_radius = radius; colors = [iszero(b) ? "white" : "red" for b ∈ desire])


# calculate_energy(best_bit_strings[9], Δ, T_max, new_graph_nodes)
# calculate_energy(desire, Δ, T_max, new_graph_nodes)