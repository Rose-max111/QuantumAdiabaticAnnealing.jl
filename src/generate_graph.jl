using Bloqade
using UnitDiskMapping
using Graphs


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


new_graph_nodes, new_graph_weights, radius = generate_factoring_new_graph(1, 1)
atoms = AtomList(new_graph_nodes)

detuning_max = (2π* 862690 / radius^6)
detuning_min = - detuning_max

T_max = 0.6
Δ = map(1:length(new_graph_nodes)) do idx
    piecewise_linear(clocks = [0.0, 0.1, 0.5, T_max], values = [detuning_min * new_graph_weights[idx], detuning_min * new_graph_weights[idx], detuning_max * new_graph_weights[idx], detuning_max * new_graph_weights[idx]])
end

Ω_max = 2π * 4
Ω = piecewise_linear(clocks = [0.0, 0.1, 0.5, T_max], values = [0, Ω_max, Ω_max, 0])

hamltonian = rydberg_h(atoms, Ω = Ω, Δ = Δ)
prob = SchrodingerProblem(zero_state(nqubits(hamltonian)), T_max, hamltonian)