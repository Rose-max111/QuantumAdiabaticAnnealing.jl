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


function generate_factoring_new_graph(x1::Int, x2::Int)
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
            if distance(mres.grid_graph.nodes[i].loc, mres.grid_graph.nodes[pin_one].loc) < radius
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


function generate_random_lattice(lattice_length::Int, lattice_constant::Float64, percentage::Float64)
    num_of_nodes = Int(floor(lattice_length * lattice_length * percentage))
    num_of_total_lattice = Int(lattice_length * lattice_length)
    list_of_nodes = randperm(num_of_total_lattice)[1:num_of_nodes] 
    # num_of_nodes = 20
    # list_of_nodes = [i for i in 1:20]

    nodes = Vector{Tuple{Float64, Float64}}()
    weights = ones(num_of_nodes)
    for i in 1:lattice_length
        for j in 1:lattice_length
            if ((i-1) * lattice_length + j) in list_of_nodes
                push!(nodes, ((i-1) * lattice_constant, (j-1)*lattice_constant))
            end
        end
    end
    return nodes, weights
end