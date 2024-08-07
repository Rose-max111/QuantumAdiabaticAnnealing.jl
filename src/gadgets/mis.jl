struct UnitDiskMISGadget
    weights::Vector{Float64}
    inputs::Vector{Int}
    outputs::Vector{Int}
    # for unit disk embedding
    coordinates::Vector{Tuple{Float64, Float64}}
    radius::Float64
end
topology(g::UnitDiskMISGadget) = unit_disk_graph(g.coordinates, g.radius)

function mis_gadget(::CellularAutomata1D{54})
    locations = [(-9, -4),
        (-9, -2),
        (-7, -2),
        (-5, -2),
        (-3, -2),
        (-1, -1),
        (-2, 1),
        (-2, 2),
        (-4, 2),
        (0, 1),
        (0, 2),
        (-1, 4),
        (1, 1),
        (3, 0),
        (1, 5),
        (3, 5),
        (3, 4),
        (4, 2),
        (4, 7),
        (5, 5),
        (5, 4),
        (6, 4),
        (3, -2),
        (5, -2),
        (7, -2),
        (7, -4),
        (9, -2),
        (11, -2),
        (13, -2),
        (7, 2),
        (9, 2),
        (10, 4),
        (8, 5),
        (12, 4),
        (12, 2),
        (14, 2),
        (14, 1),
        (16, 2),
        (16, 1),
        (15, 4),
        (17, 1),
        (19, 0),
        (17, 5),
        (19, 4),
        (19, 5),
        (21, 5),
        (21, 4),
        (22, 4),
        (20, 2),
        (19, -2),
        (21, -2),
        (23, -2),
        (23, -4),
        (15, -1),
        (20, 7),
        (8, 7),
        (5, 8),
        (7, 8),
        (6, 9),
        (8, 9),
        (7, 11)]


    weights = [1,
        2,
        2,
        2,
        2,
        2,
        4,
        4,
        1,
        4,
        4,
        2,
        2,
        3,
        2,
        4,
        4,
        2,
        2,
        4,
        4,
        2,
        2,
        2,
        3,
        1,
        2,
        2,
        2,
        1,
        1,
        2,
        2,
        2,
        2,
        4,
        4,
        4,
        4,
        2,
        2,
        3,
        2,
        4,
        4,
        4,
        4,
        1,
        2,
        2,
        2,
        2,
        1,
        2,
        1,
        2,
        2,
        2,
        2,
        2,
        1]
    return UnitDiskMISGadget(weights, [1, 26, 53], [61], locations, sqrt(5.0))
end

function mis_gadget(::CellularAutomata1D{110})
    weights = [1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1]
    locations = [(1, 3),
        (3, 5),
        (6, 4),
        (4, 0),
        (2, 4),
        (2, 3),
        (3, 4),
        (3, 2),
        (4, 4),
        (4, 2),
        (5, 4),
        (3, 7)]
    return UnitDiskMISGadget(weights, [3, 12, 1], [4], locations, sqrt(4.2))
end

function transversal_mis_gadget(::CellularAutomata1D{110})
    locations = [(1, 0),
        (3, 0),
        (5, 0),
        (6, 0),
        (9, 0),
        (7.5, 0.5),
        (8.5, 1.5),
        (6.5, 1.5),
        (7.5, 2.5),
        (9, 2.5),
        (6.5, 3),
        (3, 2),
        (0, 3),
        (1, 3),
        (2, 3),
        (3, 3),
        (4, 3),
        (5, 4),
        (4, 4),
        (3, 5),
        (2, 5),
        (2, 7),
        (2, 8)]
    weights = Float64.([1,
        3,
        2,
        2,
        1,
        4,
        4,
        4,
        4,
        1,
        2,
        2,
        1,
        1,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        1])
    # L input = 13, R input = 5
    # return locations, weights, P, O
    return UnitDiskMISGadget(weights, [13, 2, 5], [23], locations, sqrt(4.2))
end

function transversal_mis_gadget(::CellularAutomata1D{110}, ninput::Int, Time::Int; gradient=nothing)
    locations = Tuple{Float64, Float64}[]
    weights = Float64[]
    g = transversal_mis_gadget(CellularAutomata1D{110}())
    basic_loc, basic_weight = g.coordinates, g.weights

    input_id = Int[]
    inputL_id = Int[]
    inputR_id = Int[]
    output_id = Int[]
    input_layer_id = Vector{Vector{Int}}()

    if gradient === nothing
        gradient = ones(Int, Time)
    end
    padding_position = [0, 0]
    input_additional_weight = 0
    last_graph_size = 0
    for row in 1:Time
        for column in 1:ninput
            this_graph_loc = map(x -> (x[1] + padding_position[1], x[2] + padding_position[2]), basic_loc)
            this_graph_weight = copy(basic_weight)
            # top, left and right inputs
            column == 1 && push!(inputL_id,  length(locations) + g.inputs[1])
            row == 1 && push!(input_id, length(locations) + g.inputs[2])
            column == ninput && push!(inputR_id, length(locations) + g.inputs[3])

            if gradient !== nothing
                this_graph_weight .*= gradient[row]
            end
            # deal with the border between this row and last row (copy gadget)
            this_graph_weight[2] += input_additional_weight * gradient[row]

            # deal with the border between this column and last column (copy gadget)
            if column != 1
                this_graph_weight[1] += gradient[row]
                this_graph_weight[13] += gradient[row]
                weights[end-last_graph_size+5] += gradient[row]
                weights[end-last_graph_size+10] += gradient[row]
            end

            # deal with the last column specially (remove the rightest part)
            if column == ninput
                this_graph_weight[2] -= gradient[row]
                this_graph_weight[18] -= gradient[row]
                this_graph_weight = this_graph_weight[vcat(1:2, 12:end)]
                this_graph_loc = this_graph_loc[vcat(1:2, 12:end)]
                inputR = length(locations) + 9
                if column == 1
                    inputL = length(locations) + 4
                end
            end

            # deal with the border between this row and next row
            if row != Time
                this_graph_weight[end] += gradient[row+1]
                append!(locations, this_graph_loc[1:end])
                append!(weights, this_graph_weight[1:end])
                last_graph_size = length(this_graph_loc)
            else
                this_graph_weight[end-1] -= gradient[row]
                push!(output_id, length(locations) + length(this_graph_loc) - 1)
                append!(locations, this_graph_loc[1:end-1])
                append!(weights, this_graph_weight[1:end-1])
                last_graph_size = length(this_graph_loc) - 1
            end
            padding_position[1] += 10
        end
        push!(input_layer_id, [inputL, inputR])
        padding_position[1] = -row
        padding_position[2] += 10
        input_additional_weight = 1
    end
    #return locations, weights, input_id, output_id, input_layer_id
    return UnitDiskMISGadget(weights, input_id, output_id, locations, sqrt(4.2))
end