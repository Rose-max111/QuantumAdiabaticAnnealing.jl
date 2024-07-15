function rule110_generate()
    el = Edge.([(1,5),
        (1,6),
        (2,5),
        (2,7),
        (2,9),
        (3,9),
        (3,11),
        (4,10),
        (5,6),
        (5,7),
        (5,9),
        (6,7),
        (6,8),
        (7,8),
        (7,9),
        (7,11),
        (8,10),
        (9,10),
        (9,11),
        (12,2)
        ])
    weights = [1, 2, 1, 1, 2, 2, 2, 2, 2, 2, 1, 1]

    locations = [(1,3),
            (3,5),
            (6,4),
            (4,0),
            (2,4),
            (2,3),
            (3,4),
            (3,2),
            (4,4),
            (4,2),
            (5,4),
            (3,7)]
    return SimpleGraph(el), weights, locations, 3, 12, 1, 4
end

function transversal_graph()
    locations = [(1,0),
                (3,0),
                (5,0),
                (6,0),
                (9,0),
                (7.5,0.5),
                (8.5,1.5),
                (6.5,1.5),
                (7.5,2.5),
                (9,2.5),
                (6.5,3),
                (3,2),
                (0,3),
                (1,3),
                (2,3),
                (3,3),
                (4,3),
                (5,4),
                (4,4),
                (3,5),
                (2,5),
                (2,7),
                (2,8)]
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
    P = 2
    O = 23
    # L input = 13, R input = 5
    return locations, weights, P, O
end

function rule110_transverse_generate(ninput::Int, Time::Int; gradient = nothing)
    locations = []
    weights = []
    basic_loc, basic_weight, input, output = transversal_graph()
    
    input_id = []
    output_id = []
    input_layer_id = Vector{Vector{Int}}()

    if gradient == nothing
        gradient = ones(Int, Time)
    end
    padding_position = [0,0]
    input_additional_weight = 0
    last_graph_size = 0
    for row in 1:Time
        inputL = 0
        inputR = 0
        for column in 1:ninput
            this_graph_loc = map(x -> (x[1] + padding_position[1], x[2] + padding_position[2]), basic_loc)
            this_graph_weight = copy(basic_weight)
            if row == 1
                push!(input_id, length(locations) + 2)
            end
            if column == 1
                inputL = length(locations) + 13
            end
            if column == ninput
                inputR = length(locations) + 5
            end

            if gradient != nothing
                this_graph_weight .*= gradient[row]
            end
            # deal with the border between this row and last row (copy gadget)
            this_graph_weight[2] += input_additional_weight * gradient[row]

            # deal with the border between this column and last column (copy gadget)
            if column != 1
                this_graph_weight[1] += gradient[row]
                this_graph_weight[13] += gradient[row]
                weights[end - last_graph_size + 5] += gradient[row]
                weights[end - last_graph_size + 10] += gradient[row]
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
                this_graph_weight[end] += gradient[row + 1]
                append!(locations, this_graph_loc[1:end])
                append!(weights, this_graph_weight[1:end])
                last_graph_size = length(this_graph_loc)
            else
                this_graph_weight[end - 1] -= gradient[row]
                push!(output_id, length(locations) + length(this_graph_loc) - 1)
                append!(locations, this_graph_loc[1:end - 1])
                append!(weights, this_graph_weight[1:end - 1])
                last_graph_size = length(this_graph_loc) - 1
            end
            padding_position[1] += 10
        end
        push!(input_layer_id, [inputL, inputR])
        padding_position[1] = -row
        padding_position[2] += 10
        input_additional_weight = 1
    end
    return locations, weights, input_id, output_id, input_layer_id
end


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

    show_graph(G, locs = locations; format = :svg, vertex_colors = colorspace[weights], texts = map(string, 1:nv(G)), vertex_text_colors = fill("White", nv(G)))
end

function show_transversal_graph_weight(n, m; gradient = nothing)
    locations, weights, input_id, output_id, input_layer_id = rule110_transverse_generate(n,m; gradient = gradient)
    G = SimpleGraph(length(weights))
    radius_square = 4.2
    for i in 1:length(locations)
        for j in i+1:length(locations)
            if (locations[i][1] - locations[j][1])^2 + (locations[i][2]-locations[j][2])^2 <= radius_square
                add_edge!(G, i, j)
            end
        end
    end

    show_graph(G, locs = locations; format = :svg, texts = map(string, weights), vertex_text_colors = fill("Black", nv(G)))
end