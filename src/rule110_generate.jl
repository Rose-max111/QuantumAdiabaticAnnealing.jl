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
    weights = [1,
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
            1]
    P = 2
    O = 23
    return locations, weights, P, O
end

function rule110_transverse_generate(ninput::Int, Time::Int)
    locations = []
    weights = []
    basic_loc, basic_weight, input, output = transversal_graph()
    
    input_id = []
    output_id = []

    padding_position = [0,0]
    input_additional_weight = 0
    last_graph_size = 0
    for row in 1:Time
        for column in 1:ninput
            this_graph_loc = map(x -> (x[1] + padding_position[1], x[2] + padding_position[2]), basic_loc)
            this_graph_weight = copy(basic_weight)
            if row == 1
                push!(input_id, length(locations) + 2)
            end
            # deal with the border between this row and last row
            this_graph_weight[2] += input_additional_weight
            # deal with the border between this column and last column
            if column != 1
                this_graph_weight[1] += 1
                this_graph_weight[13] += 1
                weights[end - last_graph_size + 5] += 1
                weights[end - last_graph_size + 10] += 1
            end
            # deal with the border between this row and next row
            if row != Time
                this_graph_weight[end] += 1
                append!(locations, this_graph_loc[1:end])
                append!(weights, this_graph_weight[1:end])
                last_graph_size = length(this_graph_loc)
            else
                this_graph_weight[end - 1] -= 1
                push!(output_id, length(locations) + length(this_graph_loc) - 1)
                append!(locations, this_graph_loc[1:end - 1])
                append!(weights, this_graph_weight[1:end - 1])
                last_graph_size = length(this_graph_loc) - 1
            end
            padding_position[1] += 10
        end
        padding_position[1] = -row
        padding_position[2] += 10
        input_additional_weight = 1
    end
    return locations, weights, input_id, output_id
end
