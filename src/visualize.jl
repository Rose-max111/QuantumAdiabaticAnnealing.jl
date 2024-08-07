function show_transversal_graph_weight(n, m; gradient=nothing)
    locations, weights, input_id, output_id, input_layer_id = rule110_transverse_generate(n, m; gradient=gradient)
    G = SimpleGraph(length(weights))
    radius_square = 4.2
    for i in 1:length(locations)
        for j in i+1:length(locations)
            if (locations[i][1] - locations[j][1])^2 + (locations[i][2] - locations[j][2])^2 <= radius_square
                add_edge!(G, i, j)
            end
        end
    end

    show_graph(G, locs=locations; format=:svg, texts=map(string, weights), vertex_text_colors=fill("Black", nv(G)))
end

function rule110_gadget_plot()
    config = GraphDisplayConfig(vertex_size=9.0, vertex_line_width=1.3)
    graph, weights, locations, P, Q, R, Target = rule110_generate()
    locations = map(t -> (-30 * t[1], -30 * t[2]), locations)
    problem = GenericTensorNetwork(IndependentSet(graph, weights))
    max8_configs = solve(problem, SingleConfigMax(8))[]
    max8_configs = [max8_configs.orders[i].c.data for i in 1:8]
    sort!(max8_configs, lt=(a, b) -> (a[P] * 4 + a[Q] * 2 + a[R] * 1 < b[P] * 4 + b[Q] * 2 + b[R] * 1))
    show_configs(graph, locations, [max8_configs[j+4*i] for i = 0:1, j = 1:4]; format=:png, filename="gadget110.png", config=config, padding_left=20, padding_top=30, padding_bottom=30, padding_right=20)
end

function show_transversal_graph(n, m)
    locations, weights, input_id, output_id, input_layer_id = rule110_transverse_generate(n, m)
    G = SimpleGraph(length(weights))
    radius_square = 4.2
    for i in 1:length(locations)
        for j in i+1:length(locations)
            if (locations[i][1] - locations[j][1])^2 + (locations[i][2] - locations[j][2])^2 <= radius_square
                add_edge!(G, i, j)
            end
        end
    end

    weights = Int.(weights)
    colorspace = ["Blue", "Red", "Green", "Black"]

    locations = map(t -> (Int(50 * t[1]), Int(50 * t[2])), locations)
    # show_graph(G, locations; format = :svg,vertex_sizes = fill(0.2, nv(G)), vertex_colors = colorspace[weights], texts = map(string, 1:nv(G)), vertex_text_colors = fill("White", nv(G)))
    show_graph(G, locations; format=:svg, filename="rule110_transvarient.svg", vertex_sizes=fill(8, nv(G)), vertex_colors=colorspace[weights], texts=map(string, 1:nv(G)), vertex_text_colors=fill("White", nv(G)))
end

