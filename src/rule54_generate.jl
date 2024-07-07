function rule54_generate()
    el = Edge.([(1,2),
    (2,3),
    (3,4),
    (4,5),
    (5,6),
    (6,7),
    (6,10), 
    (7,8),
    (7,9),
    (7,10),
    (7,11),
    (8,9),
    (8,10),
    (8,11),
    (8,12),
    (10,11),
    (10,13),
    (11,12),
    (11,13),
    (12,15),
    (13,14),
    (14,18),
    (15,16),
    (15,17),
    (16,17),
    (16,21),
    (16,20),
    (16,19),
    (17,18),
    (17,20),
    (17,21),
    (18,21),
    (18,23),
    (19,20),
    (19,55),
    (20,21),
    (20,22),
    (21,22),
    (22,32),
    (22,29),
    (23,24),
    (24,25),
    (24,27),
    (25,26),
    (27,28),
    (28,36),
    (29,30),
    (30,31),
    (31,32),
    (31,33),
    (32,54),
    (33,34),
    (33,35),
    (34,35),
    (34,39),
    (34,38),
    (34,37),
    (35,36),
    (35,37),
    (35,38),
    (36,37),
    (37,38),
    (37,40),
    (38,39),
    (38,40),
    (39,42),
    (40,41),
    (41,48),
    (42,43),
    (42,44),
    (43,48),
    (43,46),
    (43,45),
    (43,44),
    (44,45),
    (44,46),
    (44,53),
    (45,53),
    (45,47),
    (45,46),
    (46,47),
    (46,48),
    (48,49),
    (49,50),
    (50,51),
    (51,52),
    (54,56),
    (54,57),
    (55,56),
    (55,58),
    (56,57),
    (56,58),
    (57,58),
    (57,59),
    (58,59),
    (59,60)])

    locations = [
    (-8,-6),
    (-8,-4),
    (-8,-2),
    (-6,-2),
    (-4,-2),
    (-2,-1),
    (-2,0),
    (-2,1),
    (-3,2),
    (0,0),
    (0,1),
    (0,2),
    (1,-1),
    (2,0),
    (2,3),
    (3,4),
    (3,3),
    (4,1),
    (4,6),
    (5,4),
    (5,3),
    (6,3),
    (5,-1),
    (6,-2),
    (6,-4),
    (6,-6),
    (8,-2),
    (10,-2),
    (7,1),
    (8,1),
    (9,3),
    (8,4),
    (11,2),
    (12,1),
    (12,0),
    (12,-1),
    (14,0),
    (14,1),
    (14,2),
    (15,-1),
    (16,0),
    (16,3),
    (17,3),
    (17,4),
    (19,4),
    (19,3),
    (20,3),
    (18,1),
    (19,-1),
    (20,-2),
    (20,-4),
    (20,-6),
    (18,6),
    (8,6),
    (5,7),
    (7,7),
    (8,8),
    (6,8),
    (7,10),
    (7,11)]

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
    2,
    2,
    4,
    4,
    3,
    2,
    4,
    4,
    2,
    2,
    3,
    2,
    1,
    2,
    2,
    1,
    1,
    2,
    2,
    2,
    4,
    4,
    2,
    4,
    4,
    2,
    2,
    2,
    2,
    4,
    4,
    4,
    4,
    1,
    3,
    2,
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
    1]

    return SimpleGraph(el), locations, weights, 1, 26, 52
end
# weights[1] = -50
# weights[26] = -50
# weights[52] = -50

# locations = map(t->t.*100, locations)
# locations = map(t->(t[1], -t[2]), locations)
# graph = SimpleGraph(el)

# LuxorGraphPlot.show_graph(graph, locations;format = :png)

# problem = GenericTensorNetwork(IndependentSet(graph, weights));
# max3_configs = read_config(solve(problem, SingleConfigMax(3))[])

# spectrum = solve(problem, SizeMax(3))[]

# max_config_weighted = solve(problem, SingleConfigMax())[]
# show_graph(graph, locations; format=:png, vertex_colors=
        #   [iszero(max_config_weighted.c.data[i]) ? "white" : "red" for i=1:nv(graph)])

# show_configs(graph, locations, [max3_configs[i] for i=1:3, j=1:1]; padding_top=20)