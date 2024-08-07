using QuantumAdiabaticAnnealing, Test
using QuantumAdiabaticAnnealing:generate_random_lattice

@testset "generate_random_lattice" begin
    nodes, weights = generate_random_lattice(5, 1.5, 0.8)
    for i in 1:length(nodes)
        if nodes[i][1] < 0 || nodes[i][2] < 0 || nodes[i][1] > 6 || nodes[i][2] > 6
            @test false
        end
        for j in i+1:length(nodes)
            if nodes[i][1] == nodes[j][1] && nodes[i][2] == nodes[j][2]
                @test false
            end
        end
    end
    @test true
end