using QuantumAdiabaticAnnealing, Test
using QuantumAdiabaticAnnealing: get_low_energy_state, generate_some_graph, distance
using Bloqade, KrylovKit
using Bloqade: AtomList, rydberg_h
using KrylovKit: eigsolve

@testset "distance" begin
    @test distance((0,0), (3,4)) ≈ 5.0
    @test distance((-1,-2), (3,1)) ≈ 5.0
end

@testset "get_low_energy_state" begin
    new_graph_nodes, new_graph_weights = generate_some_graph()
    # Δ = rand(50:200, length(new_graph_nodes))
    # Ω = rand(1:10, length(new_graph_nodes))

    Δ = [224.02199999999996,448.0439999999999,224.02199999999996,448.0439999999999,448.0439999999999,448.0439999999999,448.0439999999999,224.02199999999996,224.02199999999996,448.0439999999999,448.0439999999999,448.0439999999999,448.0439999999999,224.02199999999996,224.02199999999996,448.0439999999999,224.02199999999996]
    Ω = ones(length(new_graph_nodes))

    dmrg_energy,psi= get_low_energy_state(Δ, Ω, new_graph_nodes)

    H0 = rydberg_h(AtomList(new_graph_nodes), Ω = Ω, Δ = Δ)
    @show H0
    mat_H0 = mat(H0)
    eigvals, eigvecs, infos = eigsolve(mat_H0, 3, :SR)

    @show dmrg_energy[1], eigvals[1]
    @show dmrg_energy[2], eigvals[2]
    @show dmrg_energy[3], eigvals[3]
    @test abs(dmrg_energy[1] - eigvals[1]) / abs(dmrg_energy[1]) < 1e-4
    @test abs(dmrg_energy[2] - eigvals[2]) / abs(dmrg_energy[2]) < 1e-3
    @test abs(dmrg_energy[3] - eigvals[3]) / abs(dmrg_energy[3]) < 1e-2
end   