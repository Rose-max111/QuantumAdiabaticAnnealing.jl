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

    Δ = [392.00100000000003,784.0020000000001,392.00100000000003,784.0020000000001,784.0020000000001,784.0020000000001,784.0020000000001,392.00100000000003,392.00100000000003,784.0020000000001,784.0020000000001,784.0020000000001,784.0020000000001,392.00100000000003,392.00100000000003,784.0020000000001,392.00100000000003]
    Ω = ones(length(new_graph_nodes))

    dmrg_energy0, dmrg_energy1, dmrg_energy2 = get_low_energy_state(Δ, Ω, new_graph_nodes)

    H0 = rydberg_h(AtomList(new_graph_nodes), Ω = Ω, Δ = Δ)
    @show H0
    mat_H0 = mat(H0)
    eigvals, eigvecs, infos = eigsolve(mat_H0, 3, :SR)

    # println("dmrg_energy0: ", dmrg_energy0, "dmrg_energy1: ", dmrg_energy1, "dmrg_energy2: ", dmrg_energy2)
    # println("exact_0: ", eigvals[1], "exact_1: ", eigvals[2], "exact_2: ", eigvals[3])
    @show dmrg_energy0, eigvals[1]
    @show dmrg_energy1, eigvals[2]
    @show dmrg_energy2, eigvals[3]
    @test abs(dmrg_energy0 - eigvals[1]) / abs(dmrg_energy0) < 1e-4
    @test abs(dmrg_energy1 - eigvals[2]) / abs(dmrg_energy1) < 1e-4
    @test abs(dmrg_energy2 - eigvals[3]) / abs(dmrg_energy2) < 1e-4
end   