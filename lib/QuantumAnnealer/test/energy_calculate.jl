using Test
using QuantumAnnealer
using KrylovKit
using BloqadeLattices: AtomList
using BloqadeExpr: rydberg_h
using KrylovKit: eigsolve

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


@testset "get_low_energy_state" begin
    new_graph_nodes, new_graph_weights = generate_some_graph()

    Δ = fill(30 * 2π, length(new_graph_nodes))
    Ω = fill(4 * 2π, length(new_graph_nodes))

    dmrg_energy,psi= get_low_energy_state(Δ, Ω, new_graph_nodes;outputlevel = 1)
    @show dmrg_energy

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

# @testset "cpu vs gpu" begin
#     new_graph_nodes, new_graph_weights = generate_random_lattice(5, 1.6, 0.8)

#     Δ = fill(30 * 2π, length(new_graph_nodes))
#     Ω = fill(4 * 2π, length(new_graph_nodes))

#     @info "Gpu processing"
#     dmrg_energy_gpu, psi_gpu = get_low_energy_state_gpu(Δ, Ω, new_graph_nodes;outputlevel = 1)
#     @info "Cpu processing"
#     dmrg_energy_cpu, psi_cpu = get_low_energy_state(Δ, Ω, new_graph_nodes;outputlevel = 1)
#     @show dmrg_energy_cpu
#     @show dmrg_energy_gpu

#     @test abs(dmrg_energy_cpu[1] - dmrg_energy_gpu[1]) / abs(dmrg_energy_cpu[1]) < 1e-2
#     @test abs(dmrg_energy_cpu[2] - dmrg_energy_gpu[2]) / abs(dmrg_energy_cpu[2]) < 1e-2
#     @test abs(dmrg_energy_cpu[3] - dmrg_energy_gpu[3]) / abs(dmrg_energy_cpu[3]) < 1e-2
# end