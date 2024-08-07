using QuantumAdiabaticAnnealing, Test
using QuantumAdiabaticAnnealing.QuantumAnnealing: get_low_energy_state, get_low_energy_state_gpu
using QuantumAdiabaticAnnealing: generate_random_lattice
using Bloqade, KrylovKit
using Bloqade: AtomList, rydberg_h
using KrylovKit: eigsolve
using CUDA; CUDA.allowscalar(false)

@testset "distance" begin
    @test distance((0,0), (3,4)) ≈ 5.0
    @test distance((-1,-2), (3,1)) ≈ 5.0
end

@testset "get_low_energy_state" begin
    new_graph_nodes, new_graph_weights = generate_some_graph()

    Δ = fill(30 * 2π, length(new_graph_nodes))
    Ω = fill(4 * 2π, length(new_graph_nodes))

    dmrg_energy,psi= get_low_energy_state_gpu(Δ, Ω, new_graph_nodes;outputlevel = 1)
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