using Test, QuantumAdiabaticAnnealing, CUDA
using QuantumAdiabaticAnnealing: linear_to_cartesian, cartesian_to_linear, hasparent, child_nodes, SimulatedAnnealingHamiltonian,
    parent_logic, random_state, natom, track_equilibration!, step!, HeatBath, Metropolis, step_kernel!

@testset "SimulatedAnnealingHamiltonian" begin
    sa = SimulatedAnnealingHamiltonian(3, 2)
    @test sa isa SimulatedAnnealingHamiltonian
    @test natom(sa) == 6
    @test sa.n == 3
    @test sa.m == 2
    @test linear_to_cartesian(sa, 4) == (1, 2)
    @test cartesian_to_linear(sa, 1, 2) == 4
    @test hasparent(sa, 1) == false
    @test hasparent(sa, 4) == true
    @test random_state(sa, 30) isa Matrix
    @test parent_logic(sa, 4) == (3, 1, 2, 4)
    @test child_nodes(sa, 1) == (6, 4, 5)
end

@testset "sa (cpu)" begin
    sa = SimulatedAnnealingHamiltonian(4, 2)
    nbatch = 30
    @testset "Metroplis" begin
        state = random_state(sa, nbatch)
        step!(Metropolis(), sa, state, 1.0, 1)
        @test track_equilibration!(Metropolis(), sa, state) isa SimulatedAnnealingHamiltonian
    end
    @testset "HeatBath" begin
        state = random_state(sa, nbatch)
        step!(HeatBath(), sa, state, 1.0, 1)
        @test track_equilibration!(HeatBath(), sa, state) isa SimulatedAnnealingHamiltonian
    end
end

if CUDA.functional()
    @testset "sa (gpu)" begin
        sa = SimulatedAnnealingHamiltonian(4, 2)
        nbatch = 30
        @testset "Metropolis" begin
            state = CuArray(random_state(sa, nbatch))
            @test step!(Metropolis(), sa, state, 1.0, 1) isa CuMatrix
            @test track_equilibration!(Metropolis(), sa, state) isa SimulatedAnnealingHamiltonian
        end
        @testset "HeatBath" begin
            state = CuArray(random_state(sa, nbatch))
            @test step!(HeatBath(), sa, state, 1.0, 1) isa CuMatrix
            @test track_equilibration!(HeatBath(), sa, state) isa SimulatedAnnealingHamiltonian
        end
    end
end


using KrylovKit
using QuantumAdiabaticAnnealing:toy_model_transition_matrix, toy_model_state_energy
@testset "sa_energy(gpu)" begin
    n_size = 6
    
    # initiate SA condition
    sa = SimulatedAnnealingHamiltonian(n_size, 2)
    nbatch = 5000
    state = CuArray(random_state(sa, nbatch))
    anneal_time = 10000
    end_temp = 2.0
    energy_gradient = CUDA.ones(anneal_time)
    
    @btime track_equilibration!(HeatBath(), $sa, $state, energy_gradient, fill(end_temp, anneal_time))
    
    cpu_state = Array(state)

    # calculate the corresponding sample-distribution energy
    state_energy = [calculate_energy(sa, cpu_state, energy_gradient, i) for i in 1:nbatch]
    state_energy_average = 1.0 * sum(state_energy) / nbatch
    @info "sample energy is $state_energy_average"

    # # calculate corresponding sample-distribution 
    # state_10 = [sum(cpu_state[:,i] .* [2^(j-1) for j in 1:sa.n*sa.m]) for i in 1:size(cpu_state)[2]]
    # figure_count = [count(==(i), state_10) for i in 0:2^(sa.n*sa.m)-1]
    # distribution_count = 1.0 .* figure_count ./ sum(figure_count)

    # use ED to calculate the exact boltzmann distribution
    n = n_size
    Temp = end_temp
    P = toy_model_transition_matrix(HeatBath(), n, Temp; period_condition = true)
    eigvals, eigvecs, infos = eigsolve(P, rand(Float64, size(P)[1]), 2, :LR; maxiter = 5000)
    
    # calculate the exact boltzmann distribution corresponding energy
    ED_state_energy = [toy_model_state_energy(i, n; period_condition = true) for i in 0:2^(2n)-1]
    maxeigvec_energy = sum(ED_state_energy .* (abs.(Real.(eigvecs[1])))) / abs(sum(Real.(eigvecs[1])))

    # # calculate exact boltzmann distribution
    # distribution = [abs(Real.(eigvecs[1][i])) for i in 1:size(eigvecs[1])[1]] ./ abs(sum(Real.(eigvecs[1])))

    @info "exact boltzmann energy is $maxeigvec_energy"
    @test abs(maxeigvec_energy - state_energy_average) / maxeigvec_energy <= 0.01
end