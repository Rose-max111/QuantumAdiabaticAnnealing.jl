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