using Random
using CUDA
using BitBasis
using QuantumAdiabaticAnnealing: rule110

struct ParallelSimulatedAnnealingState{MT<:AbstractMatrix{Bool}}
    n::Int # number of atoms per layer
    m::Int # number of full layers (count the last layer!)
    state::MT  # state of the system, rows are atoms, columns are batches
end
function ParallelSimulatedAnnealingState(n::Int, m::Int, nbatch::Int)
    num_atoms = n * m
    state = rand(Bool, num_atoms, nbatch)
    return ParallelSimulatedAnnealingState(n, m, state)
end
nbatch(sa::ParallelSimulatedAnnealingState) = size(sa.state, 2)
natom(sa::ParallelSimulatedAnnealingState) = size(sa.state, 1)
atoms(sa::ParallelSimulatedAnnealingState) = Base.OneTo(natom(sa))
Base.copy(sa::ParallelSimulatedAnnealingState) = ParallelSimulatedAnnealingState(sa.n, sa.m, copy(sa.state))
hasparent(sa::ParallelSimulatedAnnealingState, node::Int) = node > sa.n
linear_to_cartesian(sa::ParallelSimulatedAnnealingState, node::Int) = CartesianIndices((sa.n, sa.m))[node]
cartesian_to_linear(sa::ParallelSimulatedAnnealingState, ci::CartesianIndex) = LinearIndices((sa.n, sa.m))[ci]
CUDA.cu(sa::ParallelSimulatedAnnealingState) = ParallelSimulatedAnnealingState(sa.n, sa.m, CUDA.CuArray(sa.state))

# evaluate the energy of the i-th gadget (involving atoms i and its parents)
function evaluate_parent(sa::ParallelSimulatedAnnealingState, inode::Int, ibatch::Int)
    idp = parent_logic(sa, inode)
    trueoutput = rule110(sa.state[idp[1], ibatch], sa.state[idp[2], ibatch], sa.state[idp[3], ibatch])
    return trueoutput ⊻ sa.state[idp[4], ibatch]
end
function calculate_energy(sa::ParallelSimulatedAnnealingState, ibatch::Int)
    return sum(i->evaluate_parent(sa, i, ibatch), atoms(sa))
end
function parent_logic(sa::ParallelSimulatedAnnealingState, node::Int)
    n = sa.n
    ci = linear_to_cartesian(sa, node)
    (
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]-1, n), ci[2]-1)),  # periodic boundary condition
        cartesian_to_linear(sa, CartesianIndex(ci[1], ci[2]-1)),
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+1, n), ci[2]-1)),
        node
    )
end

function child_nodes(sa::ParallelSimulatedAnnealingState, node::Int)
    n = sa.n
    ci = linear_to_cartesian(sa, node)
    (
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]-1, n), ci[2]+1)),  # periodic boundary condition
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1], n), ci[2]+1)),
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+1, n), ci[2]+1)),
    )
end

# i can be -1, 0, 1
function child_logic(sa::ParallelSimulatedAnnealingState, node::Int, i::Int)
    n = sa.n
    ci = linear_to_cartesian(sa, node)
    (
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+i-1, n), ci[2])),  # periodic boundary condition
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+i, n), ci[2])),
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+i+1, n), ci[2])),
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+i, n), ci[2]+1)),
    )
end

using Test
@testset "SimulatedAnnealingHamiltonian" begin
    sa = ParallelSimulatedAnnealingState(3, 2, 30)
    @test sa isa ParallelSimulatedAnnealingState
    @test natom(sa) == 6
    @test nbatch(sa) == 30
    @test copy(sa) isa ParallelSimulatedAnnealingState
    @test sa.n == 3
    @test sa.m == 2
    @test sa.state != 0
    @test linear_to_cartesian(sa, 4) == CartesianIndex(1, 2)
    @test cartesian_to_linear(sa, CartesianIndex(1, 2)) == 4
    @test hasparent(sa, 1) == false
    @test hasparent(sa, 4) == true
    @test parent_logic(sa, 4) == (3, 1, 2, 4)
    @test child_nodes(sa, 1) == (6, 4, 5)
    @test child_logic(sa, 1, -1) == (2, 3, 1, 6)
    @test child_logic(sa, 1, 0) == (3, 1, 2, 4)
    @test child_logic(sa, 1, 1) == (1, 2, 3, 5)
end

abstract type TransitionRule end
struct HeatBath <: TransitionRule end
struct Metroplis <: TransitionRule end

function step!(rule::TransitionRule, sa::ParallelSimulatedAnnealingState, Temp::Float64, node::Int)
    for ibatch in 1:nbatch(sa)
        step_kernel!(rule, sa, Temp, node, ibatch)
    end
end
function step!(rule::TransitionRule, sa::ParallelSimulatedAnnealingState{<:CuMatrix}, Temp::Float64, node::Int)
    function kernel(rule, sa, Temp, node)
        ibatch = threadIdx().x
        step_kernel!(rule, sa, Temp, node, ibatch)
        return nothing
    end
    @cuda threads = nbatch(sa) kernel(rule, sa, Temp, node)
end

@inline function step_kernel!(rule::TransitionRule, sa::ParallelSimulatedAnnealingState, Temp::Float64, node::Int, ibatch::Int)
    ΔE = 0
    n, m = sa.n, sa.m
    i, j = linear_to_cartesian(sa, node).I
    if j > 1 # not the first layer
        ΔE -= 2 * evaluate_parent(sa, node, ibatch)
    end
    if j < m # not the last layer
        cnodes = child_nodes(sa, node)
        for node in cnodes
            ΔE -= evaluate_parent(sa, node, ibatch)
        end
        # flip the node
        sa.state[node, ibatch] ⊻= true
        for node in cnodes
            ΔE += evaluate_parent(sa, node, ibatch)
        end
    end
    if rand() < prob_accept(rule, Temp, ΔE)
        ΔE
    else  # revert the flip
        sa.state[node, ibatch] ⊻= true
        0
    end
end
prob_accept(::Metroplis, Temp, ΔE::T) where T<:Real = ΔE < 0 ? one(T) : exp(- (ΔE) / Temp)
prob_accept(::HeatBath, Temp, ΔE::Real) = inv(1 + exp(ΔE / Temp))

function track_equilibration!(rule::TransitionRule, sa::ParallelSimulatedAnnealingState, tempscale = 4 .- (1:100 .-1) * 0.04)
    best_state = copy(sa)
    best_energy = calculate_energy.(Ref(sa), 1:nbatch(sa))
    # NOTE: do we really need niters? or just set it to 1?
    for Temp in tempscale
        # NOTE: do we really need to shuffle the nodes?
        for node in 1:natom(sa)
            eng = step!(rule, sa, Temp, node)
            mask = eng .< best_energy
            best_energy[mask] .= eng[mask]
            best_state.state[:, mask] .= sa.state[:, mask]
        end
    end
    return sa
end

@testset "sa (cpu)" begin
    sa = ParallelSimulatedAnnealingState(4, 2, 30)
    @testset "Metroplis" begin
        step!(Metroplis(), sa, 1.0, 1)
        @test track_equilibration!(Metroplis(), sa) isa ParallelSimulatedAnnealingState
    end
    @testset "HeatBath" begin
        step!(HeatBath(), sa, 1.0, 1)
        @test track_equilibration!(HeatBath(), sa) isa ParallelSimulatedAnnealingState
    end
end