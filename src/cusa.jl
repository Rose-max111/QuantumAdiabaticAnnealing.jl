struct SimulatedAnnealingHamiltonian
    n::Int # number of atoms per layer
    m::Int # number of full layers (count the last layer!)
end
natom(sa::SimulatedAnnealingHamiltonian) = sa.n * sa.m
atoms(sa::SimulatedAnnealingHamiltonian) = Base.OneTo(natom(sa))
function random_state(sa::SimulatedAnnealingHamiltonian, nbatch::Int)
    return rand(Bool, natom(sa), nbatch)
end
hasparent(sa::SimulatedAnnealingHamiltonian, node::Int) = node > sa.n
linear_to_cartesian(sa::SimulatedAnnealingHamiltonian, node::Int) = CartesianIndices((sa.n, sa.m))[node]
cartesian_to_linear(sa::SimulatedAnnealingHamiltonian, ci::CartesianIndex) = LinearIndices((sa.n, sa.m))[ci]

# evaluate the energy of the i-th gadget (involving atoms i and its parents)
function evaluate_parent(sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, inode::Int, ibatch::Int)
    idp = parent_logic(sa, inode)
    trueoutput = rule110(state[idp[1], ibatch], state[idp[2], ibatch], state[idp[3], ibatch])
    return trueoutput ⊻ state[idp[4], ibatch]
end
function calculate_energy(sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, ibatch::Int)
    return sum(i->evaluate_parent(sa, state, i, ibatch), sa.n+1:natom(sa))
end
function parent_logic(sa::SimulatedAnnealingHamiltonian, node::Int)
    n = sa.n
    ci = linear_to_cartesian(sa, node)
    (
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]-1, n), ci[2]-1)),  # periodic boundary condition
        cartesian_to_linear(sa, CartesianIndex(ci[1], ci[2]-1)),
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+1, n), ci[2]-1)),
        node
    )
end

function child_nodes(sa::SimulatedAnnealingHamiltonian, node::Int)
    n = sa.n
    ci = linear_to_cartesian(sa, node)
    (
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]-1, n), ci[2]+1)),  # periodic boundary condition
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1], n), ci[2]+1)),
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+1, n), ci[2]+1)),
    )
end

# i can be -1, 0, 1
function child_logic(sa::SimulatedAnnealingHamiltonian, node::Int, i::Int)
    n = sa.n
    ci = linear_to_cartesian(sa, node)
    (
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+i-1, n), ci[2])),  # periodic boundary condition
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+i, n), ci[2])),
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+i+1, n), ci[2])),
        cartesian_to_linear(sa, CartesianIndex(mod1(ci[1]+i, n), ci[2]+1)),
    )
end

abstract type TransitionRule end
struct HeatBath <: TransitionRule end
struct Metropolis <: TransitionRule end

function step!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, Temp::Float64, node::Int)
    for ibatch in 1:size(state, 2)
        step_kernel!(rule, sa, state, Temp, node, ibatch)
    end
    state
end
function step!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::CuMatrix, Temp::Float64, node::Int)
    function kernel(rule, sa, Temp, node)
        ibatch = threadIdx().x
        step_kernel!(rule, sa, state, Temp, node, ibatch)
        return nothing
    end
    @cuda threads = size(sa, 2) kernel(rule, sa, Temp, node)
    state
end

function step_kernel!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, Temp::Float64, node::Int, ibatch::Int)
    ΔE = 0
    i, j = linear_to_cartesian(sa, node).I
    if j > 1 # not the first layer
        ΔE -= 2 * evaluate_parent(sa, state, node, ibatch)
    end
    if j < sa.m # not the last layer
        cnodes = child_nodes(sa, node)
        for node in cnodes
            ΔE -= evaluate_parent(sa, state, node, ibatch)
        end
        # flip the node
        @inbounds state[node, ibatch] ⊻= true
        for node in cnodes
            ΔE += evaluate_parent(sa, state, node, ibatch)
        end
    end
    if rand() < prob_accept(rule, Temp, ΔE)
        ΔE
    else  # revert the flip
        @inbounds state[node, ibatch] ⊻= true
        0
    end
end
prob_accept(::Metropolis, Temp, ΔE::T) where T<:Real = ΔE < 0 ? 0.0 : exp(- (ΔE) / Temp)
prob_accept(::HeatBath, Temp, ΔE::Real) = inv(1 + exp(ΔE / Temp))

function track_equilibration!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, tempscale = 4 .- (1:100 .-1) * 0.04)
    # NOTE: do we really need niters? or just set it to 1?
    for Temp in tempscale
        # NOTE: do we really need to shuffle the nodes?
        for node in 1:natom(sa)
            step!(rule, sa, state, Temp, node)
        end
    end
    return sa
end