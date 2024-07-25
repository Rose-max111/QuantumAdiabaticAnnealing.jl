struct SimulatedAnnealingHamiltonian
    n::Int # number of atoms per layer
    m::Int # number of full layers (count the last layer!)
end
natom(sa::SimulatedAnnealingHamiltonian) = sa.n * sa.m
atoms(sa::SimulatedAnnealingHamiltonian) = Base.OneTo(natom(sa))
function random_state(sa::SimulatedAnnealingHamiltonian, nbatch::Integer)
    return rand(Bool, natom(sa), nbatch)
end
hasparent(sa::SimulatedAnnealingHamiltonian, node::Integer) = node > sa.n
function linear_to_cartesian(sa::SimulatedAnnealingHamiltonian, node::Integer)
    j, i = divrem(node-1, sa.n)
    return i + 1, j + 1
end
cartesian_to_linear(sa::SimulatedAnnealingHamiltonian, i::Integer, j::Integer) = i + (j - 1) * sa.n

# evaluate the energy of the i-th gadget (involving atoms i and its parents)
function evaluate_parent(sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, inode::Integer, ibatch::Integer)
    idp = parent_logic(sa, inode)
    trueoutput = @inbounds rule110(state[idp[1], ibatch], state[idp[2], ibatch], state[idp[3], ibatch])
    return @inbounds trueoutput ⊻ state[idp[4], ibatch]
end
function calculate_energy(sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, ibatch::Integer)
    return sum(i->evaluate_parent(sa, state, i, ibatch), sa.n+1:natom(sa))
end
function parent_logic(sa::SimulatedAnnealingHamiltonian, node::Integer)
    n = sa.n
    i, j = linear_to_cartesian(sa, node)
    (
        cartesian_to_linear(sa, mod1(i-1, n), j-1),  # periodic boundary condition
        cartesian_to_linear(sa, i, j-1),
        cartesian_to_linear(sa, mod1(i+1, n), j-1),
        node
    )
end

function child_nodes(sa::SimulatedAnnealingHamiltonian, node::Integer)
    n = sa.n
    i, j = linear_to_cartesian(sa, node)
    (
        cartesian_to_linear(sa, mod1(i-1, n), j+1),  # periodic boundary condition
        cartesian_to_linear(sa, mod1(i, n), j+1),
        cartesian_to_linear(sa, mod1(i+1, n), j+1),
    )
end

abstract type TransitionRule end
struct HeatBath <: TransitionRule end
struct Metropolis <: TransitionRule end

function step!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, Temp::Float64, node::Integer)
    for ibatch in 1:size(state, 2)
        step_kernel!(rule, sa, state, Temp, node, ibatch)
    end
    state
end
function step!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::CuMatrix, Temp::Float64, node::Integer)
    @inline function kernel(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, Temp::Float64, node::Integer)
        ibatch = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
        if ibatch <= size(state, 2)
            step_kernel!(rule, sa, state, Temp, node, ibatch)
        end
        return nothing
    end
    kernel = @cuda launch=false kernel(rule, sa, state, Temp, node)
    config = launch_configuration(kernel.fun)
    threads = min(size(state, 2), config.threads)
    blocks = cld(size(state, 2), threads)
    CUDA.@sync kernel(rule, sa, state, Temp, node; threads, blocks)
    state
end

@inline function step_kernel!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state, Temp::Float64, node::Integer, ibatch::Integer)
    ΔE = 0
    i, j = linear_to_cartesian(sa, node)
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