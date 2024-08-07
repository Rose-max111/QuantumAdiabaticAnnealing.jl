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
cartesian_to_linear(sa::SimulatedAnnealingHamiltonian, i::Integer, j::Integer) = i + (j - 1) * sa.n

# evaluate the energy of the i-th gadget (involving atoms i and its parents)
function evaluate_parent(sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, inode::Integer, ibatch::Integer)
    i, j = CartesianIndices((sa.n, sa.m))[inode].I
    idp = parent_nodes(sa, inode)
    trueoutput = @inbounds rule110(state[idp[1], ibatch], state[idp[2], ibatch], state[idp[3], ibatch])
    return @inbounds (trueoutput ⊻ state[inode, ibatch]) * (energy_gradient[ibatch] ^ (sa.m - j))
end
function calculate_energy(sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, ibatch::Integer)
    return sum(i->evaluate_parent(sa, state, energy_gradient, i, ibatch), sa.n+1:natom(sa))
end

function local_energy!(sa::SimulatedAnnealingHamiltonian, state::CuMatrix, energy_gradient::CuVector, energy::CuVector)
    @inline function kernel(sa::SimulatedAnnealingHamiltonian, state, energy_gradient, energy)
        ibatch = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
        if ibatch <= size(state, 2)
            # step_kernel!(rule, sa, state, energy_gradient, Temp, ibatch)
            for i in sa.n+1:natom(sa)
                energy[ibatch] += evaluate_parent(sa, state, energy_gradient, i, ibatch)
            end
        end
        return nothing
    end
    kernel = @cuda launch=false kernel(sa, state, energy_gradient, energy)
    config = launch_configuration(kernel.fun)
    threads = min(size(state, 2), config.threads)
    blocks = cld(size(state, 2), threads)
    CUDA.@sync kernel(sa, state, energy_gradient, energy; threads, blocks)
    energy
end

function parent_nodes(sa::SimulatedAnnealingHamiltonian, node::Integer)
    n = sa.n
    i, j = CartesianIndices((n, sa.m))[node].I
    lis = LinearIndices((n, sa.m))
    @inbounds (
        lis[mod1(i-1, n), j-1],  # periodic boundary condition
        lis[i, j-1],
        lis[mod1(i+1, n), j-1],
    )
end

function child_nodes(sa::SimulatedAnnealingHamiltonian, node::Integer)
    n = sa.n
    i, j = CartesianIndices((n, sa.m))[node].I
    lis = LinearIndices((n, sa.m))
    @inbounds (
        lis[mod1(i-1, n), j+1],  # periodic boundary condition
        lis[mod1(i, n), j+1],
        lis[mod1(i+1, n), j+1],
    )
end

function step!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, Temp)
    for ibatch in 1:size(state, 2)
        step_kernel!(rule, sa, state, energy_gradient, Temp, ibatch)
    end
    state
end
function step!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::CuMatrix, energy_gradient::AbstractArray, Temp)
    @inline function kernel(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, Temp)
        ibatch = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
        if ibatch <= size(state, 2)
            step_kernel!(rule, sa, state, energy_gradient, Temp, ibatch)
        end
        return nothing
    end
    kernel = @cuda launch=false kernel(rule, sa, state, energy_gradient, Temp)
    config = launch_configuration(kernel.fun)
    threads = min(size(state, 2), config.threads)
    blocks = cld(size(state, 2), threads)
    CUDA.@sync kernel(rule, sa, state, energy_gradient, Temp; threads, blocks)
    state
end

@inline function step_kernel!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state, energy_gradient::AbstractArray, Temp, ibatch::Integer)
    ΔE_with_next_layer = 0
    ΔE_with_previous_layer = 0
    node = rand(atoms(sa))

    i, j = CartesianIndices((sa.n, sa.m))[node].I
    if j > 1 # not the first layer
        # ΔE += energy_gradient[ibatch]^(sa.m - j) - 2 * evaluate_parent(sa, state, energy_gradient, node, ibatch)
        ΔE_with_previous_layer += energy_gradient[ibatch]^(sa.m - j) - 2 * evaluate_parent(sa, state, energy_gradient, node, ibatch)
    end
    if j < sa.m # not the last layer
        cnodes = child_nodes(sa, node)
        for node in cnodes
            ΔE_with_next_layer -= evaluate_parent(sa, state, energy_gradient, node, ibatch)
        end
        # flip the node@
        @inbounds state[node, ibatch] ⊻= true
        for node in cnodes
            ΔE_with_next_layer += evaluate_parent(sa, state, energy_gradient, node, ibatch)
        end
        @inbounds state[node, ibatch] ⊻= true
    end
    flip_max_prob = 1
    if j < sa.m
        flip_max_prob *= prob_accept(rule, Temp[ibatch][j], ΔE_with_next_layer)
    end
    if j > 1
        flip_max_prob *= prob_accept(rule, Temp[ibatch][j-1], ΔE_with_previous_layer)
    end
    if rand() < flip_max_prob
        @inbounds state[node, ibatch] ⊻= true
        ΔE_with_next_layer
    else
        0
    end
end
prob_accept(::Metropolis, Temp, ΔE::T) where T<:Real = ΔE < 0 ? 1.0 : exp(- (ΔE) / Temp)
prob_accept(::HeatBath, Temp, ΔE::Real) = inv(1 + exp(ΔE / Temp))

function track_equilibration!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, tempscale = 4 .- (1:100 .-1) * 0.04)
    # NOTE: do we really need niters? or just set it to 1?
    for Temp in tempscale
        # NOTE: do we really need to shuffle the nodes?
        for _ in 1:natom(sa)
        # for node in shuffle!(Vector(1:natom(sa)))
            step!(rule, sa, state, energy_gradient, Temp)
        end
    end
    return sa
end

function toymodel_gausspulse(sa::SimulatedAnnealingHamiltonian,
                            gausspulse_amplitude::Float64,
                            gausspulse_width::Float64,
                            middle_position::Float64,
                            gradient::Float64)
    # amplitude * e^(- (1 /width) * (x-middle_position)^2)
    eachposition = Tuple([gausspulse_amplitude * gradient^(- (1.0/gausspulse_width) * (i-middle_position)^2) + 1e-4 for i in 1:sa.m-1])
    return eachposition
end

function track_equilibration_gausspulse_cpu!(rule::TransitionRule,
                                        sa::SimulatedAnnealingHamiltonian, 
                                        state::AbstractMatrix, 
                                        energy_gradient::AbstractArray, 
                                        gausspulse_amplitude::Float64,
                                        gausspulse_width::Float64,
                                        annealing_time
                                        )
    midposition = 1.0 - sqrt(-(gausspulse_width) * log(1e-5/gausspulse_amplitude)) / log(energy_gradient[1])
    each_movement = ((1.0 - midposition) * 2 + (sa.m - 2)) / (annealing_time - 1)
    # @info "midposition = $midposition"
    # @info "each_movement = $each_movement"

    # NOTE: do we really need niters? or just set it to 1?
    for t in 1:annealing_time
        singlebatch_temp = toymodel_gausspulse(sa, gausspulse_amplitude, gausspulse_width, midposition, energy_gradient[1])
        Temp = fill(singlebatch_temp, size(state, 2))
        for _ in 1:natom(sa)
            step!(rule, sa, state, energy_gradient, Temp)
        end
        midposition += each_movement
        # if t == annealing_time / 2
        #     @info "medium process temp = $singlebatch_temp"
        # end
    end
    return sa
end

function track_equilibration_gausspulse_gpu!(rule::TransitionRule,
                                        sa::SimulatedAnnealingHamiltonian, 
                                        state::AbstractMatrix, 
                                        energy_gradient::AbstractArray, 
                                        gausspulse_amplitude::Float32,
                                        gausspulse_width::Float32,
                                        annealing_time
                                        )
    midposition = 1.0 - sqrt(-(gausspulse_width) * log(1e-5/gausspulse_amplitude)) / log(energy_gradient[1])
    each_movement = ((1.0 - midposition) * 2 + (sa.m - 2)) / (annealing_time - 1)
    @info "each_movement = $each_movement"

    # NOTE: do we really need niters? or just set it to 1?
    for t in 1:annealing_time
        singlebatch_temp = toymodel_gausspulse(sa, gausspulse_amplitude, gausspulse_width, midposition, energy_gradient[1])
        Temp = CuArray(fill(Float32.(singlebatch_temp), size(state, 2)))
        for _ in 1:natom(sa)
            step!(rule, sa, state, energy_gradient, Temp)
        end
        midposition += each_movement
    end
    return sa
end