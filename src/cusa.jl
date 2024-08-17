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

# evaluate the energy of the i-th gadget (involving atoms i and its parents)
function evaluate_parent(sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, inode::Integer, ibatch::Integer)
    i, j = CartesianIndices((sa.n, sa.m))[inode].I
    idp = parent_nodes(sa, inode)
    # trueoutput = @inbounds rule110(state[idp[1], ibatch], state[idp[2], ibatch], state[idp[3], ibatch])
    trueoutput = @inbounds automatarule(state[idp[1], ibatch], state[idp[2], ibatch], state[idp[3], ibatch], 150)
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

function step!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, Temp, node=nothing)
    for ibatch in 1:size(state, 2)
        step_kernel!(rule, sa, state, energy_gradient, Temp, ibatch, node)
    end
    state
end
function step!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::CuMatrix, energy_gradient::AbstractArray, Temp, node=nothing)
    @inline function kernel(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, Temp, node=nothing)
        ibatch = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
        if ibatch <= size(state, 2)
            step_kernel!(rule, sa, state, energy_gradient, Temp, ibatch, node)
        end
        return nothing
    end
    kernel = @cuda launch=false kernel(rule, sa, state, energy_gradient, Temp, node)
    config = launch_configuration(kernel.fun)
    threads = min(size(state, 2), config.threads)
    blocks = cld(size(state, 2), threads)
    CUDA.@sync kernel(rule, sa, state, energy_gradient, Temp, node; threads, blocks)
    state
end

@inline function step_kernel!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state, energy_gradient::AbstractArray, Temp, ibatch::Integer, node=nothing)
    ΔE_with_next_layer = 0
    ΔE_with_previous_layer = 0
    # node = rand(atoms(sa))
    if node === nothing
        node = rand(atoms(sa))
    end
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
    if j == sa.m
        flip_max_prob *= prob_accept(rule, Temp[ibatch][j-1], ΔE_with_previous_layer)
        # flip_max_prob = prob_accept(HeatBath(), Temp[ibatch][j-1], ΔE_with_previous_layer)
    elseif j == 1
        flip_max_prob *= prob_accept(rule, Temp[ibatch][j], ΔE_with_next_layer)
        # flip_max_prob *= prob_accept(HeatBath(), Temp[ibatch][j], ΔE_with_next_layer)
    else
        flip_max_prob = 1.0 / (1.0 + exp(ΔE_with_previous_layer / Temp[ibatch][j-1] + ΔE_with_next_layer / Temp[ibatch][j]))
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



function step_parallel!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, Temp, flip_id)
    # @info "flip_id = $flip_id"
    for ibatch in 1:size(state, 2)
        for this_time_flip in flip_id
            step_kernel!(rule, sa, state, energy_gradient, Temp, ibatch, this_time_flip)
        end
    end
    state
end
function step_parallel!(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::CuMatrix, energy_gradient::AbstractArray, Temp, flip_id)
    @inline function kernel(rule::TransitionRule, sa::SimulatedAnnealingHamiltonian, state::AbstractMatrix, energy_gradient::AbstractArray, Temp, flip_id)
        id = (blockIdx().x - Int32(1)) * blockDim().x + threadIdx().x
        stride = blockDim().x * gridDim().x
        Nx = size(state, 2)
        Ny = length(flip_id)
        cind = CartesianIndices((Nx, Ny))
        for k in id:stride:Nx*Ny
            ibatch = cind[k][1]
            id = cind[k][2]
            step_kernel!(rule, sa, state, energy_gradient, Temp, ibatch, flip_id[id])
        end
        return nothing
    end
    kernel = @cuda launch=false kernel(rule, sa, state, energy_gradient, Temp, flip_id)
    config = launch_configuration(kernel.fun)
    threads = min(size(state, 2) * length(flip_id), config.threads)
    blocks = cld(size(state, 2) * length(flip_id), threads)
    CUDA.@sync kernel(rule, sa, state, energy_gradient, Temp, flip_id; threads, blocks)
    state
end


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
    eachposition = Tuple([gausspulse_amplitude * gradient^(- (1.0/gausspulse_width) * abs(i-middle_position)) + 1e-5 for i in 1:sa.m-1])
    return eachposition
end

function get_parallel_flip_id(sa)
    ret = Vector{Vector{Int}}()
    for cnt in 1:6
        temp = Vector{Int}()
        for layer in 1+div(cnt - 1, 3):2:sa.m
            for position in mod1(cnt, 3):3:(sa.n - sa.n % 3)
                push!(temp, LinearIndices((sa.n, sa.m))[position, layer])
            end
        end
        push!(ret, temp)
    end
    if sa.n % 3 >= 1
        push!(ret, Vector(sa.n:2*sa.n:sa.n*sa.m))
        push!(ret, Vector(2*sa.n:2*sa.n:sa.n*sa.m))
    end
    if sa.n % 3 >= 2
        push!(ret, Vector(sa.n-1:2*sa.n:sa.n*sa.m))
        push!(ret, Vector(2*sa.n-1:2*sa.n:sa.n*sa.m))
    end
    return ret
end


function track_equilibration_gausspulse_cpu!(rule::TransitionRule,
                                        sa::SimulatedAnnealingHamiltonian, 
                                        state::AbstractMatrix, 
                                        energy_gradient, 
                                        gausspulse_amplitude::Float64,
                                        gausspulse_width::Float64,
                                        annealing_time; accelerate_flip = false
                                        )
    midposition = 1.0 - (-(gausspulse_width) * log(1e-5/gausspulse_amplitude) / log(energy_gradient))
    each_movement = ((1.0 - midposition) * 2 + (sa.m - 2)) / (annealing_time - 1)
    # @info "midposition = $midposition"
    # @info "each_movement = $each_movement"

    # NOTE: do we really need niters? or just set it to 1?
    single_layer_temp = []
    for t in 1:annealing_time
        # @info "midposition = $midposition"
        singlebatch_temp = toymodel_gausspulse(sa, gausspulse_amplitude, gausspulse_width, midposition, energy_gradient)
        Temp = fill((singlebatch_temp), size(state, 2))
        if accelerate_flip == false
            for thisatom in 1:natom(sa)
                step!(rule, sa, state, fill(1.0, size(state, 2)), Temp, thisatom)
            end
        else
            flip_list = get_parallel_flip_id(sa)
            for eachflip in flip_list
                step_parallel!(rule, sa, state, fill(1.0, size(state, 2)), Temp, eachflip)
            end
        end
        midposition += each_movement
        # push!(single_layer_temp, singlebatch_temp[1])
    end
    # return single_layer_temp
end

function track_equilibration_gausspulse_gpu!(rule::TransitionRule,
                                        sa::SimulatedAnnealingHamiltonian, 
                                        state::AbstractMatrix, 
                                        energy_gradient, 
                                        gausspulse_amplitude,
                                        gausspulse_width,
                                        annealing_time; accelerate_flip = false
                                        )
    midposition = 1.0 - (-(gausspulse_width) * log(1e-5/gausspulse_amplitude) / log(energy_gradient))
    each_movement = ((1.0 - midposition) * 2 + (sa.m - 2)) / (annealing_time - 1)
    @info "each_movement = $each_movement"

    # NOTE: do we really need niters? or just set it to 1?
    for t in 1:annealing_time
        singlebatch_temp = toymodel_gausspulse(sa, gausspulse_amplitude, gausspulse_width, midposition, energy_gradient)
        Temp = CuArray(fill(Float32.(singlebatch_temp), size(state, 2)))
        if accelerate_flip == false
            for thisatom in 1:natom(sa)
                step!(rule, sa, state, CuArray(fill(1.0, size(state, 2))), Temp, thisatom)
            end
        else
            flip_list = get_parallel_flip_id(sa)
            for eachflip in flip_list
                step_parallel!(rule, sa, state, CuArray(fill(1.0, size(state, 2))), Temp, CuArray(eachflip))
            end
        end
        midposition += each_movement
    end
    return sa
end


function track_equilibration_gausspulse_reverse_cpu!(rule::TransitionRule,
                                        sa::SimulatedAnnealingHamiltonian, 
                                        state::AbstractMatrix, 
                                        energy_gradient, 
                                        gausspulse_amplitude::Float64,
                                        gausspulse_width::Float64,
                                        annealing_time;
                                        accelerate_flip = false
                                        )
    midposition = 1.0 - (-(gausspulse_width) * log(1e-5/gausspulse_amplitude) / log(energy_gradient))
    each_movement = ((1.0 - midposition) * 2 + (sa.m - 2)) / (annealing_time - 1)
    midposition = sa.m - 1.0 + 1.0 - midposition

    single_layer_temp = []
    for t in 1:annealing_time
        @info "midposition = $midposition"
        singlebatch_temp = toymodel_gausspulse(sa, gausspulse_amplitude, gausspulse_width, midposition, energy_gradient)
        Temp = fill(singlebatch_temp, size(state, 2))
        if accelerate_flip == false
            for thisatom in 1:natom(sa)
                step!(rule, sa, state, fill(1.0, size(state, 2)), Temp, thisatom)
            end
        else
            flip_list = get_parallel_flip_id(sa)
            for eachflip in flip_list
                step_parallel!(rule, sa, state, fill(1.0, size(state, 2)), Temp, eachflip)
            end
        end
        midposition -= each_movement
        # push!(single_layer_temp, singlebatch_temp[1])
    end
    # return single_layer_temp
end