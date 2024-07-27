using QuantumAdiabaticAnnealing: step!, track_equilibration!, SimulatedAnnealingHamiltonian, random_state, Metropolis, HeatBath
using QuantumAdiabaticAnnealing: calculate_energy, toy_model_state_energy, toy_model_transition_matrix, TransitionRule
using KrylovKit
using CUDA; CUDA.allowscalar(false)
using BenchmarkTools

# sa = SimulatedAnnealingHamiltonian(40, 100)
# nbatch = 10000

# state = CuArray(random_state(sa, nbatch))
# @btime CUDA.@sync step!(HeatBath(), $sa, $state, 1.0, 1)

# cpu_state = random_state(sa, nbatch)
# @btime step!(HeatBath(), $sa, $cpu_state, 1.0, 1)
#@test track_equilibration!(Metropolis(), sa, state) isa SimulatedAnnealingHamiltonian

function evaluate_50percent_time(rule::TransitionRule, width::Integer, depth::Integer, energy_gradient)
    sa = SimulatedAnnealingHamiltonian(width, depth)
    nbatch = 5000
    Temp_start = 5 * energy_gradient^(sa.m-2)
    Temp_end = 1e-5
    energy_gradient_sa = CuArray(fill(Float32(energy_gradient), nbatch))

    anneal_time = 0
    max_try = 1
    while max_try > 0
        next_try = anneal_time + max_try
        η = (Temp_end / Temp_start) ^ (1/(next_try-1))
        Temp_sa = Float64.([Temp_start * η^i for i in 0:next_try-1])
        
        state = CuArray(random_state(sa, nbatch))
        # @info "Stage1, max_try = $max_try, next_try = $next_try, begin annealing"
        # @info "$Temp_sa"
        @time CUDA.@sync track_equilibration!(HeatBath(), sa, state, energy_gradient_sa, Temp_sa)

        cpu_state = Array(state)
        state_energy = [calculate_energy(sa, cpu_state, Array(energy_gradient_sa), i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        # @info "Stage 1, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            break
        else
            max_try *= 2
        end
    end
    while max_try != 1
        max_try /= 2

        next_try = anneal_time + max_try
        η = (Temp_end / Temp_start) ^ (1/(next_try-1))
        Temp_sa = Float64.([Temp_start * η^i for i in 0:next_try-1])
        
        @info "Stage2, max_try = $max_try, next_try = $next_try, begin annealing"
        state = CuArray(random_state(sa, nbatch))
        CUDA.@sync track_equilibration!(HeatBath(), sa, state, energy_gradient_sa, Temp_sa)
        @info "finish annealing"
        
        cpu_state = Array(state)
        state_energy = [calculate_energy(sa, cpu_state, Array(energy_gradient_sa), i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        # @info "Stage 2, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            continue
        else
            anneal_time = next_try
        end
    end
    return anneal_time
end

width = ARGS[1]
depth = ARGS[2]
λ = ARGS[3]
device = ARGS[4]

width = parse(Int, width)
depth = parse(Int, depth)
λ = parse(Float32, λ)
device = parse(Int, device)

CUDA.device!(device)
evaluate_time = evaluate_50percent_time(HeatBath(), width, depth, λ)

@info "width = $width, depth = $depth, λ = $λ, evaluate_time = $evaluate_time"

filepath = joinpath(@__DIR__, "data_toymodel/W=$(width)_D=$(depth)_E=$(λ).txt")
open(filepath,"w") do file
    println(file, evaluate_time)
end

# depth = Vector(2:2:20)
# evaluate_time = Vector{Int}()
# for D in depth
#     push!(evaluate_time, evaluate_50percent_time(HeatBath(), 20, D, 3))
#     @info "depth = $D, time = $(evaluate_time[end])"
# end

# time_1o5 = copy(evaluate_time)
# depth_1o5 = depth[1:7]

# time_3 = copy(evaluate_time)
# depth_3 = copy(depth)


# using CairoMakie
# f = Figure()
# ax = Axis(f[1, 1], xlabel = "depth", ylabel = "scan time",
#     title = "Time v.s. Depth(50% success probability)")
# scatter!(ax, depth_1o5, time_1o5, label = "λ = 1.5")
# scatter!(ax, depth_3, time_3 .* 40, label = "λ = 3(rescale with 40)")
# f[1, 2] = Legend(f, ax)
# f