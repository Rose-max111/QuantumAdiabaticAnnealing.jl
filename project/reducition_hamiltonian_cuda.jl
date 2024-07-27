using QuantumAdiabaticAnnealing: step!, track_equilibration!, SimulatedAnnealingHamiltonian, random_state, Metropolis, HeatBath
using QuantumAdiabaticAnnealing: calculate_energy, toy_model_state_energy, toy_model_transition_matrix, TransitionRule
using KrylovKit
using CUDA
using BenchmarkTools

# sa = SimulatedAnnealingHamiltonian(40, 100)
# nbatch = 10000

# state = CuArray(random_state(sa, nbatch))
# @btime CUDA.@sync step!(HeatBath(), $sa, $state, 1.0, 1)

# cpu_state = random_state(sa, nbatch)
# @btime step!(HeatBath(), $sa, $cpu_state, 1.0, 1)
#@test track_equilibration!(Metropolis(), sa, state) isa SimulatedAnnealingHamiltonian


sa = SimulatedAnnealingHamiltonian(3, 4)
nbatch = 5000

state = CuArray(random_state(sa, nbatch))
# state = random_state(sa, nbatch)
anneal_time = 1
# begin_temp = 20.0
end_temp = 2.0
energy_gradient = CuArray(fill(3.0f0, anneal_time))
@btime track_equilibration!(HeatBath(), $sa, $state, energy_gradient, fill(end_temp, anneal_time))

cpu_state = Array(state)
state_10 = [sum(cpu_state[:,i] .* [2^(j-1) for j in 1:sa.n*sa.m]) for i in 1:size(cpu_state)[2]]
figure_count = [count(==(i), state_10) for i in 0:2^(sa.n*sa.m)-1]
distribution_count = 1.0 .* figure_count ./ sum(figure_count)
state_energy = [calculate_energy(sa, cpu_state, energy_gradient, i) for i in 1:nbatch]
state_energy_average = 1.0 * sum(state_energy) / nbatch


n = 3
m = 4
Temp = end_temp
P = toy_model_transition_matrix(HeatBath(), n, m, Temp; period_condition = true, energy_gradient = 3)
eigvals, eigvecs, infos = eigsolve(P, rand(Float64, size(P)[1]), 2, :LR; maxiter = 5000)

# calculate the exact boltzmann distribution corresponding energy
ED_state_energy = [toy_model_state_energy(i, n, m; period_condition = true, energy_gradient = 3) for i in 0:2^(n*m)-1]
maxeigvec_energy = sum(ED_state_energy .* (abs.(Real.(eigvecs[1])))) / abs(sum(Real.(eigvecs[1])))


function evaluate_50percent_time(rule::TransitionRule, width::Integer, depth::Integer, energy_gradient)
    sa = SimulatedAnnealingHamiltonian(width, depth)
    nbatch = 5000
    Temp_start = 5 * energy_gradient^(sa.m-2)
    Temp_end = 1e-5

    anneal_time = 0
    max_try = 1
    while max_try > 0
        next_try = anneal_time + max_try
        η = (Temp_end / Temp_start) ^ (1/next_try)
        energy_gradient_sa = CuArray(fill(Float32(energy_gradient), next_try))
        Temp_sa = Float64.([Temp_start * η^i for i in 0:next_try-1])
        
        state = CuArray(random_state(sa, nbatch))
        @info "Stage1, max_try = $max_try, next_try = $next_try, begin annealing"
        CUDA.@sync track_equilibration!(HeatBath(), sa, state, energy_gradient_sa, Temp_sa)
        
        cpu_state = Array(state)
        state_energy = [calculate_energy(sa, cpu_state, energy_gradient, i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        @info "Stage 1, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            break
        else
            continue
            # anneal_time = next_try
        end
        max_try *= 2
    end
    while max_try != 1
        max_try /= 2

        next_try = anneal_time + max_try
        η = (Temp_end / Temp_start) ^ (1/next_try)
        energy_gradient_sa = CuArray(fill(Float32(energy_gradient), next_try))
        Temp_sa = Float64.([Temp_start * η^i for i in 0:next_try-1])
        
        @info "Stage2, max_try = $max_try, next_try = $next_try, begin annealing"
        state = CuArray(random_state(sa, nbatch))
        CUDA.@sync track_equilibration!(HeatBath(), sa, state, energy_gradient_sa, Temp_sa)
        
        cpu_state = Array(state)
        state_energy = [calculate_energy(sa, cpu_state, energy_gradient, i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        @info "Stage 2, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            continue
        else
            anneal_time = next_try
        end
    end
    return anneal_time
end

evaluate_50percent_time(HeatBath(), 3, 4, 2.0)