using QuantumAdiabaticAnnealing
using QuantumAdiabaticAnnealing: random_state, calculate_energy
using QuantumAdiabaticAnnealing: track_equilibration_gausspulse_reverse_cpu!

function evaluate_50percent_time_cpu(width::Integer, depth::Integer, gauss_width, energy_gradient)
    sa = SimulatedAnnealingHamiltonian(width, depth)
    nbatch = 300
    energy_gradient_sa = fill(Float64(energy_gradient), nbatch)

    anneal_time = 0
    max_try = 1
    while max_try > 0
        next_try = anneal_time + max_try
        
        state = random_state(sa, nbatch)
        @info "Stage1, max_try = $max_try, next_try = $next_try, begin annealing"
        # @info "$Temp_sa"
        @time track_equilibration_gausspulse_cpu!(HeatBath(), sa, state, energy_gradient_sa, 10.0, gauss_width, next_try)
        @info "finish annealing"

        state_energy = [calculate_energy(sa, state, fill(1.0, nbatch), i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        @info "Stage 1, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            break
        else
            max_try *= 2
        end
    end
    anneal_time += max_try / 2
    max_try /= 2
    if max_try <= 1
        return anneal_time
    end
    while max_try != 1
        next_try = anneal_time + max_try
        
        state = random_state(sa, nbatch)
        @info "Stage2, max_try = $max_try, next_try = $next_try, begin annealing"
        # @info "$Temp_sa"
        @time track_equilibration_gausspulse_cpu!(HeatBath(), sa, state, energy_gradient_sa, 10.0, gauss_width, next_try)
        @info "finish annealing"

        state_energy = [calculate_energy(sa, state, fill(1.0, nbatch), i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        # @info "Stage 2, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            continue
        else
            anneal_time = next_try
        end
        max_try /=2
    end
    return anneal_time
end


function evaluate_50percent_time_gpu(width::Integer, depth::Integer, gauss_width, energy_gradient)
    sa = SimulatedAnnealingHamiltonian(width, depth)
    nbatch = 300
    energy_gradient_sa = CuArray(fill(Float32(energy_gradient), nbatch))

    anneal_time = 0
    max_try = 1
    while max_try > 0
        next_try = anneal_time + max_try
        
        state = CuArray(random_state(sa, nbatch))
        @info "Stage1, max_try = $max_try, next_try = $next_try, begin annealing"
        # @info "$Temp_sa"
        @time track_equilibration_gausspulse_gpu!(HeatBath(), sa, state, energy_gradient_sa, 10.0, gauss_width, next_try)
        @info "finish annealing"

        cpu_state = Array(state)
        state_energy = [calculate_energy(sa, cpu_state, fill(1.0, nbatch), i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        @info "Stage 1, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            break
        else
            max_try *= 2
        end
    end
    anneal_time += max_try / 2
    max_try /= 2
    if max_try <= 1
        return anneal_time
    end
    while max_try != 1
        next_try = anneal_time + max_try
        
        state = CuArray(random_state(sa, nbatch))
        @info "Stage2, max_try = $max_try, next_try = $next_try, begin annealing"
        # @info "$Temp_sa"
        @time track_equilibration_gausspulse_gpu!(HeatBath(), sa, state, energy_gradient_sa, 10.0, gauss_width, next_try)
        @info "finish annealing"

        cpu_state = Array(state)
        state_energy = [calculate_energy(sa, cpu_state, fill(1.0, nbatch), i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        # @info "Stage 2, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            continue
        else
            anneal_time = next_try
        end
        max_try /=2
    end
    return anneal_time
end


# sa = SimulatedAnnealingHamiltonian(6, 12)
# nbatch = 100
# energy_gradient_sa = fill(1.2, nbatch)
# state = random_state(sa, nbatch)
# state_energy = [calculate_energy(sa, state, energy_gradient_sa, i) for i in 1:nbatch]

# # track_equilibration_gausspulse_reverse_cpu!(HeatBath(), sa, state, energy_gradient_sa, 10.0, 2.0, 500)

# track_equilibration_gausspulse_cpu!(HeatBath(), sa, state, energy_gradient_sa, 10.0, 20.0, 10000)

# state_energy = [calculate_energy(sa, state, fill(1.0, nbatch), i) for i in 1:nbatch]
# success = count(x -> x == 0, state_energy)