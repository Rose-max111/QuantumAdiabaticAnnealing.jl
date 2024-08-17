using QuantumAdiabaticAnnealing
using QuantumAdiabaticAnnealing: random_state, calculate_energy
using QuantumAdiabaticAnnealing: track_equilibration_gausspulse_gpu!, track_equilibration_gausspulse_reverse_cpu!
using CUDA
using QuantumAdiabaticAnnealing: get_parallel_flip_id

function evaluate_50percent_time_cpu(width::Integer, depth::Integer, gauss_width, energy_gradient)
    sa = SimulatedAnnealingHamiltonian(width, depth)
    nbatch = 5000

    anneal_time = 0
    max_try = 2
    while max_try > 0
        next_try = anneal_time + max_try
        
        state = random_state(sa, nbatch)
        @info "Stage1, max_try = $max_try, next_try = $next_try, begin annealing"
        # @info "$Temp_sa"
        @time track_equilibration_gausspulse_cpu!(HeatBath(), sa, state, energy_gradient, 10.0, gauss_width, next_try; accelerate_flip = true)
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
        @time track_equilibration_gausspulse_cpu!(HeatBath(), sa, state, energy_gradient, 10.0, gauss_width, next_try; accelerate_flip = true)
        @info "finish annealing"

        state_energy = [calculate_energy(sa, state, fill(1.0, nbatch), i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        # @info "Stage 2, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            max_try /=2
            continue
        else
            max_try /=2
            anneal_time = next_try
        end
    end
    return anneal_time
end


function evaluate_50percent_time_gpu(width::Integer, depth::Integer, gauss_width, energy_gradient)
    sa = SimulatedAnnealingHamiltonian(width, depth)
    nbatch = 5000

    anneal_time = 0.0
    max_try = 2.0
    while max_try > 0
        next_try = anneal_time + max_try
        
        state = CuArray(random_state(sa, nbatch))
        @info "Stage1, max_try = $max_try, next_try = $next_try, begin annealing"
        # @info "$Temp_sa"
        @time track_equilibration_gausspulse_gpu!(HeatBath(), sa, state, energy_gradient, 10.0, gauss_width, next_try; accelerate_flip = true)
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
        @time track_equilibration_gausspulse_gpu!(HeatBath(), sa, state, energy_gradient, 10.0, gauss_width, next_try; accelerate_flip = true)
        @info "finish annealing"

        cpu_state = Array(state)
        state_energy = [calculate_energy(sa, cpu_state, fill(1.0, nbatch), i) for i in 1:nbatch]
        success = count(x -> x == 0, state_energy)
        # @info "Stage 2, now max_try = $max_try, success time = $success, anneal_time = $anneal_time"
        if 1.0 * success / nbatch >= 0.49
            max_try /=2
            continue
        else
            max_try /=2
            anneal_time = next_try
        end
    end
    return anneal_time + 1
end

# evaluate_50percent_time_cpu(12, 8, 1.0, 1.5)

sa = SimulatedAnnealingHamiltonian(10, 10)
nbatch = 1000
energy_gradient = 1.4
state = random_state(sa, nbatch)
# state_energy = [calculate_energy(sa, state, fill(1.0, nbatch), i) for i in 1:nbatch]

# single_layer_temp = track_equilibration_gausspulse_reverse_cpu!(HeatBath(), sa, state, energy_gradient, 10.0, 1.0, 800; accelerate_flip = true)

single_layer_temp = track_equilibration_gausspulse_cpu!(HeatBath(), sa, state, energy_gradient, 10.0, 1.0, 800; accelerate_flip = true)

state_energy = [calculate_energy(sa, state, fill(1.0, nbatch), i) for i in 1:nbatch]
success = count(x -> x == 0, state_energy)
