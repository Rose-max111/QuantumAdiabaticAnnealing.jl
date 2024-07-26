using QuantumAdiabaticAnnealing: step!, track_equilibration!, SimulatedAnnealingHamiltonian, random_state, Metropolis, HeatBath
using QuantumAdiabaticAnnealing: calculate_energy
using CUDA
using BenchmarkTools

# sa = SimulatedAnnealingHamiltonian(40, 100)
# nbatch = 10000

# state = CuArray(random_state(sa, nbatch))
# @btime CUDA.@sync step!(HeatBath(), $sa, $state, 1.0, 1)

# cpu_state = random_state(sa, nbatch)
# @btime step!(HeatBath(), $sa, $cpu_state, 1.0, 1)
#@test track_equilibration!(Metropolis(), sa, state) isa SimulatedAnnealingHamiltonian


sa = SimulatedAnnealingHamiltonian(8, 2)
nbatch = 5000

state = CuArray(random_state(sa, nbatch))
# state = random_state(sa, nbatch)
anneal_time = 10000
# begin_temp = 20.0
end_temp = 2.0
energy_gradient = CUDA.ones(anneal_time)
@btime track_equilibration!(HeatBath(), $sa, $state, energy_gradient, fill(end_temp, anneal_time))

cpu_state = Array(state)
state_10 = [sum(cpu_state[:,i] .* [2^(j-1) for j in 1:sa.n*sa.m]) for i in 1:size(cpu_state)[2]]
figure_count = [count(==(i), state_10) for i in 0:2^(sa.n*sa.m)-1]
distribution_count = 1.0 .* figure_count ./ sum(figure_count)
state_energy = [calculate_energy(sa, cpu_state, energy_gradient, i) for i in 1:nbatch]
state_energy_average = 1.0 * sum(state_energy) / nbatch
