using QuantumAdiabaticAnnealing: step!, track_equilibration!, SimulatedAnnealingHamiltonian, random_state, Metropolis, HeatBath
using QuantumAdiabaticAnnealing: calculate_energy
using CUDA
using BenchmarkTools

sa = SimulatedAnnealingHamiltonian(40, 100)
nbatch = 10000

state = CuArray(random_state(sa, nbatch))
@btime CUDA.@sync step!(HeatBath(), $sa, $state, 1.0, 1)

cpu_state = random_state(sa, nbatch)
@btime step!(HeatBath(), $sa, $cpu_state, 1.0, 1)
#@test track_equilibration!(Metropolis(), sa, state) isa SimulatedAnnealingHamiltonian


sa = SimulatedAnnealingHamiltonian(5, 2)
nbatch = 10000

state = CuArray(random_state(sa, nbatch))
anneal_time = 100000
begin_temp = 10.0
end_temp = 2.0
anneal_temp = begin_temp .- (0:anneal_time - 1) .* (begin_temp - end_temp) / anneal_time
@btime CUDA.@sync track_equilibration!(HeatBath(), $sa, $state, anneal_temp)

cpu_state = Array(state)
state_energy = [calculate_energy(sa, cpu_state, i) for i in 1:nbatch]
state_energy_average = 1.0 * sum(state_energy) / nbatch


using KrylovKit
n = 5
Temp = end_temp
P = toy_model_transition_matrix(HeatBath(), n, Temp; period_condition = true)
eigvals, eigvecs, infos = eigsolve(P, rand(Float64, size(P)[1]), 2, :LR; maxiter = 5000)
ED_state_energy = [toy_model_state_energy(i, n; period_condition = true) for i in 0:2^(2n)-1]
maxeigvec_energy = sum(ED_state_energy .* (abs.(Real.(eigvecs[1])))) / abs(sum(Real.(eigvecs[1])))