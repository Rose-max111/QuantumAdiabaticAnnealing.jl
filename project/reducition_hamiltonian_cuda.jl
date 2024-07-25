using QuantumAdiabaticAnnealing: step!, track_equilibration!, SimulatedAnnealingHamiltonian, random_state, Metropolis, HeatBath
using CUDA
using BenchmarkTools

sa = SimulatedAnnealingHamiltonian(40, 100)
nbatch = 10000

state = CuArray(random_state(sa, nbatch))
@btime CUDA.@sync step!(HeatBath(), $sa, $state, 1.0, 1)

cpu_state = random_state(sa, nbatch)
@btime step!(HeatBath(), $sa, $cpu_state, 1.0, 1)
#@test track_equilibration!(Metropolis(), sa, state) isa SimulatedAnnealingHamiltonian