using QuantumAdiabaticAnnealing: generate_some_graph, state_energy_calculation, pulse_energy_plot, Hamiltonian_energy_plot
using QuantumAdiabaticAnnealing: generate_random_lattice
using Bloqade, CairoMakie

# new_graph_nodes, new_graph_weights = generate_some_graph()
new_graph_nodes, new_graph_weights = generate_random_lattice(4, 4.5, 0.8)

atoms = AtomList(new_graph_nodes)

T_max = 10
detuning_min = 3 * 2π
detuning_max = 30 * 2π
Δ = map(1:length(new_graph_nodes)) do idx
    piecewise_linear(clocks = [0, T_max], values = [detuning_min * new_graph_weights[idx], detuning_max * new_graph_weights[idx]])
end

Ω = piecewise_linear(clocks = [0, T_max], values = [4 * 2π, 4 * 2π])

hamiltonian = rydberg_h(atoms, Ω = Ω, Δ = Δ)
# subspace = blockade_subspace(atoms, radius)

eigvals, times = Hamiltonian_energy_plot(hamiltonian, T_max, 0.01, 3)

DMRGeigvals, times = pulse_energy_plot(Δ, Ω, new_graph_nodes, T_max, 0.01)


tdelta_over_2pi = [Δ[1](t) / 2π for t in times]
fig = Figure(resolution = (1500, 800))
ax = Axis(fig[1, 1], title = "Energy Line", xlabel = "Δ  /  2π", ylabel = "Energy gap δ / 2π")

lines!(ax, delta_over_2pi, [(val[2] - val[1])/2π for val in eigvals], color = :blue, linewidth = 2, label = "ED FE vs GS")
# lines!(ax, delta_over_omega, [(val[3]) for val in eigvals], color = :green, linewidth = 2, label = "second excitation vs ground state")
lines!(ax, delta_over_2pi, [(val[2] - val[1])/2π for val in DMRGeigvals], color = :brown, linewidth = 2, label = "DMRG FE vs GS")
# lines!(ax, delta_over_omega, [(val[3]) for val in DMRGeigvals], color = :red, linewidth = 2, label = "DMRG SE vs GS")
axislegend(ax; position = :lb, labelsize = 15)


fig