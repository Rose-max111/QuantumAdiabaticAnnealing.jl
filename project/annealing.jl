using QuantumAdiabaticAnnealing: generate_some_graph, state_energy_calculation, pulse_energy_plot, Hamiltonian_energy_plot
using Bloqade, CairoMakie

new_graph_nodes, new_graph_weights = generate_some_graph()

atoms = AtomList(new_graph_nodes)

radius = 4.8

detuning_max = 400
detuning_min = - detuning_max

T_max = 10.0
# empty_time = 2.0
# Δ = map(1:length(new_graph_nodes)) do idx
#     piecewise_linear(clocks = [0.0, empty_time, T_max - empty_time, T_max], values = [detuning_min * new_graph_weights[idx], detuning_min * new_graph_weights[idx], detuning_max * new_graph_weights[idx], detuning_max * new_graph_weights[idx]])
# end

# Ω_max = 2π * 4
# Ω = piecewise_linear(clocks = [0.0, empty_time, T_max - empty_time, T_max], values = [0, Ω_max, Ω_max, 0])

Δ = map(1:length(new_graph_nodes)) do idx
    piecewise_linear(clocks = [0, T_max], values = [0.05 * new_graph_weights[idx], detuning_max * new_graph_weights[idx]])
end

Ω = piecewise_linear(clocks = [0, T_max], values = [1, 1])

hamiltonian = rydberg_h(atoms, Ω = Ω, Δ = Δ)
subspace = blockade_subspace(atoms, radius)

eigvals, times = Hamiltonian_energy_plot(hamiltonian, T_max, 0.1, 3)

DMRGeigvals, times = pulse_energy_plot(Δ, Ω, new_graph_nodes, T_max, 0.1)


delta_over_omega = [Δ[1](t) / Ω(t) for t in times]
fig = Figure(resolution = (1500, 800))
ax = Axis(fig[1, 1], title = "Energy Line", xlabel = "Δ  /  Ω", ylabel = "log(Energy minus ground state energy)")

# lines!(ax, delta_over_omega, [(val[1]) for val in eigvals], color = :blue, linewidth = 2, label = "first excitation vs ground state")
# lines!(ax, delta_over_omega, [(val[3]) for val in eigvals], color = :green, linewidth = 2, label = "second excitation vs ground state")
lines!(ax, delta_over_omega, [(val[1]) for val in DMRGeigvals], color = :brown, linewidth = 2, label = "DMRG FE vs GS")
lines!(ax, delta_over_omega, [(val[3]) for val in DMRGeigvals], color = :red, linewidth = 2, label = "DMRG SE vs GS")
axislegend(ax; position = :lb, labelsize = 15)


fig