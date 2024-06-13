using logic_cooling: generate_some_graph, calculate_energy, generate_energy_plot
using Bloqade, CairoMakie

new_graph_nodes, new_graph_weights = generate_some_graph()

atoms = AtomList(new_graph_nodes)

radius = 4.8

detuning_max = 400
detuning_min = - detuning_max

T_max = 10.0
empty_time = 2.0
Δ = map(1:length(new_graph_nodes)) do idx
    piecewise_linear(clocks = [0.0, empty_time, T_max - empty_time, T_max], values = [detuning_min * new_graph_weights[idx], detuning_min * new_graph_weights[idx], detuning_max * new_graph_weights[idx], detuning_max * new_graph_weights[idx]])
end

Ω_max = 2π * 4
Ω = piecewise_linear(clocks = [0.0, empty_time, T_max - empty_time, T_max], values = [0, Ω_max, Ω_max, 0])

hamiltonian = rydberg_h(atoms, Ω = Ω, Δ = Δ)
subspace = blockade_subspace(atoms, radius)

eigvals, times = generate_energy_plot(hamiltonian, T_max, 0.1, 3, subspace)

PlotTime = [t for t in times]
fig = Figure(resolution = (1500, 800))
ax = Axis(fig[1, 1], title = "Energy Line", xlabel = "Time", ylabel = "log(Energy minus ground state energy)")

lines!(ax, PlotTime, [log(val[2] - val[1]) for val in eigvals], color = :blue, linewidth = 2, label = "first excitation")
lines!(ax, PlotTime, [log(val[3] - val[1]) for val in eigvals], color = :green, linewidth = 2, label = "second excitation")
# lines!(ax, PlotTime, [log(val[4] - val[1]) for val in eigvals], color = :brown, linewidth = 2, label = "third excitation")
# lines!(ax, PlotTime, [log(val[5] - val[1]) for val in eigvals], color = :red, linewidth = 2, label = "fourth excitation")
axislegend(ax; position = :lb, labelsize = 15)


fig
# nsites = length(new_graph_nodes)
# clocks = 0.0:1e-2:T_max
# prob2 = KrylovEvolution(zero_state(subspace), clocks, hamltonian)
# emulate!(prob2)

# bitstring_hist(prob2.reg; nlargest = 20)

# best_bit_strings = most_probable(prob2.reg, 10)
# desire = [0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0]
# Bloqade.plot(atoms, blockade_radius = radius; colors = [iszero(b) ? "white" : "red" for b ∈ desire])


# calculate_energy(best_bit_strings[9], Δ, T_max, new_graph_nodes)
# calculate_energy(desire, Δ, T_max, new_graph_nodes)