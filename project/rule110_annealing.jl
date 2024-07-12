using QuantumAdiabaticAnnealing: rule110_transverse_generate
using QuantumAdiabaticAnnealing: Hamiltonian_energy_plot
using BloqadeExpr
using Bloqade, CairoMakie

locations, weights, inputs_id, outputs_id, input_layer_id = rule110_transverse_generate(1, 1)
locations = map(t -> (Float64(t[1]), Float64(t[2])), locations)

atoms = Bloqade.AtomList(locations)

for i in inputs_id
    weights[i] = -50
end
for i in input_layer_id
    weights[i[1]] = weights[i[2]] = -50
end

# for i in outputs_id
#     weights[i] = -50
# end

T_max = 4
T_pure = 0.8
Δ_min= 0.2 * 2π
Δ_max = 5 * 2π
# Δ = map(1:length(locations)) do idx
#     piecewise_linear(clocks = [0, T_pure, T_max - T_pure, T_max], values = [Δ_min * weights[idx], Δ_min * weights[idx], Δ_max * weights[idx], Δ_max * weights[idx]])
# end
Δ = map(1:length(locations)) do idx
    piecewise_linear(clocks = [0, T_max], values = [Δ_min * weights[idx], Δ_max * weights[idx]])
end


Ω_min = 0 
Ω_max = 1 * 2π 
# Ω = piecewise_linear(clocks = [0, T_pure, T_max - T_pure, T_max], values = [Ω_min,  Ω_max, Ω_max, Ω_min])
Ω = piecewise_linear(clocks = [0, T_max], values = [Ω_max, Ω_max])

hamiltonian = rydberg_h(atoms, Ω = Ω, Δ = Δ)
space = blockade_subspace(atoms, 2.05)
eigvals, times = Hamiltonian_energy_plot(hamiltonian, T_max, 0.02, 3; subspace = space)

delta_over_omega = [Ω(t) / Δ[1](t) for t in times]
fig = Figure(resolution = (1500, 800))
ax = Axis(fig[1, 1], title = "Energy Line", xlabel = "Ω  /  Δ", ylabel = "Energy gap δ")

lines!(ax, delta_over_omega, [(val[2] - val[1]) for val in eigvals], color = :blue, linewidth = 2, label = "ED FE vs GS")
# lines!(ax, delta_over_omega, [(val[3] - val[1]) for val in eigvals], color = :red, linewidth = 2, label = "ED SE vs GS")
# lines!(ax, delta_over_omega, [(val[4] - val[1]) for val in eigvals], color = :brown, linewidth = 2, label = "ED TE vs GS")
axislegend(ax; position = :lb, labelsize = 15)


fig