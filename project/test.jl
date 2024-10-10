using QuantumAdiabaticAnnealing
using QuantumAdiabaticAnnealing:spectral
using QuantumAdiabaticAnnealing:midposition_calculate

1.0 - midposition_calculate(Exponentialtype(), 10.0, 1.0, 1.5)

ans = []
for n in 3:5
    Mat = toy_model_transition_matrix(CellularAutomata1D{110}(), HeatBath(), n, 4, 10; on_site_energy=nothing, period_condition=true, energy_gradient=1.0)
    @info "finish calculating transition matrix"
    eigvals, eigvecs = spectral(Mat, 2)
    @info "n = $n, Δ = $(eigvals[1] - eigvals[2])"
    push!(ans, eigvals[1] - eigvals[2])
end
b = ans .* Vector(3:5)