using QuantumAdiabaticAnnealing, Test
using SparseArrays
using KrylovKit
using QuantumAdiabaticAnnealing: rule110, toy_model_state_energy, toy_model_transition_matrix, calculate_state_energy, spectral_gap, transition_matrix, plot_spectral_gap, plot_spectrum


@testset "rule110" begin
    @test rule110(0, 0, 0) == 0
    @test rule110(0, 0, 1) == 1
    @test rule110(0, 1, 0) == 1
    @test rule110(0, 1, 1) == 1
    @test rule110(1, 0, 0) == 0
    @test rule110(1, 0, 1) == 1
    @test rule110(1, 1, 0) == 1
    @test rule110(1, 1, 1) == 0
end

@testset "toy_model_transition_matrix" begin
    n, m = 4, 2
    Temp = 0.1
    ca = CellularAutomata1D{110}()
    rule = HeatBath()
    P = toy_model_transition_matrix(ca, rule, n, m, Temp; period_condition = false, on_site_energy=nothing, energy_gradient=1)
    eigvals, eigvecs, infos = eigsolve(P, rand(Float64, 2^(2*n-2)), 2, :LM; maxiter = 5000)
    total = 0
    for i in 1:2^(2*n-2)
        if abs(Real(eigvecs[1][i])) > 1e-4
            println(reverse(bitstring(i-1)[end-n+1:end]), " ", bitstring(i-1)[end-(2*n-2)+1:end-n], " ", Real(eigvecs[1][i]))
            total += 1
            @test toy_model_state_energy(i-1, n; period_condition=false, on_site_energy=nothing, energy_gradient=1) == 0
        end
    end
    @test total == 2^n

    n = 3
    Temp = 0.1
    P = toy_model_transition_matrix(n, Temp; period_condition = true, on_site_energy=nothing, energy_gradient=1)
    eigvals, eigvecs, infos = eigsolve(P, rand(Float64, 2^(2*n)), 2, :LM; maxiter = 5000)
    total = 0
    for i in 1:2^(2*n)
        if abs(Real(eigvecs[1][i])) > 1e-4
            println(reverse(bitstring(i-1)[end-n+1:end]), " ", bitstring(i-1)[end-(2*n)+1:end-n], " ", Real(eigvecs[1][i]))
            total += 1
            @test toy_model_state_energy(i-1, n; period_condition = true) == 0
        end
    end
    @test total == 2^n
end

@testset "state energy" begin
    @test calculate_state_energy(0, fill(0.5, 4)) == -0.5 * 4 + 2
end

@testset "spectral gap" begin
    @test isapprox(spectral_gap(Metropolis(), 1e10, fill(0.0, 3)), 0; atol=1e-5)
    @test spectral_gap(Metropolis(), 1.0, fill(0.0, 3)) ≈ 0.21850742086556896
    @test spectral_gap(HeatBath(), 1e10, fill(0.0, 3)) ≈ 0.24999999999458677
    @test spectral_gap(HeatBath(), 1.0, fill(0.0, 3)) ≈ 0.1573320844237781
end

@testset "transition matrix" begin
    P = transition_matrix(Metropolis(), 2.0, fill(0.5, 4))
    @test size(P) == (2^6, 2^6)
    @test all(≈(1), sum(P; dims=1))

    P = transition_matrix(HeatBath(), 2.0, fill(0.5, 4))
    @test size(P) == (2^6, 2^6)
    @test all(≈(1), sum(P; dims=1))
end

@testset "visualize" begin
    plot_spectral_gap(HeatBath(), 10)
    plot_spectral_gap(Metropolis(), 10)
    plot_spectrum(Metropolis(), 0.3, fill(1.0, 3), 6)
    plot_spectrum(HeatBath(), 0.3, fill(1.0, 4), 16)
end