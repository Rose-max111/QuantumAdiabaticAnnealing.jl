using Test, QuantumAdiabaticAnnealing
using QuantumAdiabaticAnnealing: calculate_state_energy, spectral_gap, transition_matrix, plot_spectral_gap, plot_spectrum

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