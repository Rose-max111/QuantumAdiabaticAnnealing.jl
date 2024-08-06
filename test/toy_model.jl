using QuantumAdiabaticAnnealing, Test
using SparseArrays
using KrylovKit
using QuantumAdiabaticAnnealing: rule110, toy_model_state_energy, toy_model_transition_matrix

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

@testset "toy_model_state_energy" begin
    @test toy_model_state_energy(0b10100011, 5) == 1
    @test toy_model_state_energy(0b1011011001110111, 9) == 5
    @test toy_model_state_energy(0b1101101001110111, 9) == 3
    @test toy_model_state_energy(0b1111, 2; period_condition = true) == 2
    @test toy_model_state_energy(0b10111001, 4; period_condition = true) == 2
end

@testset "toy_model_transition_matrix" begin
    n = 4
    Temp = 0.001
    P = toy_model_transition_matrix(n, Temp)
    eigvals, eigvecs, infos = eigsolve(P, rand(Float64, 2^(2*n-2)), 2, :LM; maxiter = 5000)
    total = 0
    for i in 1:2^(2*n-2)
        if abs(Real(eigvecs[1][i])) > 1e-4
            println(reverse(bitstring(i-1)[end-n+1:end]), " ", bitstring(i-1)[end-(2*n-2)+1:end-n], " ", Real(eigvecs[1][i]))
            total += 1
            @test toy_model_state_energy(i-1, n) == 0
        end
    end
    @test total == 2^n

    n = 3
    Temp = 0.001
    P = toy_model_transition_matrix(n, Temp; period_condition = true)
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
