using QuantumAdiabaticAnnealing, Test

@testset "count bits" begin
    @test QuantumAdiabaticAnnealing.countbits(3) == 2
    @test QuantumAdiabaticAnnealing.countbits(4) == 3
end

@testset "distance" begin
    @test QuantumAdiabaticAnnealing.distance((0,0), (3,4)) ≈ 5.0
    @test QuantumAdiabaticAnnealing.distance((-1,-2), (3,1)) ≈ 5.0
end