using QuantumAdiabaticAnnealing, Test

@testset "count bits" begin
    @test QuantumAdiabaticAnnealing.countbits(3) == 2
    @test QuantumAdiabaticAnnealing.countbits(4) == 3
end

