using Test, QuantumAdiabaticAnnealing

@testset "rule110" begin
    @test count(!iszero, [rule110(x, y, z) for x=0:1, y=0:1, z=0:1]) == 5
end
