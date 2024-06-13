using logic_cooling, Test

@testset "count bits" begin
    @test logic_cooling.countbits(3) == 2
    @test logic_cooling.countbits(4) == 3
end