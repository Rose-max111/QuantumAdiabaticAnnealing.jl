using logic_cooling, Test

@testset "count bits" begin
    @test logic_cooling.count_bits(53) == count_ones(53)
end