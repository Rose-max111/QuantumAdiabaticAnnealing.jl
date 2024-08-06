using QuantumAdiabaticAnnealing, Test

@testset "Point" begin
    p1 = Point(1.0, 2.0)
    p2 = Point(3.0, 4.0)
    @test zero(p1) ≈ Point(0.0, 0.0)
    @test p1 + p2 ≈ Point(4.0, 6.0)
    @test p1 - p2 ≈ Point(-2.0, -2.0)
    @test p1 * 2 ≈ Point(2.0, 4.0)
    @test p1 / 2 ≈ Point(0.5, 1.0)
    @test distance(p1, p2) ≈ 2.8284271247461903
end

