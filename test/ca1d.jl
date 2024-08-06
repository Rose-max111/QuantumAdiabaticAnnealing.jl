using Test, QuantumAdiabaticAnnealing

@testset "rule110" begin
    function rule110(p, q, r)
        return (q + r + q*r + p*q*r) % 2
    end
    for p in 0:1, q in 0:1, r in 0:1
        @test rule110(p, q, r) == logic_gate(CellularAutomata1D{110}(), p, q, r)
    end
end
