using Test, QuantumAdiabaticAnnealing

@testset "rule110" begin
    function rule110(p, q, r)
        return (q + r + q*r + p*q*r) % 2
    end
    for p in 0:1, q in 0:1, r in 0:1
        @test rule110(p, q, r) == logic_gate(CellularAutomata1D{110}(), p, q, r)
    end

    @test QuantumAdiabaticAnnealing.simulate(CellularAutomata1D{110}(), Bool[0, 1, 0, 0], 2) == Bool[0 1 1; 1 1 1; 0 0 0; 0 0 1]
end
