# e.g. rule 110 could be specified as `CellularAutomata1D{110}()`
struct CellularAutomata1D{INT} end

function logic_gate(::CellularAutomata1D{N}, p, q, r) where N
    return (N >> (p << 2 | q << 1 | r)) & 1
end

function generic_logic_grate(p, q, r, N)
    return (N >> (p << 2 | q << 1 | r)) & 1
end
rule110(p, q, r) = logic_gate(CellularAutomata1D{110}(), p, q, r)
rule000(p, q, r) = logic_gate(CellularAutomata1D{0}(), p, q, r)
automatarule(p, q, r, id) = generic_logic_grate(p, q, r, id)