using QuantumAdiabaticAnnealing
using Enzyme

sp=spinglass_mapping(3,4)
spin_M = Vector{Float64}()
for i in 1:length(sp.onsite)
    append!(spin_M, [sp.M[i][1], sp.M[i][2], sp.M[i][3]])
end
# @info "length of spin_M is $(length(spin_M))"
spin_M_grad = fill(0.0, 3*length(sp.onsite))
Vtrans=fill(1.0, length(sp.onsite))
t=4.0
T=10.0
@info "instantaneous H is $(spinglass_hamiltonian(spin_M, Vtrans, sp.gradient, sp.n, sp.m, t, T))"
autodiff(Reverse, 
        spinglass_hamiltonian,
        Duplicated(spin_M, spin_M_grad),
        Const(Vtrans),
        Const(sp.gradient),
        Const(sp.n),
        Const(sp.m),
        Const(t),
        Const(T))
ret = Vector{Vector{Float64}}()
for i in 1:length(sp.onsite)
    push!(ret, [spin_M_grad[3*(i-1)+1], spin_M_grad[3*(i-1)+2], spin_M_grad[3*(i-1)+3]])
end