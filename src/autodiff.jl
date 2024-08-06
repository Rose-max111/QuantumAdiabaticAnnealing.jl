function spinglass_hamiltonian(spin_M, Vtrans, gradient, n, m, t, Tmax)
    sp = spinglass_mapping(n, m;gradient=gradient)
    # @info "length(sp.onsite) = $(length(sp.onsite))"
    ret = 0.0
    # calculate spinglass Hz component
    for i in 1:length(sp.onsite)
        ret += 1.0 * spin_M[3*(i-1) + 3] * sp.onsite[i] * (t/Tmax)
    end
    
    # calculate the Hx component
    for i in 1:length(sp.onsite)
        ret += 1.0 * spin_M[3*(i-1) + 1] * Vtrans[i] * (1.0-t/Tmax)
    end

    # calculate the interaction term
    for edge in sp.edges
        u, v, w = edge[1], edge[2], edge[3]
        ret += 1.0 * w * spin_M[3*(u-1) + 3] * spin_M[3*(v-1) + 3] * (t/Tmax)
    end
    return ret
end

function spinglass_hamiltonian(sp, t, Tmax, Vtrans)
    spin_M = Vector{Float64}()
    for i in 1:length(sp.onsite)
        append!(spin_M, [sp.M[i][1], sp.M[i][2], sp.M[i][3]])
    end
    return spinglass_hamiltonian(spin_M, Vtrans, sp.gradient, sp.n, sp.m, t, Tmax)
end


function instantaneous_field_autodiff(sp, t, T, Vtrans::Vector{Float64})
    spin_M = Vector{Float64}()
    for i in 1:length(sp.onsite)
        append!(spin_M, [sp.M[i][1], sp.M[i][2], sp.M[i][3]])
    end
    # @info "length of spin_M is $(length(spin_M))"
    spin_M_grad = fill(0.0, 3*length(sp.onsite))
    # @info "instantaneous H is $(spinglass_hamiltonian(spin_M, Vtrans, sp.gradient, sp.n, sp.m, t, T))"
    autodiff(Enzyme.Reverse, 
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
    return -ret
end