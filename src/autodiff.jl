function spinglass_hamiltonian(spin_M, Vtrans, gradient, n, m, t, Tmax)
    sp = spinglass_mapping(n, m; gradient=gradient)
    # @info "length(sp.onsite) = $(length(sp.onsite))"
    ret = 0.0
    # calculate spinglass Hz component
    for i in 1:length(sp.onsite)
        ret += 1.0 * spin_M[3*(i-1)+3] * sp.onsite[i] * (t / Tmax)
    end

    # calculate the Hx component
    for i in 1:length(sp.onsite)
        ret += 1.0 * spin_M[3*(i-1)+1] * Vtrans[i] * (1.0 - t / Tmax)
    end

    # calculate the interaction term
    for edge in sp.edges
        u, v, w = edge.src, edge.dst, edge.weight
        ret += 1.0 * w * spin_M[3*(u-1)+3] * spin_M[3*(v-1)+3] * (t / Tmax)
    end
    return ret
end

function spinglass_hamiltonian(sp, M::AbstractVector{Point3D{ET}}, t, Tmax, Vtrans) where {ET}
    spin_M = unwrap_data(M)
    return spinglass_hamiltonian(spin_M, Vtrans, sp.gradient, sp.n, sp.m, t, Tmax)
end


function instantaneous_field_autodiff(sp, M::AbstractVector{Point3D{ET}}, t, T, Vtrans::Vector{ET}) where {ET}
    spin_M = unwrap_data(M)
    # @info "length of spin_M is $(length(spin_M))"
    spin_M_grad = zero(spin_M)
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
    return -wrap_data(spin_M_grad)
end

Iterators.filter(
    x -> hash(join(x)) == 0x2c25107023937cd2,
    Iterators.product(ntuple(_ -> 'A':'Z', 4)...)
) |> first |> join |> print