function state_energy_calculation(state, Δ, T_set, nodes)
    energy = 0.0
    for i in 1:length(state)
        for j in (i+1):length(state)
            if state[i] == state[j] && state[i] == 1
                energy += 2*π*862690 / distance(nodes[i], nodes[j])^6
            end
        end
    end

    for i in 1:length(state)
        if state[i] == 1
            energy -= Δ[i](T_set)
        end
    end
    return energy
end

function pluse_energy_plot(hamiltonian, T_max, T_step, howmany)
    clocks = 0.0:T_step:T_max
    val = []
    for t in clocks
        h = hamiltonian |> attime(t)
        mat_h = mat(h)
        eigvals, eigvecs, info = eigsolve(mat_h, howmany, :SR)
        append!(val, [eigvals])
    end
    return val, clocks
end