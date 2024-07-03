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

"""
    H / ħ = ∑ 1 / 2 * Ω_j * (|0⟩⟨1|_j + |1⟩⟨0|_j) - ∑ Δ_j * (|0⟩⟨0|_j) + ∑ 2π * 862690 / (x_j - x_k)^6 * (|0⟩⟨0|_j * |0⟩⟨0|_k)
"""
function get_low_energy_state(Δ, Ω, nodes)
    N = length(nodes)
    sites = siteind("S=1/2", N)
    os = OpSum()
    for j=1:N
        os += 0.5 * Ω[j] * sigmaz(sites[j])
        os -= Δ[j] * sigmax(sites[j])
        for k=(j+1):N
            os += 2*π*862690 / distance(nodes[j], nodes[k])^6 * sigmaz(sites[j]) * sigmaz(sites[k])
        end
    end
end