function distance(u::Tuple, v::Tuple)
    return sqrt((u[1] - v[1])^2 + (u[2] - v[2])^2)
end

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

function Hamiltonian_energy_plot(hamiltonian, T_max, T_step, howmany)
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

function pulse_energy_plot(Δ_T, Ω_T, nodes, T_max, T_step)
    clocks = 0.0:T_step:T_max
    val = []
    for t in clocks
        Δ = [Δ_T[i](t) for i in 1:length(nodes)]
        Ω = fill(Ω_T(t), length(nodes))
        eigvals = get_low_energy_state(Δ, Ω, nodes)
        append!(val, [eigvals])
    end
    return val, clocks
end

"""
    H / ħ = ∑ 1 / 2 * Ω_j * (|0⟩⟨1|_j + |1⟩⟨0|_j) - ∑ Δ_j * (|0⟩⟨0|_j) + ∑ 2π * 862690 / (x_j - x_k)^6 * (|0⟩⟨0|_j * |0⟩⟨0|_k)
"""
function get_low_energy_state(Δ, Ω, nodes)
    N = length(nodes)
    sites = siteinds("S=1/2", N)
    
    os = OpSum()
    for j=1:N
        os += Ω[j], "S+", j
        os += Ω[j], "S-", j
        os += -Δ[j], "Proj0", j
        for k in (j+1):N
            os += 2π * 862690 / (distance(nodes[j], nodes[k])^6), "Proj0", j, "Proj0", k
        end
    end
    H = MPO(os, sites)

    h = 10.0
    weight = 1000
    nsweeps = 40
    maxdim = [10,10,10,20,20,40,80,100,200,200]
    cutoff = [1E-8]
    noise = [1E-6]

    #
    # Compute the ground state psi0
    #

    psi0_init = random_mps(sites;linkdims = 4)
    energy0, psi0 = dmrg(H, psi0_init; nsweeps, maxdim, cutoff, noise)

    psi1_init = random_mps(sites;linkdims = 4)
    energy1, psi1 = dmrg(H, [psi0], psi1_init; nsweeps, maxdim, cutoff, noise, weight)

    psi2_init = random_mps(sites;linkdims = 8)
    energy2, psi2 = dmrg(H, [psi0, psi1], psi2_init; nsweeps, maxdim, cutoff, noise, weight)

    real_energy1 = inner(psi1', H, psi1)
    real_energy2 = inner(psi2', H, psi2)
    return [energy0, real_energy1, real_energy2]
end