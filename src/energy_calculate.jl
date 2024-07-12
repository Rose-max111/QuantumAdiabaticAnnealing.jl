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

function Hamiltonian_energy_plot(hamiltonian, T_max, T_step, howmany; subspace = nothing, outputwhich = nothing)
    clocks = 0.0:T_step:T_max
    val = []
    for t in clocks
        h = hamiltonian |> attime(t)
        if subspace == nothing
            mat_h = mat(h)
        else
            mat_h = mat(h, subspace)
        end
        eigvals, eigvecs, info = eigsolve(mat_h, howmany, :SR; tol = 1e-3, maxiter = 5000)
        append!(val, [eigvals])
        if outputwhich == nothing
            @info t, eigvals[1], eigvals[2], eigvals[2] - eigvals[1], info.converged, info.numiter
        else
            @info t, eigvals[2] - eigvals[1], eigvals[3] - eigvals[2], eigvals[4] - eigvals[3], info.converged, info.numiter
        end

    end
    return val, clocks
end

function pulse_energy_plot(Δ_T, Ω_T, nodes, T_max, T_step)
    clocks = 0.0:T_step:T_max
    val = []
    last_psi = [nothing, nothing, nothing]
    for t in clocks
        Δ = [Δ_T[i](t) for i in 1:length(nodes)]
        Ω = fill(Ω_T(t), length(nodes))
        eigvals, new_psi = get_low_energy_state(Δ, Ω, nodes;)
        eigvals = sort(eigvals)
        push!(val, eigvals)
        @show t, eigvals
    end
    return val, clocks
end

"""
    H / ħ = ∑ 1 / 2 * Ω_j * (|0⟩⟨1|_j + |1⟩⟨0|_j) - ∑ Δ_j * (|0⟩⟨0|_j) + ∑ 2π * 862690 / (x_j - x_k)^6 * (|0⟩⟨0|_j * |0⟩⟨0|_k)
"""
function get_low_energy_state(Δ, Ω, nodes; given_psi = [nothing, nothing, nothing], outputlevel = 0)
    N = length(nodes)
    sites = siteinds("S=1/2", N)
    
    os = OpSum()
    for j=1:N
        os += Ω[j], "Sx", j
        os += -Δ[j], "Proj0", j
        for k in (j+1):N
            os += 2π * 862690 / (distance(nodes[j], nodes[k])^6), "Proj0", j, "Proj0", k
        end
    end
    H = MPO(os, sites)
    # H_gpu = mtl(H)

    weight = 2000
    nsweeps = 50
    maxdim = [10,10,20,40,80,100,200,400,400,800]
    cutoff = [1E-10]
    noise = [1E-7, 1E-7, 1E-8, 1E-8, 1E-9, 1E-9, 1E-10, 1E-11, 1E-11, 0.0]


    psi_init = random_mps(sites;linkdims = 2)
    # psi_init_gpu = mtl(psi_init)
    energy0, psi0 = dmrg(H, psi_init; nsweeps, maxdim, cutoff, noise, outputlevel = outputlevel)

    psi_pre = Vector{typeof(psi0)}()
    energy = []
    howmany = 5
    append!(energy, energy0)
    push!(psi_pre, psi0)

    for T in 2:howmany
        psi_init = random_mps(sites;linkdims = 2)
        # psi_init_gpu = mtl(psi_init)
        energy0, psi0 = dmrg(H, psi_pre, psi_init; nsweeps, maxdim, cutoff, noise, weight, outputlevel = outputlevel)
        append!(energy, energy0)
        push!(psi_pre, psi0)
    end

    sort!(energy)
    return energy, psi_pre
end



function get_low_energy_state_gpu(Δ, Ω, nodes; given_psi = [nothing, nothing, nothing], outputlevel = 0)
    N = length(nodes)
    sites = siteinds("S=1/2", N)
    
    os = OpSum()
    for j=1:N
        os += Ω[j], "Sx", j
        os += -Δ[j], "Proj0", j
        for k in (j+1):N
            os += 2π * 862690 / (distance(nodes[j], nodes[k])^6), "Proj0", j, "Proj0", k
        end
    end
    H = MPO(os, sites)
    H_gpu = cu(H)

    weight = 2000
    nsweeps = 50
    maxdim = [10,10,20,40,80,100,200,400,400,800]
    cutoff = [1E-10]
    noise = [1E-7, 1E-7, 1E-8, 1E-8, 1E-9, 1E-9, 1E-10, 1E-11, 1E-11, 0.0]

    psi_init = random_mps(sites;linkdims = 4)
    psi_init_gpu = cu(psi_init)
    energy0, psi0 = dmrg(H_gpu, psi_init_gpu; nsweeps, maxdim, cutoff, noise, outputlevel = outputlevel)

    psi_pre = Vector{typeof(psi0)}()
    energy = []
    howmany = 5
    append!(energy, energy0)
    push!(psi_pre, psi0)

    for T in 2:howmany
        psi_init = random_mps(sites;linkdims = 4)
        psi_init_gpu = cu(psi_init)
        energy0, psi0 = dmrg(H_gpu, psi_pre, psi_init_gpu; nsweeps, maxdim, cutoff, noise, weight, outputlevel = outputlevel)
        append!(energy, energy0)
        push!(psi_pre, psi0)
    end

    sort!(energy)
    return energy, psi_pre
end