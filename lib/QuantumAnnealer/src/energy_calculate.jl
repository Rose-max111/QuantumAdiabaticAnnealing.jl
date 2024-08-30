function van_der_waals_potential(u, v; C=862690)
    d2 = (u[1] - v[1])^2 + (u[2] - v[2])^2
    return 2*π*C / d2^3
end

function state_energy_calculation(state, Δ, T_set, nodes)
    energy = 0.0
    for i in eachindex(state)
        for j in (i+1):length(state)
            if state[i] == state[j] == 1
                energy += van_der_waals_potential(nodes[i], nodes[j])
            end
        end
    end

    for i in eachindex(state)
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
        mat_h = subspace === nothing ? mat(h) : mat(h, subspace)
        eigvals, eigvecs, info = eigsolve(mat_h, rand(Float64, size(mat_h, 1)), howmany, :SR; tol = 1e-7, maxiter = 5000)
        append!(val, [eigvals])
        if outputwhich === nothing
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
function get_low_energy_state(Δ, Ω, nodes; outputlevel = 0)
    N = length(nodes)
    sites = siteinds("S=1/2", N)
    
    os = OpSum()
    for j=1:N
        os += Ω[j], "Sx", j
        os += -Δ[j], "Proj0", j
        for k in (j+1):N
            os += van_der_waals_potential(nodes[j], nodes[k]), "Proj0", j, "Proj0", k
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

function get_low_energy_state_gpu(Δ, Ω, nodes; outputlevel = 0)
    N = length(nodes)
    sites = siteinds("S=1/2", N)
    
    os = OpSum()
    for j=1:N
        os += Ω[j], "Sx", j
        os += -Δ[j], "Proj0", j
        for k in (j+1):N
            os += van_der_waals_potential(nodes[j], nodes[k]), "Proj0", j, "Proj0", k
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

