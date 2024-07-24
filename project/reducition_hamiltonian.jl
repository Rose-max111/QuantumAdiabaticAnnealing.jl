using QuantumAdiabaticAnnealing
using SparseArrays
using KrylovKit
using Random
using CurveFit
using CairoMakie

function rule110(p, q, r)
    return (q + r + q*r + p*q*r) % 2
end

function calculate_state_energy(state, n; on_site_energy = nothing)
    ret = 0
    if on_site_energy == nothing
        on_site_energy = fill(0, n)
    end
    for l in 1:n
        ret += on_site_energy[l] * ((state & (2^(l-1))) > 0 ? 1 : -1)
    end
    for l in 1:n-2
        p = (state & (2^(l-1))) > 0 ? 1 : 0
        q = (state & (2^(l))) > 0 ? 1 : 0
        r = (state & (2^(l+1))) > 0 ? 1 : 0
        next = (state & (2^(l + n - 1))) > 0 ? 1 : 0
        ret += (rule110(p, q, r) ⊻ next)
    end
    return ret
end

function transition_matrix(n, Temp; on_site_energy = nothing)
    total_atoms = 2 * n - 2
    if on_site_energy == nothing
        on_site_energy = fill(0, n)
    end
    state_energy = [calculate_state_energy(i, n; on_site_energy) for i in 0:(2^total_atoms - 1)]


    row = Vector{Int}()
    col = Vector{Int}()
    val = Vector{Float64}()
    @info "finish calculating state energy"
    # P = spzeros(2^total_atoms, 2^total_atoms)
    for i in 0:(2^total_atoms - 1)
        other_prob = 0
        for j in 1:total_atoms
            reverse_mask = i ⊻ (2^(j-1))
            push!(row, reverse_mask + 1)
            push!(col, i+1)
            if state_energy[reverse_mask + 1] < state_energy[i + 1] # accept transition
                push!(val, 1.0 / total_atoms)
                # P[reverse_mask + 1, i + 1] = 1.0 / total_atoms
            else
                push!(val, exp(- (state_energy[reverse_mask + 1] - state_energy[i + 1]) / Temp) * 1.0 / total_atoms)
                # P[reverse_mask + 1, i + 1] = exp(- (state_energy[reverse_mask + 1] - state_energy[i + 1]) / Temp)
            end
            other_prob += val[end]
        end
        push!(row, i+1)
        push!(col, i+1)
        push!(val ,1 - other_prob)
    end
    # @info val
    return sparse(row, col, val)
end

mutable struct SimulatedAnnealingHamiltonian
    const b::Float64 # violate hamiltonian energy
    const n::Int # number of atoms per layer
    const m::Int # number of full layers (don't count the last layer)

    best_state::Vector{Bool}
    now_state::Vector{Bool}
    now_energy::Float64
    best_energy::Float64
    function SimulatedAnnealingHamiltonian(n::Int, m::Int;b::Float64 = 1)
        # internal states of the annealing process
        best_state = zeros(Bool, n * m + n)
        now_state = rand(Bool, n * m + n)
        new(b, n, m, best_state, now_state, 0.0, Inf)
    end
end

function evaluate(sa::SimulatedAnnealingHamiltonian, idp::Vector{Int}) # 1, 2, 3, 4->left, mid, right, next
    total_atoms = sa.n * sa.m + sa.n
    for i in 1:4
        if idp[i] < 1 || idp[i] > total_atoms
            return 0
        end
    end
    layers = [Int(floor((idp[i] - 1) / sa.n)) for i in 1:4]
    if layers[1]!=layers[2] || layers[2]!=layers[3] || layers[3]!=layers[4] - 1
        return 0
    end
    trueoutput = rule110(sa.now_state[idp[1]], sa.now_state[idp[2]], sa.now_state[idp[3]])
    return trueoutput ⊻ sa.now_state[idp[4]]
end

function calculate_energy(sa::SimulatedAnnealingHamiltonian)
    ret = 0
    for i in 1:sa.n*sa.m
        ret += evaluate(sa, [i, i + 1, i + 2, i + sa.n + 1])
    end
    return ret    
end

function step!(sa::SimulatedAnnealingHamiltonian, T::Float64, node::Int)
    minus = 0
    plus = 0
    if node == sa.n*sa.m+1 || node == sa.n * (sa.m+1)
        return 
    end
    if node > sa.n*sa.m # last layer
        minus += evaluate(sa, [node - sa.n - 1, node - sa.n, node - sa.n + 1, node])
        sa.now_state[node] = sa.now_state[node] ⊻ 1
        plus += evaluate(sa, [node - sa.n - 1, node - sa.n, node - sa.n + 1, node])
    else
        for i in 0:2
            minus += evaluate(sa, [node - 2 + i, node - 1 + i, node + i, node - 1 + i + sa.n])
        end
        minus += evaluate(sa, [node - sa.n - 1, node - sa.n, node - sa.n + 1, node])
        
        sa.now_state[node] = sa.now_state[node] ⊻ 1
        
        for i in 0:2
            plus += evaluate(sa, [node - 2 + i, node - 1 + i, node + i, node - 1 + i + sa.n])
        end
        plus += evaluate(sa, [node - sa.n - 1, node - sa.n, node - sa.n + 1, node])
    end
    energy_change = plus - minus

    p_acc = 0
    if energy_change < 0 
        p_acc = 1
    else
        p_acc = exp(-energy_change / T)
    end
    if rand() < p_acc
        sa.now_energy += energy_change
    else
        sa.now_state[node] = sa.now_state[node] ⊻ 1
    end

    # if sa.now_energy != calculate_energy(sa)
    #     @info "state = $(sa.now_state), node = $node, now_energy = $(sa.now_energy), calculate_energy = $(calculate_energy(sa)), sa.n = $(sa.n), sa.m = $(sa.m)"
    #     @info "plus = $plus, minus = $minus"
    #     for i in 0:2
    #         @info evaluate(sa, [node - 2 + i, node - 1 + i, node + i, node - 1 + i + sa.n])
    #     end
    #     @info evaluate(sa, [node - sa.n - 1, node - sa.n, node - sa.n + 1, node])
    #     error("error")
    # end
end

function track_equilibration!(sa::SimulatedAnnealingHamiltonian, tempscale = 4 .- (1:100 .-1) * 0.04, niters = 4000)
    nodelists = [i for i in 1:(sa.n*sa.m+sa.n)]
    exit = false
    for T in tempscale
        for runtime in 1:niters
            Random.shuffle!(nodelists)
            for node in nodelists
                step!(sa, T, node)
                if sa.now_energy < sa.best_energy
                    sa.best_energy = sa.now_energy
                    sa.best_state = copy(sa.now_state)
                end
                # if sa.now_energy == 0
                #     exit = true
                #     break
                # end
            end
            if exit == true
                break
            end
        end
        if exit == true
            break
        end
    end
    return sa
end

function annealing(n, m, run_annealing_iters, tempscale, niters)
    obserable = 0
    SA = SimulatedAnnealingHamiltonian(n, m; b = 1.0)
    SA.now_energy = calculate_energy(SA)
    for i in 1:run_annealing_iters
        SA = SimulatedAnnealingHamiltonian(n, m; b = 1.0)
        SA.now_energy = calculate_energy(SA)
        track_equilibration!(SA, tempscale, niters)
        # @info "runtime = $i, this_time_best_obj = $(SA.best_obj), success_time = $(success_time)"
        obserable += SA.now_energy
    end
    return 1.0 * obserable / run_annealing_iters
end

# n = 3
# m = 1
# obserable_average = []
# high_temp = 1
# for run_time = 1:1:20
#     obs = annealing(n, m, 20000, fill(1.0 * high_temp, run_time), 1)
#     push!(obserable_average,obs)
#     @info "runtime = $run_time, obserable_average = $obs"
# end


Temp = 1
ED_average_slope = []
sample_w = 3:10
f = Figure()
ax = Axis(f[1,1], title = "log(Energy Gap) vs log(Width)", xlabel = "log(Width)", ylabel = "log(Energy Gap)")
for average_field in 1:1:8
    ED_gap = []
    for n in sample_w
        on_site_energy = nothing
        if average_field != nothing
            on_site_energy = fill(average_field, n)
        end
        # @info on_site_energy
        total_atoms = 2*n-2
        P = transition_matrix(n, Temp; on_site_energy)
        eigvals, eigvecs, infos = eigsolve(P, rand(Float64, 2^(total_atoms)), 2, :LM; maxiter = 5000)
        push!(ED_gap, eigvals[1] - eigvals[2])
        # @info "n = $n, average_field = $(average_field), gap = $(eigvals[1] - eigvals[2]), converged = $(infos.converged), numiter = $(infos.numiter)"
    end
    x = [i for i in sample_w]
    logx = log.(x)
    logy = log.(Real.(ED_gap))
    scatter!(ax, logx, logy, label = "average_field = $average_field")
    # fit = curve_fit(LinearFit, logx, logy)
    # push!(ED_average_slope, fit.coefs[2])
    @info "average_field = $(average_field)"
end

axislegend(position = :rb)
f
# end

# x = [i for i in 0.1:0.1:5]
# ED_average_slope = Float64.(ED_average_slope)
# ax = Axis(f[1,1], title = "Slope vs Average Field", xlabel = "Average Field", ylabel = "Slope")
# scatter!(ax, x, ED_average_slope)
# f
# save("slope_vs_average_field.png", f)
# for val in ED_average_slope
#     println(val)
# end