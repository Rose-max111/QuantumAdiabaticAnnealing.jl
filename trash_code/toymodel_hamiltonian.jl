using QuantumAdiabaticAnnealing
using SparseArrays
using KrylovKit
using Random
using CurveFit
using CairoMakie


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
sample_w = 3:9
sample_average_field = 0:1:0
f = Figure()
ax = Axis(f[1,1], title = "log(Energy Gap) vs log(Width) (period_condition)", xlabel = "log(Width)", ylabel = "log(Energy Gap)")
for average_field in sample_average_field
    ED_gap = []
    for n in sample_w
        on_site_energy = nothing
        if average_field != nothing
            on_site_energy = fill(average_field, n)
        end
        # @info on_site_energy
        total_atoms = 2*n
        P = toy_model_transition_matrix(n, Temp; on_site_energy, period_condition = true)
        eigvals_positive, eigvecs_positive, infos = eigsolve(P, rand(Float64, 2^(total_atoms)), 2, :LR; maxiter = 5000)
        eigvals_negetive, eigvecs_negetive, infos = eigsolve(P, rand(Float64, 2^(total_atoms)), 2, :SR; maxiter = 5000)
        second_largest_eigenval = max(abs(Real(eigvals_positive[2])), abs(Real(eigvals_negetive[1])))
        push!(ED_gap, Real(eigvals_positive[1]) - second_largest_eigenval)
        # @info "n = $n, average_field = $(average_field), gap = $(eigvals[1] - eigvals[2]), converged = $(infos.converged), numiter = $(infos.numiter)"
    end
    x = [i for i in sample_w]
    logx = log.(x)
    logy = log.(ED_gap)
    scatter!(ax, logx, logy, label = "average_field = $average_field")
    fit = curve_fit(LinearFit, logx, logy)
    push!(ED_average_slope, fit.coefs[2])
    @info "average_field = $(average_field)"
end

# axislegend(position = :lb)
f
# end

x = [i for i in sample_average_field]
# ax = Axis(f[1,1], title = "log(Energy Gap) vs log(Width) (period_condition)", xlabel = "log(Width)", ylabel = "log(Energy Gap)")
ED_average_slope = Float64.(ED_average_slope)
ax = Axis(f[1,1], title = "Slope vs Average Field", xlabel = "Average Field", ylabel = "Slope")
scatter!(ax, x, ED_average_slope)
f
# save("slope_vs_average_field.png", f)
# for val in ED_average_slope
#     println(val)
# end