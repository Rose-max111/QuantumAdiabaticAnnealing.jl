using QuantumAdiabaticAnnealing
using SparseArrays
using KrylovKit
using Random
using CurveFit
using CairoMakie
using Yao.BitBasis

function rule110(p, q, r)
    return (q + r + q*r + p*q*r) % 2
end

using Test
@testset "rule110" begin
    @test count(!iszero, [rule110(x, y, z) for x=0:1, y=0:1, z=0:1]) == 5
end

# NOTE: we'd better use periodic boundary condition
function calculate_state_energy(state, on_site_energy::AbstractVector{T}) where T
    ret = zero(T)
    n = length(on_site_energy)
    for l in 1:n
        ret += on_site_energy[l] * (isone(readbit(state, l)) ? 1 : -1)
    end
    for l in 1:n-2
        p = readbit(state, l)
        q = readbit(state, l+1)
        r = readbit(state, l+2)
        next = readbit(state, l + n)
        ret += (rule110(p, q, r) == next)
    end
    return ret
end

@testset "state energy" begin
    @test calculate_state_energy(0, fill(0.5, 4)) == -0.5 * 4 + 2
end

abstract type TransitionRule end
struct HeatBath <: TransitionRule end
struct Metroplis <: TransitionRule end

function transition_matrix(rule::TransitionRule, Temp, on_site_energy::AbstractVector{T}) where {T}
    n = length(on_site_energy)
    total_atoms = 2 * n - 2
    state_energy = [calculate_state_energy(i, on_site_energy) for i in 0:(2^total_atoms - 1)]

    row = Vector{Int}()
    col = Vector{Int}()
    val = Vector{Float64}()
    @info "finish calculating state energy"
    # P = spzeros(2^total_atoms, 2^total_atoms)
    for i in 0:(2^total_atoms - 1)
        other_prob = 0
        for j in 1:total_atoms
            reverse_mask = i ⊻ bmask(Int, j)
            push!(row, reverse_mask + 1)  # the transition from i to reverse_mask
            push!(col, i+1)
            ΔE = state_energy[reverse_mask + 1] - state_energy[i + 1]
            push!(val, update(rule, Temp, ΔE, 1.0 / total_atoms))
            other_prob += val[end]
        end
        push!(row, i+1)
        push!(col, i+1)
        push!(val ,1 - other_prob)
    end
    # @info val
    return sparse(row, col, val)
end

# Metropolis-Hastings algorithm
function update(::Metroplis, Temp, ΔE, prior)
    if ΔE < 0
        1.0
    else
        exp(- (ΔE) / Temp)
    end * prior
end

function update(::HeatBath, Temp, ΔE, prior)
    exp(-ΔE / Temp) / (1 + exp(-ΔE / Temp)) * prior
end


@testset "transition matrix" begin
    P = transition_matrix(Metroplis(), 2.0, fill(0.5, 4))
    @test size(P) == (2^6, 2^6)
    @test all(≈(1), sum(P; dims=1))

    P = transition_matrix(HeatBath(), 2.0, fill(0.5, 4))
    @test size(P) == (2^6, 2^6)
    @test all(≈(1), sum(P; dims=1))
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

using LinearAlgebra
function spectral(rule::TransitionRule, Temp, on_site_energy, k::Int)
    n = length(on_site_energy)
    total_atoms = 2 * n - 2
    P = transition_matrix(rule, Temp, on_site_energy)
    eigvals, eigvecs, info = eigsolve(P, rand(Float64, 2^(total_atoms)), k, :LM; maxiter = 5000)
    @assert info.converged >= k "not converged"
    @assert eigvals[1] ≈ 1
    #@assert all(isapprox(ev, real(ev); atol=1e-5) for ev in eigvals) "complex eigvals: $eigvals"
    return real.(eigvals), eigvecs
end

function spectral_gap(rule::TransitionRule, Temp, on_site_energy)
    evals, _ = spectral(rule, Temp, on_site_energy, 2)
    return abs(abs(evals[1]) - abs(evals[2]))
end

@testset "spectral gap" begin
    @test isapprox(spectral_gap(Metroplis(), 1e10, fill(0.0, 3)), 0; atol=1e-5)
    @test spectral_gap(Metroplis(), 1.0, fill(0.0, 3)) ≈ 0.21850742086556896
    @test spectral_gap(HeatBath(), 1e10, fill(0.0, 3)) ≈ 0.24999999999458677
    @test spectral_gap(HeatBath(), 1.0, fill(0.0, 3)) ≈ 0.1573320844237781
end

# NOTE: maybe we should switch to using probability to get correct result?

function plot_spectral_gap(rule::TransitionRule, Temp;
            on_site_energies = 2 .^ (0:8),
            sample_w = 3:9,
        )
    f = Figure()
    ax = Axis(f[1,1], title = "Spectral Gap v.s. Width", xscale=log10, yscale=log10, xlabel = "Width", ylabel = "Spectral Gap")
    for average_field in on_site_energies
        ED_gap = []
        for n in sample_w
            gap = spectral_gap(rule, Temp, fill(average_field, n))
            push!(ED_gap, gap)
            # @info "n = $n, average_field = $(average_field), gap = $(eigvals[1] - eigvals[2]), converged = $(infos.converged), numiter = $(infos.numiter)"
        end
        x = [i for i in sample_w]
        logx = x
        logy = real.(ED_gap)
        scatter!(ax, logx, logy, label = "h = $average_field")
        # fit = curve_fit(LinearFit, logx, logy)
        # push!(ED_average_slope, fit.coefs[2])
        @info "average_field = $(average_field)"
    end

    axislegend(position = :rb)
    f
end
plot_spectral_gap(HeatBath(), 10)
plot_spectral_gap(Metroplis(), 10)
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

tostring(n::Int, idx::Int) = bitstring(idx)[end-2n+3:end]
function plot_spectrum(rule::TransitionRule, Temp, on_site_energy, k::Int)
    n = length(on_site_energy)
    evals, evecs = spectral(rule, Temp, on_site_energy, k)
    @info evals
    barplot(real.(evecs[k]); axis = (xticks=(1:2^(2n-2), tostring.(length(on_site_energy), 0:(2^(2n-2) - 1))), title=""))
end

plot_spectrum(Metroplis(), 0.3, fill(1.0, 3), 6)
plot_spectrum(HeatBath(), 0.3, fill(1.0, 4), 16)