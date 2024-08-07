abstract type TransitionRule end
struct HeatBath <: TransitionRule end
struct Metropolis <: TransitionRule end

# Metropolis-Hastings algorithm
function update(::Metropolis, temperature, ΔE, prior)
    if ΔE < 0
        1.0
    else
        exp(- (ΔE) / temperature)
    end * prior
end

function update(::HeatBath, temperature, ΔE, prior)
    exp(-ΔE / temperature) / (1 + exp(-ΔE / temperature)) * prior
end

function toy_model_state_energy(rule::CellularAutomata1D, state, n::Int, m::Int; on_site_energy, period_condition::Bool, energy_gradient)
    ret = 0
    on_site_energy !== nothing && for l in 1:n
        ret += on_site_energy[l] * (readbit(state, l-1) > 0 ? 1 : -1)
    end
    lis = LinearIndices((n, m))
    for row in 1:m-1, mid in (period_condition ? (1:n) : (2:n-1))
        p, q, r, output = lis[mod1(mid-1, n), row], lis[mid, row], lis[mod1(mid+1, n), row], lis[mid, row+1]
        local_energy = logic_gate(rule, readbit(state, p), readbit(state, q), readbit(state, r)) ⊻ readbit(state, output)
        ret += local_energy * energy_gradient^(m - row - 1)
    end
    return ret
end

function toy_model_transition_matrix(ca::CellularAutomata1D, rule::TransitionRule, n::Int, m::Int, temperature; on_site_energy, period_condition::Bool, energy_gradient)
    total_atoms = period_condition ? n * m : n * m - 2   # ?
    state_energy = [toy_model_state_energy(ca, i, n, m; on_site_energy, period_condition, energy_gradient) for i in 0:(2^total_atoms - 1)]

    row = Vector{Int}()
    col = Vector{Int}()
    val = Vector{Float64}()
    @info "finish calculating state energy"
    for i in 0:(2^total_atoms - 1)
        other_prob = 0
        for j in 1:total_atoms
            reverse_mask = i ⊻ (2^(j-1))
            push!(row, reverse_mask + 1)
            push!(col, i+1)

            ΔE = state_energy[reverse_mask + 1] - state_energy[i + 1]
            push!(val, update(rule, temperature, ΔE, 1.0 / total_atoms))
            other_prob += val[end]
        end
        push!(row, i+1)
        push!(col, i+1)
        push!(val, 1 - other_prob)
    end
    return sparse(row, col, val)
end

function spectral(tmatrix, k::Int)
    eigvals, eigvecs, info = eigsolve(tmatrix, rand(Float64, size(tmatrix, 2)), k, :LM; maxiter = 5000)
    @assert info.converged >= k "not converged"
    @assert eigvals[1] ≈ 1
    @assert all(isapprox(ev, real(ev); atol=1e-5) for ev in eigvals) "complex eigvals: $eigvals"
    return real.(eigvals), eigvecs
end

function plot_spectral_gap(ca::CellularAutomata1D, rule::TransitionRule, temperature;
            on_site_energies = 2 .^ (0:8),
            sample_w = 3:9,
            m = 3
        )
    f = Figure()
    ax = Axis(f[1,1], title = "Spectral Gap v.s. Width", xscale=log10, yscale=log10, xlabel = "Width", ylabel = "Spectral Gap")
    for average_field in on_site_energies
        ED_gap = []
        for n in sample_w
            tmatrix = toy_model_transition_matrix(ca, rule, n, m, temperature, on_site_energy = fill(average_field, n))
            evals, evecs = spectral(tmatrix, 2)
            gap = abs(abs(evals[1]) - abs(evals[2]))
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

tostring(n::Int, idx::Int) = bitstring(idx)[end-2n+3:end]
function plot_spectrum(ca::CellularAutomata1D, rule::TransitionRule, temperature, on_site_energy, k::Int, m::Int)
    n = length(on_site_energy)
    tmatrix = toy_model_transition_matrix(ca, rule, n, m, temperature, on_site_energy = on_site_energy)
    evals, evecs = spectral(tmatrix, k)
    @info evals
    barplot(real.(evecs[k]); axis = (xticks=(1:2^(2n-2), tostring.(length(on_site_energy), 0:(2^(2n-2) - 1))), title=""))
end