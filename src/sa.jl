function rule110(p, q, r)
    return (q + r + q*r + p*q*r) % 2
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
@testset "transition matrix" begin
    P = transition_matrix(Metroplis(), 2.0, fill(0.5, 4))
    @test size(P) == (2^6, 2^6)
    @test all(≈(1), sum(P; dims=1))

    P = transition_matrix(HeatBath(), 2.0, fill(0.5, 4))
    @test size(P) == (2^6, 2^6)
    @test all(≈(1), sum(P; dims=1))
end

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

tostring(n::Int, idx::Int) = bitstring(idx)[end-2n+3:end]
function plot_spectrum(rule::TransitionRule, Temp, on_site_energy, k::Int)
    n = length(on_site_energy)
    evals, evecs = spectral(rule, Temp, on_site_energy, k)
    @info evals
    barplot(real.(evecs[k]); axis = (xticks=(1:2^(2n-2), tostring.(length(on_site_energy), 0:(2^(2n-2) - 1))), title=""))
end

plot_spectrum(Metroplis(), 0.3, fill(1.0, 3), 6)
plot_spectrum(HeatBath(), 0.3, fill(1.0, 4), 16)