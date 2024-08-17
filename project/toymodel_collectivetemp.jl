using QuantumAdiabaticAnnealing
using QuantumAdiabaticAnnealing: random_state, calculate_energy
using KrylovKit
using SparseArrays
using CairoMakie
using CurveFit
using Random


function generic_logic_grate(p, q, r, N)
    return (N >> (p << 2 | q << 1 | r)) & 1
end

function myenergy(width, depth, state, onsite, inputs, outputs)
    sa = SimulatedAnnealingHamiltonian(width, depth)
    statearray = zeros(Bool, (width * (depth), 1))
    for i in 1:width*depth
        statearray[i, 1] = Bool((state >> (i-1))&1)
    end
    # for i in 1:width
    #     statearray[i, 1] = Bool((inputs >> (i-1))&1)
    # end
    # for i in 1:width
    #     statearray[i+width*(depth+1), 1] = Bool((outputs >> (i-1))&1)
    # end
    # @info statearray
    ret = calculate_energy(sa, statearray, ones((1,1)), 1)
    for i in 1:width
        ret += onsite * (statearray[i, 1] ? 1 : -1)
    end
    return ret
end

function trueoutput(width, depth, inputs) # inputs is the first layer calculate depth layer
    statearray = zeros(Bool, (width, 1))
    next_array = zeros(Bool, (width, 1))
    for i in 1:width
        statearray[i] = Bool((inputs >> (i-1))&1)
    end
    # @info statearray
    for i in 2:depth
        for mid in 1:width
            pre = mod1(mid-1, width)
            suf = mod1(mid+1, width)
            next_array[mid] = generic_logic_grate(statearray[pre], statearray[mid], statearray[suf], 110)
        end
        # @info next_array
        statearray = copy(next_array)
    end
    # @info statearray
    return sum(statearray .* [2^(i-1) for i in 1:width])
end

function transition_matrix(width, depth, onsite, temp)
    row = Vector{Int}()
    col = Vector{Int}()
    val = Vector{Float64}()
    # inputs = rand(Vector(0:2^width-1))
    inputs = 3
    outputs = trueoutput(width, depth+2, inputs)
    # inputs = 0
    allenergy = [myenergy(width, depth, i, onsite, inputs, outputs) for i in 0:2^(width * depth)-1]
    # @info allenergy
    for pre_state in 0:2^(width * depth)-1 
        other_prob = 0.0
        for change_id in 1:width*depth
            next_state = pre_state ⊻ (1 << (change_id - 1))
            push!(row, next_state + 1)
            push!(col, pre_state + 1)
            push!(val, 1.0 / (1.0 + exp((allenergy[next_state + 1] - allenergy[pre_state + 1]) / temp)) * 1.0 / (width * depth))
            other_prob += val[end]
        end
        # @info other_prob
        push!(row, pre_state + 1)
        push!(col, pre_state + 1)
        push!(val, 1 - other_prob)
    end
    return sparse(row, col, val)
end

function spectral_gap(width, depth, onsite, temp)
    @info "now calculating spectral gap for width = $width, depth = $depth, onsite = $onsite, temp = $temp"
    tmatrix = transition_matrix(width, depth, onsite, temp)
    # @info size(tmatrix)
    @info "finish calculate transition_matrix"
    # @info round.(Matrix(tmatrix), digits=2)
    eigvals, eigvecs, infos = eigsolve(tmatrix, rand(Float64, size(tmatrix, 2)), 2, :LM; maxiter = 5000)
    @info eigvals
    @info infos.numiter, infos.converged
    @assert infos.converged >= 2 "not converged"
    @assert eigvals[1] ≈ 1
    return abs(eigvals[1]) - abs(eigvals[2])
end

function plot(ydata, xdata)
    ydata = (log.(ydata ./ xdata))
    xdata = log.(Float64.(xdata))
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "log(depth)", ylabel = "log(gap)",
        title = "Log(Gap) v.s. log(Depth)")
    scatter!(ax, xdata, ydata)

    fit = curve_fit(LinearFit, xdata, ydata)
    lines!(ax, xdata, fit.(xdata))
    # ylims!(ax, low=0)
    f
end
gap = [spectral_gap(4, i, 0.0, 0.2) for i in 2:5]
plot(gap, Vector(2:5))

