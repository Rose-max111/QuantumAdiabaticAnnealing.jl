using GenericTensorNetworks

mutable struct spinglassmodel
    # n::Int
    const n::Int
    const m::Int
    const edges::Vector{Tuple{Int,Int,Float64}}
    const onsite::Vector{Float64}
    M::Vector{Vector{Float64}}
    function spinglassmodel(n::Int, m::Int, edges::Vector{Tuple{Int,Int,Float64}}, onsite::Vector{Float64}, M::Vector{Vector{Float64}} = fill([1.0, 0.0, 0.0], length(onsite)))
        new(n, m, edges, onsite, M)
    end
end

# m layer, per layer with n gadgets
# consist output layer
# period condition
function cartesian_to_linear(i::Int, j::Int, n::Int)
    return i + (j-1) * n
end

function spinglass_mapping(n::Int, m::Int)
    single_onsite = [1.0, 2.0, 2.0, 2.0, 5.0] # last one is ancilla
    single_edge_weights = [1.0, 1.0, 2.0, 3.0, 2.0, 2.0, 5.0, 2.0, 5.0, 6.0]
    # ancilla id is n*m+1 to n*m+n*(m-1)
    edges = Vector{Tuple{Int,Int,Float64}}()
    onsites = zeros(n*m + n*(m-1))
    for layer in 1:m-1
        for mid in 1:n
            pre = mod1(mid-1, n)
            nxt = mod1(mid+1, n)
            this_id = [cartesian_to_linear(pre, layer, n), cartesian_to_linear(mid, layer, n), cartesian_to_linear(nxt, layer, n),
                     cartesian_to_linear(mid, layer+1, n), cartesian_to_linear(mid, layer, n) + n*m]
            @info "this_id = $this_id"
            cnt = 0
            for i in 1:5
                for j in i+1:5
                    cnt += 1
                    push!(edges, (this_id[i], this_id[j], single_edge_weights[cnt]))
                end
            end
            for i in 1:5
                onsites[this_id[i]] += single_onsite[i]
            end
        end
    end
    return spinglassmodel(n, m, edges, onsites)
end


# H(i) = -hx̂ + ∑_j J_{i,j} M_j^z ẑ + onsite ẑ
# ̇M = H \cross M
function integrator(sp::spinglassmodel, t, T, Vtrans)
    Mxdot = - sp.onsite .* [u[2] for u in sp.M] .* (t/T)
    Mydot = sp.onsite .* [u[1] for u in sp.M] .* (t/T)
    for edge in sp.edges
        i, j, w = edge
        Mydot[i] += w * sp.M[j][3] * sp.M[i][1] * (t/T)
        Mydot[j] += w * sp.M[i][3] * sp.M[j][1] * (t/T)

        Mxdot[i] -= w * sp.M[j][3] * sp.M[i][2] * (t/T)
        Mxdot[j] -= w * sp.M[i][3] * sp.M[j][2] * (t/T)
    end
    Mydot -= - Vtrans .* [u[3] for u in sp.M] .* ((T-t)/T)
    Mzdot = - Vtrans .* [u[2] for u in sp.M] .* ((T-t)/T)

    return Mxdot, Mydot, Mzdot
end

function renorm(point)
    norm = sqrt(sum(point .* point))
    return point ./ norm
end

ss = spinglassmodel(1, 1, Vector{Tuple{Int64,Int64,Float64}}(), [2.0])
tt = 0
TT = round(π/4 , digits = 4)
while tt < TT
    delta_t = min(1e-4, TT-tt)
    Mxdot, Mydot, Mzdot = integrator(ss, 1, 1, zeros(1))
    # Mxdot1, Mydot1, Mzdot1 = integrator(sp, t/5, T/5, Vtrans)
    # @assert round.((Mxdot1 .-  Mxdot), digits=2) == zeros(length(sp.onsite))
    # @assert round.((Mydot1 .-  Mydot), digits=2) == zeros(length(sp.onsite))
    # @assert round.((Mzdot1 .-  Mzdot), digits=2) == zeros(length(sp.onsite))
    d_Mx = Mxdot .* delta_t
    d_My = Mydot .* delta_t
    d_Mz = Mzdot .* delta_t
    ss.M = [renorm([ss.M[i][1]+d_Mx[i], ss.M[i][2]+d_My[i], ss.M[i][3]+d_Mz[i]]) for i in 1:length(ss.onsite)]
    tt += delta_t
end
println(ss.M)



function sp_energy(sp::spinglassmodel, t, T, Vtrans)
    ret = 0.0
    ret += sum([t[3] for t in sp.M] .* sp.onsite) * (t/T)
    ret -= sum([t[1] for t in sp.M] .* Vtrans) * ((T-t)/T)
    for edge in sp.edges
        i, j, w = edge
        ret += w * sp.M[i][3] * sp.M[j][3] * (t/T)
    end
    return ret
end


function integrate!(sp::spinglassmodel, dt, T, Vtrans)
    t=0
    while(t<T)
        delta_t = min(dt, T-t)
        Mxdot, Mydot, Mzdot = integrator(sp, t, T, Vtrans)
        # Mxdot1, Mydot1, Mzdot1 = integrator(sp, t/5, T/5, Vtrans)
        # @assert round.((Mxdot1 .-  Mxdot), digits=2) == zeros(length(sp.onsite))
        # @assert round.((Mydot1 .-  Mydot), digits=2) == zeros(length(sp.onsite))
        # @assert round.((Mzdot1 .-  Mzdot), digits=2) == zeros(length(sp.onsite))
        d_Mx = Mxdot .* delta_t
        d_My = Mydot .* delta_t
        d_Mz = Mzdot .* delta_t
        sp.M = [renorm([sp.M[i][1]+d_Mx[i], sp.M[i][2]+d_My[i], sp.M[i][3]+d_Mz[i]]) for i in 1:length(sp.onsite)]
        t += delta_t
        # @info "t = $t, energy = $(sp_energy(sp, t, T, Vtrans))"
    end
end


sp = spinglass_mapping(3, 2)

# hyperedges = [[t[1], t[2]] for t in sp.edges]
# for t in 1:length(sp.onsite)
#     push!(hyperedges, [t])
# end
# hyperweights = [t[3] for t in sp.edges]
# append!(hyperweights, sp.onsite)

# spgls = SpinGlass(length(sp.onsite), hyperedges, hyperweights)
# problem = GenericTensorNetwork(spgls)
# configs = solve(problem, ConfigsMin())[]

Vtrans = fill(1, length(sp.onsite))

integrate!(sp, 0.00001, 1, Vtrans)

println(sp.M)
