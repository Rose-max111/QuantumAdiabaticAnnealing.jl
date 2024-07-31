using GenericTensorNetworks

mutable struct spinglassmodel
    # n::Int
    const n::Int
    const m::Int
    const edges::Vector{Tuple{Int,Int,Float64}}
    const onsite::Vector{Float64}
    theta::Vector{Float64}
    theta_dot::Vector{Float64}
    function spinglassmodel(n::Int, m::Int, edges::Vector{Tuple{Int,Int,Float64}}, onsite::Vector{Float64}, theta::Vector{Float64} = fill(π/2, length(onsite)), theta_dot::Vector{Float64} = zeros(length(onsite)))
        new(n, m, edges, onsite, theta, theta_dot)
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

function integrator(sp::spinglassmodel, t, T, Vtrans)
    delta_theta_dot = - sin.(sp.theta) .* sp.onsite
    for edge in sp.edges
        i, j, w = edge
        delta_theta_dot[i] -= w * sin(sp.theta[i]) * cos(sp.theta[j])
        delta_theta_dot[j] -= w * cos(sp.theta[i]) * sin(sp.theta[j])
    end
    delta_theta_dot = (delta_theta_dot .* (t/T)) .- (Vtrans .* cos.(sp.theta) .* ((T-t)/T))
    return delta_theta_dot
end

function sp_energy(sp::spinglassmodel, t, T, Vtrans)
    ret = 0.0
    ret += sum(sp.onsite .* (cos.(sp.theta))) * (t/T)
    ret -= sum(Vtrans .* sin.(sp.theta)) * (T-t)/T
    for edge in sp.edges
        i, j, w = edge
        ret += w * cos(sp.theta[i]) * cos(sp.theta[j]) * (t/T)
    end
    return ret
end

function integrate!(sp::spinglassmodel, dt, T, Vtrans)
    t=0
    while(t<T)
        delta_t = min(dt, T-t)
        delta_theta_dot = integrator(sp, t, T, Vtrans)
        sp.theta = sp.theta + sp.theta_dot .* delta_t
        sp.theta_dot = sp.theta_dot .+ delta_theta_dot .* delta_t
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
# configs = solve(problem, SingleConfigMin())[]

Vtrans = fill(2.0, length(sp.onsite))

integrate!(sp, 0.00001, 1, Vtrans)

cos.(sp.theta)
# sp.theta ./ 2π
