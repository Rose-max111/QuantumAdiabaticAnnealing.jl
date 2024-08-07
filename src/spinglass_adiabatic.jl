struct WeightedEdge{T<:Real} <: Graphs.AbstractEdge{Int}
    src::Int
    dst::Int
    weight::T
end

struct SpinGlassModel{T<:Real}
    n::Int
    m::Int
    gradient::T
    edges::Vector{WeightedEdge{T}}
    onsite::Vector{T}
    function SpinGlassModel(n::Int, m::Int, gradient::T, edges::Vector{WeightedEdge{T}}, onsite::Vector{T}) where T
        @assert length(onsite) == 2*n*m-n "onsite energy length should be $(2*m*n-n), but got $(length(onsite))"
        new{T}(n, m, gradient, edges, onsite)
    end
end

nspin(sp::SpinGlassModel) = 2 * sp.n * sp.m - sp.n

struct SpinGlassModelCache{T<:Real}
    M::Vector{Point3D{T}}  # magnetic momentum for each site
    Mc::Vector{Point3D{T}}  # cache for magnetic momentum
    inst_H::Vector{Point3D{T}}   # instantaneous field for each site
    ks::Vector{Vector{Point3D{T}}}  # ODE cache, e.g. k1, k2, k3, k4 for Runge-Kutta
end

function init_cache(sp::SpinGlassModel, nk::Int)
    M = vcat([Point(0.0, 0.0, 1.0) for i in 1:sp.n], [Point(-1.0, 0.0, 0.0) for i in 1:nspin(sp)-sp.n])
    return SpinGlassModelCache(M, similar(M), similar(M), [similar(M) for i in 1:nk])
end

function spinglass_random_mapping(n::Int, interaction_part::Vector{Tuple{Int, Int}}, interaction_weight::Vector{T}, onsite_part=fill(0.0, n)) where T
    edges = Vector{WeightedEdge{Int,Int,T}}()
    for i in 1:length(interaction_part)
        push!(edges, (interaction_part[i][1], interaction_part[i][2], interaction_weight[i]))
    end|
    return SpinGlassModel(n, 1, one(T), edges, onsite_part)
end

# m layer, per layer with n gadgets
# consist output layer
# period condition
function spinglass_mapping(n::Int, m::Int; gradient = 1.0)
    single_onsite = [1.0, 2.0, 2.0, 2.0, 5.0] # last one is ancilla
    single_edge_weights = [1.0, 1.0, 2.0, 3.0, 2.0, 2.0, 5.0, 2.0, 5.0, 6.0]
    # ancilla id is n*m+1 to n*m+n*(m-1)
    edges = WeightedEdge{Float64}[]
    onsites = zeros(n*m + n*(m-1))
    lis = LinearIndices((n, m))
    for layer in 1:m-1
        for mid in 1:n
            pre = mod1(mid-1, n)
            nxt = mod1(mid+1, n)
            this_id = [lis[pre, layer], lis[mid, layer], lis[nxt, layer],
                     lis[mid, layer+1], lis[mid, layer] + n*m]
            # @info "this_id = $this_id"
            cnt = 0
            for i in 1:5
                for j in i+1:5
                    cnt += 1
                    push!(edges, WeightedEdge(this_id[i], this_id[j], single_edge_weights[cnt] * (gradient^(m-1 - layer))))
                end
            end
            for i in 1:5
                onsites[this_id[i]] += single_onsite[i] * (gradient^(m-1 - layer))
            end
        end
    end
    return SpinGlassModel(n, m, gradient, edges, onsites)
end

function field!(out::AbstractVector{Point3D{T}}, M::Vector{Point3D{T}}, H::Vector{Point3D{T}}) where T
    # H(i) = -hx̂ + ∑_j J_{i,j} M_j^z ẑ + onsite ẑ
    # ̇Ṁ = H \cross M
    out .= cross.(M, H)
end

function instantaneous_field!(H::AbstractVector{Point3D{ET}}, sp::SpinGlassModel, M, t, T, Vtrans::Vector{ET}) where ET
    H .= Point.(Vtrans .* (-(T-t)/T), 0.0, sp.onsite .* (-t/T))  # onsite field
    for edge in sp.edges
        H[edge.src] += Point(0.0, 0.0, -edge.weight * M[edge.dst][3] * (t/T))
        H[edge.dst] += Point(0.0, 0.0, -edge.weight * M[edge.src][3] * (t/T))
    end
    return H
end

function field!(Mdot::AbstractVector{Point3D{ET}}, sp::SpinGlassModel, cache::SpinGlassModelCache, t, T, Vtrans::Vector{Float64}; pin_input::Bool = true) where ET # evaluate F(t, y)
    pin_input && fill!(view(cache.M, Base.OneTo(sp.n)), Point(0.0, 0.0, 1.0))

    instantaneous_field!(cache.inst_H, sp, cache.M, t, T, Vtrans)
    field!(Mdot, cache.M, cache.inst_H)

    if pin_input == true   # ??? why 0, 0, 0?
        Mdot[1:sp.n] .= Ref(Point(0.0, 0.0, 0.0))
    end

    return Mdot
end

function runge_kutta_singlejump!(sp::SpinGlassModel, cache::SpinGlassModelCache, t0, delta_t, T, Vtrans::Vector{Float64}; pin_input = true)
    cache.Mc .= cache.M
    field!(cache.ks[1], sp, cache, t0, T, Vtrans; pin_input)
    cache.Mc .= normalize.(cache.M .+ cache.ks[1] .* (delta_t / 2))
    field!(cache.ks[2], sp, cache, t0 + delta_t / 2, T, Vtrans; pin_input)
    cache.Mc .= normalize.(cache.M .+ cache.ks[2] .* (delta_t / 2))
    field!(cache.ks[3], sp, cache, t0 + delta_t / 2, T, Vtrans; pin_input)
    cache.Mc .= normalize.(cache.M .+ cache.ks[3] .* delta_t)
    field!(cache.ks[4], sp, cache, t0 + delta_t, T, Vtrans; pin_input)
    cache.M .= normalize.(cache.M .+ (cache.ks[1] .+ 2 .* cache.ks[2] .+ 2 .* cache.ks[3] .+ cache.ks[4]) .* (delta_t / 6))
end

function runge_kutta_integrate(sp::SpinGlassModel, dt::Float64, T::Float64, Vtrans::Vector{Float64}; T_end = nothing, pin_input = true)
    t0 = 0.0
    if T_end === nothing
        T_end = T
    end
    sp_print = Vector{Point3D{Float64}}[]
    H_print = Vector{Point3D{Float64}}[]
    cache = init_cache(sp, 4)
    while t0 < T_end
        delta_t = min(dt, T_end-t0)
        # save data
        push!(sp_print, copy(cache.M))
        instantaneous_field!(cache.inst_H, sp, cache.M, t0, T, Vtrans)
        push!(H_print, copy(cache.inst_H))

        # evolve
        runge_kutta_singlejump!(sp, cache, t0, delta_t, T, Vtrans; pin_input)
        t0 += delta_t
    end
    return sp_print, H_print
end

function euclidean_integrate(sp::SpinGlassModel, dt, T, Vtrans::Vector{Float64}; pin_input = true)
    t=0
    cache = init_cache(sp, 0)
    while t<T
        delta_t = min(dt, T-t)
        field!(cache.Mc, sp, cache, t, T, Vtrans; pin_input)
        cache.M .+= cache.Mc .* delta_t
        t += delta_t
    end
    return cache.M
end

function sp_ground_state(sp::SpinGlassModel)
    hyperedges = [[t.src, t.dst] for t in sp.edges]
    hyperweights = [t.weight for t in sp.edges]
    for i in 1:length(sp.onsite)
        push!(hyperedges, [i])
        push!(hyperweights, sp.onsite[i])
    end
    spgls = SpinGlass(length(sp.onsite), hyperedges, hyperweights)
    spproblem = GenericTensorNetwork(spgls)
    gs = solve(spproblem, SingleConfigMin())[]
    return gs
end

function sp_sa_step!(sp::SpinGlassModel, M, Temp, node)
    theta = round(rand() * π,digits=4)
    ori_M = copy(M[node])
    pre_E = sp_energy(sp, M, 1, 1, zeros(length(sp.onsite)))
    M[node] = [sin(theta), 0.0, cos(theta)]
    post_E = sp_energy(sp, M, 1, 1, zeros(length(sp.onsite)))
    
    delta_E = post_E - pre_E
    accept_prob = delta_E < 0 ? 1.0 : exp(-delta_E / Temp)
    if rand() <= accept_prob
        return true
    else
        M[node] = ori_M
        return false
    end
end
    
function sp_groud_state_sa_single(sp::SpinGlassModel)
    Temp = 15.0
    for i in 1:length(sp.onsite)
        theta = round(rand() * π, digits=4)
        M[i] = [sin(theta), 0.0, cos(theta)]
    end
    for t in Temp:-0.1:1e-2
        for epoch in 1:1000
            this_node = rand(Vector(1:length(sp.onsite)))
            sp_sa_step!(sp, M, t, this_node)
        end
    end
    return sp_energy(sp, M, 1, 1, zeros(length(sp.onsite)))
end

function sp_ground_state_sa(n, m)
    minn = Inf
    for epoch in 1:1000
        sp = spinglass_mapping(n, m)
        minn = min(minn, sp_groud_state_sa_single(sp))
        @info "epoch = $epoch, minn = $minn"
    end
    return minn
end

function sp_energy(sp::SpinGlassModel, M::AbstractVector{Point3D{ET}}, t, T, Vtrans) where ET
    ret = zero(ET)
    ret += sum([t[3] for t in M] .* sp.onsite) * (t/T)
    ret += sum([t[1] for t in M] .* Vtrans) * ((T-t)/T)
    for edge in sp.edges
        ret += edge.weight * M[edge.src][3] * M[edge.dst][3] * (t/T)
    end
    return ret
end

function sp_check_valid(rule::CellularAutomata1D, sp::SpinGlassModel)
    for i in 1:length(sp.onsite)
        if 1 - abs(sp.M[i][3]) > 1e-1
            @info "not satisfy reach equilibrium state"
            return false
        end
    end
    # spin > 0 means 0 in generictensorwork, mean 1 in rule110
    # spin < 0 means 1 in generictensorwork, mean 0 in rule110
    state = [sp.M[i][3] > 0 ? 1 : 0 for i in 1:length(sp.onsite)]
    lis = LinearIndices((sp.n, sp.m))
    for j in 1:sp.m-1
        for i in 1:sp.n
            pre_l = lis[mod1(i-1, sp.n), j]
            mid_l = lis[i, j]
            suf_l = lis[mod1(i+1, sp.n), j]
            @info "i = $i, j = $j"
            if logic_gate(rule, state[pre_l], state[mid_l], state[suf_l]) != state[lis[i, j+1]]
                return false
            end
        end
    end
    return true
end

wrap_data(vec::AbstractVector{T}) where T = wrap_data!(Vector{Point3D{T}}(undef, length(vec) ÷ 3), vec)
function wrap_data!(M::Vector{Point3D{T}}, vec::AbstractVector{T}) where T
    for i=1:length(M)
        M[i] = Point(vec[3i-2], vec[3i-1], vec[3i])
    end
    return M
end
unwrap_data(M::Vector{Point3D{T}}) where T = unwrap_data!(Vector{T}(undef, 3*length(M)), M)
function unwrap_data!(vec::AbstractVector{T}, M::Vector{Point3D{T}}) where T
    for i=1:length(M)
        vec[3i-2] = M[i][1]
        vec[3i-1] = M[i][2]
        vec[3i] = M[i][3]
    end
    return vec
end

function spinglass_adiabatic_dp8(sg::SpinGlassModel, Tmax; pin_input = true)
    cache = init_cache(sg, 0)
    init_dp8 = unwrap_data(cache.M)
    solver8 = DP8Solver(0.0, init_dp8; atol=1e-10, rtol=1e-10, maximum_allowed_steps=5000000) do t, y, field
        vectorgradient!(field, sg, cache, y, t; pin_input, Tmax)
    end
    integrate!(solver8, Tmax)
    return wrap_data(get_current_state(solver8))
end

function vectorgradient!(out::AbstractVector, sp::SpinGlassModel, cache::SpinGlassModelCache, vec::AbstractVector{T}, this_t; pin_input, Tmax) where T
    wrap_data!(cache.M, vec)
    field!(cache.Mc, sp, cache, this_t, Tmax, sp.onsite; pin_input)
    return unwrap_data!(out, cache.Mc)
end