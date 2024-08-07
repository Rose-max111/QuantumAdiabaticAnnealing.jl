struct SpinGlassModel
    n::Int
    m::Int
    gradient::Float64
    edges::Vector{Tuple{Int,Int,Float64}}
    onsite::Vector{Float64}
    M::Vector{Point3D{Float64}}
    function SpinGlassModel(n::Int, m::Int, gradient::Float64, edges::Vector{Tuple{Int,Int,Float64}}, onsite::Vector{Float64}, M::Vector{Point3D{Float64}} = cat([Point(0.0, 0.0, 1.0) for i in 1:n], [Point(-1.0, 0.0, 0.0) for i in 1:length(onsite)-n];dims=1))
        new(n, m, gradient, edges, onsite, M)
    end
end

function spinglass_random_mapping(n, interaction_part::Vector{Tuple{Int, Int}}, interaction_weight::Vector{Float64}, onsite_part=fill(0.0, n))
    edges = Vector{Tuple{Int,Int,Float64}}()
    for i in 1:length(interaction_part)
        push!(edges, (interaction_part[i][1], interaction_part[i][2], interaction_weight[i]))
    end|
    return SpinGlassModel(n, 1, 1.0, edges, onsite_part, [Point(-1.0, 0.0, 0.0) for i in 1:n])
end
        

# m layer, per layer with n gadgets
# consist output layer
# period condition
function cartesian_to_linear(i::Int, j::Int, n::Int)
    return i + (j-1) * n
end

function spinglass_mapping(n::Int, m::Int; gradient = 1.0)
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
            # @info "this_id = $this_id"
            cnt = 0
            for i in 1:5
                for j in i+1:5
                    cnt += 1
                    push!(edges, (this_id[i], this_id[j], single_edge_weights[cnt] * (gradient^(m-1 - layer))))
                end
            end
            for i in 1:5
                onsites[this_id[i]] += single_onsite[i] * (gradient^(m-1 - layer))
            end
        end
    end
    return SpinGlassModel(n, m, gradient, edges, onsites)
end

# H(i) = -hx̂ + ∑_j J_{i,j} M_j^z ẑ + onsite ẑ
# ̇M = H \cross M
function crossdot(A::Point3D{Float64}, B::Point3D{Float64})
    return Point(A[2]*B[3] - A[3]*B[2], A[3]*B[1] - A[1]*B[3], A[1]*B[2] - A[2]*B[1])
end

function integrator(sp::SpinGlassModel, H::Vector{Point3D{Float64}})
    Mdot = crossdot.(sp.M, H)
    return Mdot
end

function instantaneous_field(sp::SpinGlassModel, t, T, Vtrans::Vector{Float64})
    H = [Point(-Vtrans[i] .* ((T-t)/T), 0.0, -sp.onsite[i] .* (t/T)) for i in 1:length(sp.onsite)]
    for edge in sp.edges
        i, j, w = edge
        H[i] += Point(0.0, 0.0, -w * sp.M[j][3] * (t/T))
        H[j] += Point(0.0, 0.0, -w * sp.M[i][3] * (t/T))
    end
    return H
end

function integrator(sp::SpinGlassModel, t, T, Vtrans::Vector{Float64}; pin_input = true) # evaluate F(t, y)
    if pin_input == true
        for i in 1:sp.n
            sp.M[i] = Point(0.0, 0.0, 1.0)
        end # used set-way
    end

    inst_H = instantaneous_field(sp, t, T, Vtrans)
    exact_Mdot = integrator(sp, inst_H)

    Mxdot = [exact_Mdot[i][1] for i in 1:length(sp.onsite)]
    Mydot = [exact_Mdot[i][2] for i in 1:length(sp.onsite)]
    Mzdot = [exact_Mdot[i][3] for i in 1:length(sp.onsite)]

    if pin_input == true
        Mxdot[1:sp.n] .= 0.0
        Mydot[1:sp.n] .= 0.0
        Mzdot[1:sp.n] .= 0.0
    end

    return [Mxdot, Mydot, Mzdot]
end

function runge_kutta_singlejump!(sp::SpinGlassModel, t0, delta_t, T, Vtrans::Vector{Float64}; pin_input = true)
    origin_M = copy(sp.M)
    k1 = integrator(sp, t0, T, Vtrans;pin_input=pin_input)
    sp.M .= [normalize(Point(origin_M[i][1] + k1[1][i] * delta_t / 2, origin_M[i][2] + k1[2][i] * delta_t / 2, origin_M[i][3] + k1[3][i] * delta_t / 2)) for i in 1:length(sp.onsite)]
    k2 = integrator(sp, t0 + delta_t / 2, T, Vtrans; pin_input=pin_input)
    sp.M .= [normalize(Point(origin_M[i][1] + k2[1][i] * delta_t / 2, origin_M[i][2] + k2[2][i] * delta_t / 2, origin_M[i][3] + k2[3][i] * delta_t / 2)) for i in 1:length(sp.onsite)]
    k3 = integrator(sp, t0 + delta_t / 2, T, Vtrans; pin_input=pin_input)
    sp.M .= [normalize(Point(origin_M[i][1] + k3[1][i] * delta_t, origin_M[i][2] + k3[2][i] * delta_t, origin_M[i][3] + k3[3][i] * delta_t)) for i in 1:length(sp.onsite)]
    k4 = integrator(sp, t0 + delta_t, T, Vtrans; pin_input=pin_input)
    real_k = (k1 .+ 2*k2 .+ 2*k3 .+ k4) ./ 6 * delta_t
    sp.M .= [normalize(Point(origin_M[i][1] .+ real_k[1][i], origin_M[i][2]+real_k[2][i], origin_M[i][3]+real_k[3][i])) for i in 1:length(sp.onsite)]
end

function runge_kutta_integrate!(sp::SpinGlassModel, dt::Float64, T::Float64, Vtrans::Vector{Float64}; T_end = nothing, pin_input = true)
    t0 = 0.0
    if T_end === nothing
        T_end = T
    end
    sp_print = Vector{Point3D{Float64}}()
    H_print = Vector{Point3D{Float64}}()
    while t0 < T_end
        delta_t = min(dt, T_end-t0)
        append!(sp_print, copy(sp.M))
        append!(H_print, instantaneous_field(sp, t0, T, Vtrans))

        runge_kutta_singlejump!(sp, t0, delta_t, T, Vtrans; pin_input = pin_input)
        t0 += delta_t
    end
    return sp_print, H_print
end

function euclidean_integrate!(sp::SpinGlassModel, dt, T, Vtrans::Vector{Float64})
    t=0
    while t<T
        delta_t = min(dt, T-t)
        Mxdot, Mydot, Mzdot = integrator(sp, t, T, Vtrans)
        d_Mx = Mxdot .* delta_t
        d_My = Mydot .* delta_t
        d_Mz = Mzdot .* delta_t
        sp.M .= [normalize(Point(sp.M[i][1]+d_Mx[i], sp.M[i][2]+d_My[i], sp.M[i][3]+d_Mz[i])) for i in 1:length(sp.onsite)]
        t += delta_t
    end
end

function sp_ground_state(sp::SpinGlassModel)
    hyperedges = [[t[1],t[2]] for t in sp.edges]
    hyperweights = [t[3] for t in sp.edges]
    for i in 1:length(sp.onsite)
        push!(hyperedges, [i])
        push!(hyperweights, sp.onsite[i])
    end
    spgls = SpinGlass(length(sp.onsite), hyperedges, hyperweights)
    spproblem = GenericTensorNetwork(spgls)
    gs = solve(spproblem, SingleConfigMin())[]
    return gs
end

function sp_sa_step!(sp::SpinGlassModel, Temp, node)
    theta = round(rand() * π,digits=4)
    ori_M = copy(sp.M[node])
    pre_E = sp_energy(sp, 1, 1, zeros(length(sp.onsite)))
    sp.M[node] = [sin(theta), 0.0, cos(theta)]
    post_E = sp_energy(sp, 1, 1, zeros(length(sp.onsite)))
    
    delta_E = post_E - pre_E
    accept_prob = delta_E < 0 ? 1.0 : exp(-delta_E / Temp)
    if rand() <= accept_prob
        return true
    else
        sp.M[node] = ori_M
        return false
    end
end
    
function sp_groud_state_sa_single(sp::SpinGlassModel)
    Temp = 15.0
    for i in 1:length(sp.onsite)
        theta = round(rand() * π, digits=4)
        sp.M[i] = [sin(theta), 0.0, cos(theta)]
    end
    for t in Temp:-0.1:1e-2
        for epoch in 1:1000
            this_node = rand(Vector(1:length(sp.onsite)))
            sp_sa_step!(sp, t, this_node)
        end
    end
    return sp_energy(sp, 1, 1, zeros(length(sp.onsite)))
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

function sp_energy(sp::SpinGlassModel, t, T, Vtrans)
    ret = 0.0
    ret += sum([t[3] for t in sp.M] .* sp.onsite) * (t/T)
    ret += sum([t[1] for t in sp.M] .* Vtrans) * ((T-t)/T)
    for edge in sp.edges
        i, j, w = edge
        ret += w * sp.M[i][3] * sp.M[j][3] * (t/T)
    end
    return ret
end

function sp_check_vaild(rule::CellularAutomata1D, sp::SpinGlassModel)
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


function sp_check_vaild_time(rule::CellularAutomata1D, sp::SpinGlassModel, T, Vtrans)
    @info "begin running runge_kutta"
    runge_kutta_integrate!(sp, min(1e-2, T/1e5), T, Vtrans)
    # @info "sp.M = $(sp.M)"
    for i in 1:length(sp.onsite)
        if 1 - abs(sp.M[i][3]) > 1e-1
            return false
        end
    end
    @info "successfully turn into final hamiltonian"
    return sp_check_vaild(rule, sp)
end

function vectorgradient(vec, this_t)
    n = Int(vec[end-2])
    m = Int(vec[end-1])
    gradient = Float64(vec[end])
    sp = spinglass_mapping(n, m;gradient=gradient)
    @assert length(vec) == 3*length(sp.onsite)+1+length(sp.onsite)+3
    sp.M .= [Point(vec[i], vec[i+1], vec[i+2]) for i in 1:3:3*length(sp.onsite)]
    Vtrans = Vector(vec[3*length(sp.onsite)+2:3*length(sp.onsite)+2+length(sp.onsite)-1])
    Tmax = Float64(vec[3*length(sp.onsite)+1])
    Mxdot, Mydot, Mzdot = integrator(sp, this_t, Tmax, Vtrans)
    # @assert Mxdot[1] == 0.0
    # @assert Mydot[1] == 0.0
    # @assert Mzdot[1] == 0.0
    ret = Vector{Float64}()
    for i in 1:length(sp.onsite)
        append!(ret, [Mxdot[i], Mydot[i], Mzdot[i]])
    end
    append!(ret, zeros(length(vec) - 3*length(sp.onsite)))
    return ret
end

function vector2sp(vec)
    n = Int(vec[end-2])
    m = Int(vec[end-1])
    gradient = Float64(vec[end])
    sp = spinglass_mapping(n, m;gradient=gradient)
    @assert length(vec) == 3*length(sp.onsite)+1+length(sp.onsite)+3
    sp.M .= [Point(vec[i], vec[i+1], vec[i+2]) for i in 1:3:3*length(sp.onsite)]
    return sp
end

function initialvector(Tmax, n, m; gradient = 1.0)
    sp = spinglass_mapping(n, m; gradient=gradient)
    Vtrans = fill(1.0, length(sp.onsite))
    ret = Vector{Float64}()
    for i in 1:length(sp.onsite)
        append!(ret, [sp.M[i][1], sp.M[i][2], sp.M[i][3]])
    end
    push!(ret, Float64(Tmax))
    append!(ret, Vtrans)
    push!(ret, n)
    push!(ret, m)
    push!(ret, gradient)
    return ret
end

function fcn(x, y, f)
    g = vectorgradient(y, x)
    copyto!(f, g)
end

function spingls!(du, u, p, t)
    g = vectorgradient(u, t)
    copyto!(du, g)
end

function printsp(sp::SpinGlassModel)
    for i in 1:length(sp.onsite)
        for j in 3:3
            print(sp.M[i][j]," ")
        end
        println("")
    end
end

function spinglass_adiabatic_dp8(n, m, Tmax, gradient)
    init_dp8 = initialvector(Tmax, n, m; gradient=gradient)
    solver8 = DP8Solver(fcn, 0.0, init_dp8; atol=1e-10, rtol=1e-10, maximum_allowed_steps=5000000)
    @time integrate!(solver8, Tmax)
    sp_ret = vector2sp(get_current_state(solver8))
    return sp_ret
end