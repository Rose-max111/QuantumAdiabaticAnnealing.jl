mutable struct spinglassmodel
    # n::Int
    const n::Int
    const m::Int
    const edges::Vector{Tuple{Int,Int,Float64}}
    const onsite::Vector{Float64}
    M::Vector{Vector{Float64}}
    function spinglassmodel(n::Int, m::Int, edges::Vector{Tuple{Int,Int,Float64}}, onsite::Vector{Float64}, M::Vector{Vector{Float64}} = fill([-1.0, 0.0, 0.0], length(onsite)))
        new(n, m, edges, onsite, M)
    end
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
    return spinglassmodel(n, m, edges, onsites)
end

function freeze_input!(sp::spinglassmodel)
    for i in 1:sp.n
        sp.onsite[i] += 50.0
    end
end


# H(i) = -hx̂ + ∑_j J_{i,j} M_j^z ẑ + onsite ẑ
# ̇M = H \cross M
function crossdot(A::Vector{Float64}, B::Vector{Float64})
    return [A[2]*B[3] - A[3]*B[2], A[3]*B[1] - A[1]*B[3], A[1]*B[2] - A[2]*B[1]]
end

function integrator(sp::spinglassmodel, H::Vector{Vector{Float64}})
    Mdot = crossdot.(H, sp.M)
    return Mdot
end

function instantaneous_field(sp::spinglassmodel, t, T, Vtrans::Vector{Float64})
    H = [[-Vtrans[i] .* ((T-t)/T), 0.0, sp.onsite[i] .* (t/T)] for i in 1:length(sp.onsite)]
    for edge in sp.edges
        i, j, w = edge
        H[i][3] += w * sp.M[j][3] * (t/T)
        H[j][3] += w * sp.M[i][3] * (t/T)
    end
    return H
end

function integrator(sp::spinglassmodel, t, T, Vtrans::Vector{Float64}) # evaluate F(t, y)
    Mxdot = - sp.onsite .* [u[2] for u in sp.M]
    Mydot = sp.onsite .* [u[1] for u in sp.M]
    for edge in sp.edges
        i, j, w = edge
        Mydot[i] += w * sp.M[j][3] * sp.M[i][1]
        Mydot[j] += w * sp.M[i][3] * sp.M[j][1]

        Mxdot[i] -= w * sp.M[j][3] * sp.M[i][2]
        Mxdot[j] -= w * sp.M[i][3] * sp.M[j][2]
    end
    Mxdot .*= (t/T)
    Mydot .*= (t/T)
    Mydot -= - Vtrans .* [u[3] for u in sp.M] .* ((T-t)/T)
    Mzdot = - Vtrans .* [u[2] for u in sp.M] .* ((T-t)/T)

    return [Mxdot, Mydot, Mzdot]
end

function runge_kutta_singlejump!(sp::spinglassmodel, t0, delta_t, T, Vtrans::Vector{Float64})
    origin_M = copy(sp.M)
    k1 = integrator(sp, t0, T, Vtrans)
    sp.M = [renorm([origin_M[i][1] + k1[1][i] * delta_t / 2, origin_M[i][2] + k1[2][i] * delta_t / 2, origin_M[i][3] + k1[3][i] * delta_t / 2]) for i in 1:length(sp.onsite)]
    k2 = integrator(sp, t0 + delta_t / 2, T, Vtrans)
    sp.M = [renorm([origin_M[i][1] + k2[1][i] * delta_t / 2, origin_M[i][2] + k2[2][i] * delta_t / 2, origin_M[i][3] + k2[3][i] * delta_t / 2]) for i in 1:length(sp.onsite)]
    k3 = integrator(sp, t0 + delta_t / 2, T, Vtrans)
    sp.M = [renorm([origin_M[i][1] + k3[1][i] * delta_t, origin_M[i][2] + k3[2][i] * delta_t, origin_M[i][3] + k3[3][i] * delta_t]) for i in 1:length(sp.onsite)]
    k4 = integrator(sp, t0 + delta_t, T, Vtrans)
    real_k = (k1 .+ 2*k2 .+ 2*k3 .+ k4) ./ 6 * delta_t
    sp.M = [renorm([origin_M[i][1] .+ real_k[1][i], origin_M[i][2]+real_k[2][i], origin_M[i][3]+real_k[3][i]]) for i in 1:length(sp.onsite)]
end

function runge_kutta_integrate!(sp::spinglassmodel, dt::Float64, T::Float64, Vtrans::Vector{Float64}; T_end = nothing)
    t0 = 0.0
    if T_end == nothing
        T_end = T
    end
    sp_error_moniter = spinglass_mapping(sp.n, sp.m)
    total_error = 0.0
    while t0 < T_end
        delta_t = min(dt, T_end-t0)
        # now calculate jump with 2-step
        # sp_error_moniter.M = copy(sp.M)
        # runge_kutta_singlejump!(sp_error_moniter, t0, delta_t/2, T, Vtrans)
        # runge_kutta_singlejump!(sp_error_moniter, t0, delta_t/2, T, Vtrans)
        # now calculate jump with 1-step

        # inst_H = instantaneous_field(sp, t0, T, Vtrans)
        # inst_norm = [sqrt(sum(inst_H[i] .* inst_H[i])) for i in 1:length(sp.onsite)]
        # inst_dot_norm = [sum(inst_H[i] .* sp.M[i]) for i in 1:length(sp.onsite)]
        # inst_cos_theta = inst_dot_norm ./ inst_norm
        # inst_min_cos_theta = minimum(inst_cos_theta)
        # @info "t = $t0, inst_min_cos_theta = $inst_min_cos_theta"
        # if inst_min_cos_theta < 0
        #     @info "inst_H = $inst_H"
        #     @info "sp.M = $(sp.M)"
        #     error("err")
        # end

        runge_kutta_singlejump!(sp, t0, delta_t, T, Vtrans)
        # now calculate local truncation error
        # local_error = maximum(sp.M .- sp_error_moniter.M) / (1 - 2^-4)
        # total_error += maximum(local_error)
        # @assert maximum(local_error) <= 1e-10
        # @info "t = $t0, local_error = $local_error"
        t0 += delta_t
    end
    # @info "total_error = $total_error"
end

function euclidean_integrate!(sp::spinglassmodel, dt, T, Vtrans::Vector{Float64})
    t=0
    while t<T
        delta_t = min(dt, T-t)
        Mxdot, Mydot, Mzdot = integrator(sp, t, T, Vtrans)
        d_Mx = Mxdot .* delta_t
        d_My = Mydot .* delta_t
        d_Mz = Mzdot .* delta_t
        sp.M = [renorm([sp.M[i][1]+d_Mx[i], sp.M[i][2]+d_My[i], sp.M[i][3]+d_Mz[i]]) for i in 1:length(sp.onsite)]
        t += delta_t
        # @info "t = $t, energy = $(sp_energy(sp, t, T, Vtrans))"
    end
end

function renorm(point)
    norm = sqrt(sum(point .* point))
    return point ./ norm
end

function sp_groud_state(sp::spinglassmodel)
    hyperedges = [[t[1],t[2]] for t in sp.edges]
    hyperweights = [t[3] for t in sp.edges]
    for i in 1:length(sp.onsite)
        push!(hyperedges, [i])
        push!(hyperweights, sp.onsite[i])
    end
    spgls = SpinGlass(length(sp.onsite), hyperedges, hyperweights)
    spproblem = GenericTensorNetwork(spgls)
    gs = solve(spproblem, SizeMin())[]
    return gs
end

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

function sp_check_vaild(sp::spinglassmodel)
    for i in 1:length(sp.onsite)
        if 1 - abs(sp.M[i][3]) > 1e-1
            return false
        end
    end
    state = [sp.M[i][3] > 0 ? 1 : 0 for i in 1:length(sp.onsite)]
    for j in 1:sp.m-1
        for i in 1:sp.n
            pre_l = cartesian_to_linear(mod1(i-1, sp.n), j, sp.n)
            mid_l = cartesian_to_linear(i, j, sp.n)
            suf_l = cartesian_to_linear(mod1(i+1, sp.n), j, sp.n)
            @info "i = $i, j = $j"
            if rule110(state[pre_l], state[mid_l], state[suf_l]) != state[cartesian_to_linear(i, j+1, sp.n)]
                return false
            end
        end
    end
    return true
end


function sp_check_vaild_time(sp::spinglassmodel, T, Vtrans)
    @info "begin running runge_kutta"
    runge_kutta_integrate!(sp, min(1e-2, T/1e5), T, Vtrans)
    # @info "sp.M = $(sp.M)"
    for i in 1:length(sp.onsite)
        if 1 - abs(sp.M[i][3]) > 1e-1
            return false
        end
    end
    @info "successfully turn into final hamiltonian"
    return sp_check_vaild(sp)
end

function vectorgradient(vec, this_t; freeze=false)
    n = Int(vec[end-2])
    m = Int(vec[end-1])
    gradient = Float64(vec[end])
    sp = spinglass_mapping(n, m;gradient=gradient)
    @assert length(vec) == 3*length(sp.onsite)+1+length(sp.onsite)+3
    if freeze == true
        freeze_input!(sp)
    end
    sp.M = [Vector(vec[i:i+2]) for i in 1:3:3*length(sp.onsite)]
    Vtrans = Vector(vec[3*length(sp.onsite)+2:3*length(sp.onsite)+2+length(sp.onsite)-1])
    Tmax = Float64(vec[3*length(sp.onsite)+1])
    Mxdot, Mydot, Mzdot = integrator(sp, this_t, Tmax, Vtrans)
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
    sp.M = [vec[i:i+2] for i in 1:3:3*length(sp.onsite)]
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
    g = vectorgradient(y, x;freeze=true)
    for i in 1:length(g)
        f[i] = g[i]
    end
end

function spingls!(du, u, p, t)
    g = vectorgradient(u, t;freeze=true)
    for i in 1:length(g)
        du[i] = g[i]
    end
end

function printsp(sp::spinglassmodel)
    for i in 1:length(sp.onsite)
        for j in 1:3
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