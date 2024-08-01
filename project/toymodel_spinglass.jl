using GenericTensorNetworks
# using CairoMakie
# using GLMakie
using Test

using DormandPrince
using DifferentialEquations

function rule110(p, q, r)
    return (q + r + q*r + p*q*r) % 2
end


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
            # @info "this_id = $this_id"
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
function crossdot(A::Vector{Float64}, B::Vector{Float64})
    return [A[2]*B[3] - A[3]*B[2], A[3]*B[1] - A[1]*B[3], A[1]*B[2] - A[2]*B[1]]
end

function integrator(sp::spinglassmodel, H::Vector{Vector{Float64}})
    Mdot = crossdot.(H, sp.M)
    return Mdot
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

    # H = [[-Vtrans[i] .* ((T-t)/T), 0.0, sp.onsite[i] .* (t/T)] for i in 1:length(sp.onsite)]
    # for edge in sp.edges
    #     i, j, w = edge
    #     H[i][3] += w * sp.M[j][3] * (t/T)
    #     H[j][3] += w * sp.M[i][3] * (t/T)
    # end
    # Mdot_bf = integrator(sp, H)
    # for i in 1:length(sp.onsite)
    #     # @info "i = $i, Mdot = $(Mdot_bf[i]), Mxdot = $(Mxdot[i]), Mydot = $(Mydot[i]), Mzdot = $(Mzdot[i])"
    #     @assert abs(Mxdot[i] - Mdot_bf[i][1]) <= 1e-5
    #     @assert abs(Mydot[i] - Mdot_bf[i][2]) <= 1e-5
    #     @assert abs(Mzdot[i] - Mdot_bf[i][3]) <= 1e-5
    # end

    return [Mxdot, Mydot, Mzdot]
end

function runge_kutta_integrate!(sp::spinglassmodel, dt, T, Vtrans::Vector{Float64})
    t0 = 0
    while t0 < T
        delta_t = min(dt, T-t0)
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
        t0 += delta_t
    end
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

@testset "runge_kutta_with_euclidean" begin
    sp_r = spinglass_mapping(3, 2)
    sp_e = spinglass_mapping(3, 2)
    Vtrans = fill(1.0, length(sp_r.onsite))
    runge_kutta_integrate!(sp_r, 1e-2, 10, Vtrans)
    euclidean_integrate!(sp_e, 1e-6, 10, Vtrans)

    for i in 1:length(sp_r.onsite)
        for j in 1:3
            @test abs(sp_r.M[i][j] - sp_e.M[i][j]) <= 1e-2
        end
    end
end

@testset "runge_kutta_error" begin
    sp_r1 = spinglass_mapping(3, 3)
    sp_r2 = spinglass_mapping(3, 3)
    Vtrans = fill(1.0, length(sp_r1.onsite))
    runge_kutta_integrate!(sp_r1, 1e-5, 60, Vtrans)
    runge_kutta_integrate!(sp_r2, 1e-4, 60, Vtrans)
    for i in 1:length(sp_r1.onsite)
        @info "sp_r1.M = $(sp_r1.M[i]), sp_r2.M = $(sp_r2.M[i])"
        @info "delta = $(sp_r1.M[i] .- sp_r2.M[i])"
    end

    for i in 1:length(sp_r1.onsite)
        for j in 1:3
            @test abs(sp_r1.M[i][j] - sp_r2.M[i][j]) <= 1e-2
        end
    end
end

function check_vaild_time(sp::spinglassmodel, T, Vtrans)
    @info "begin running runge_kutta"
    runge_kutta_integrate!(sp, min(1e-2, T/1e5), T, Vtrans)
    # @info "sp.M = $(sp.M)"
    for i in 1:length(sp.onsite)
        if 1 - abs(sp.M[i][3]) > 1e-1
            return false
        end
        # @test 1 - abs(sp.M[i][3]) <= 1e-1
    end
    @info "successfully turn into final hamiltonian"
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

function evaluate_adiabatic_time(n, m)
    Vtrans = fill(1.0, n*m + n*(m-1))
    adiabatic_time = 0
    max_try = 1
    while max_try > 0
        next_try = max_try
        sp = spinglass_mapping(n, m)
        @info "Stage1, max_try = $max_try, next_try = $next_try, checking if vaild"
        flag = check_vaild_time(sp, next_try, Vtrans)
        if flag == true
            break
        end
        max_try *= 2
    end
    if max_try == 1
        return 1
    end
    adiabatic_time += max_try / 2
    max_try /= 2
    while max_try != 1
        max_try /= 2
        next_try = adiabatic_time + max_try
        sp = spinglass_mapping(n, m)
        @info "Stage2, max_try = $max_try, next_try = $next_try, checking if vaild"
        flag = check_vaild_time(sp, next_try, Vtrans)
        if flag == false
            adiabatic_time = next_try
        end
    end
    return adiabatic_time + 1
end

# evaluate_adiabatic_time(3, 3)
# check_vaild_time(spinglass_mapping(3, 3), 1, fill(1.0, 3*3+3*2))



function vectorgradient(vec, this_t, n, m)
    sp = spinglass_mapping(n, m)
    sp.M = [Vector(vec[i:i+2]) for i in 1:3:3*length(sp.onsite)]
    Vtrans = Vector(vec[end-length(sp.onsite)+1:end])
    Tmax = Float64(vec[3*length(sp.onsite) + 1])
    Mxdot, Mydot, Mzdot = integrator(sp, this_t, Tmax, Vtrans)
    ret = Vector{Float64}()
    for i in 1:length(sp.onsite)
        append!(ret, [Mxdot[i], Mydot[i], Mzdot[i]])
    end
    append!(ret, zeros(length(sp.onsite) + 1))
    return ret
end

function vector2sp(vec, n, m)
    sp = spinglass_mapping(n, m)
    sp.M = [vec[i:i+2] for i in 1:3:3*length(sp.onsite)]
    return sp
end

function initialvector(Tmax, n, m)
    sp = spinglass_mapping(n, m)
    Vtrans = fill(1.0, length(sp.onsite))
    ret = Vector{Float64}()
    for i in 1:length(sp.onsite)
        append!(ret, [sp.M[i][1], sp.M[i][2], sp.M[i][3]])
    end
    push!(ret, Float64(Tmax))
    append!(ret, Vtrans)
    return ret
end

function fcn(x, y, f)
    g = vectorgradient(y, x, 3, 3)
    for i in 1:length(g)
        f[i] = g[i]
    end
end

function spingls!(du, u, p, t)
    g = vectorgradient(u, t, 3, 3)
    for i in 1:length(g)
        du[i] = g[i]
    end
end

function gproject(resid, u, p, t)
    for i in 1:3:3*(3*3+3*2)
        resid[i] = 1 - (u[i]^2 + u[i+1]^2 + u[i+2]^2)
        resid[i+1] = 0
        resid[i+2] = 0
    end
end

@testset "runge_kutta_vs_dmprince" begin
    Tmax = 100.0
    init_dp5 = initialvector(Tmax, 3, 3)
    init_dp8 = initialvector(Tmax, 3, 3)
    solver5 = DP5Solver(fcn, 0.0, init_dp5; atol=1e-12, rtol=1e-12, maximum_allowed_steps=500000)
    integrate!(solver5, Tmax)
    solver8 = DP8Solver(fcn, 0.0, init_dp8; atol=1e-12, rtol=1e-12, maximum_allowed_steps=500000)
    integrate!(solver8, Tmax)

    init_de = initialvector(Tmax, 3, 3)
    cb = ManifoldProjection(gproject)
    tspan = (0.0, Tmax)
    prob = ODEProblem(spingls!, init_de, tspan)
    sol = DifferentialEquations.solve(prob, Vern7(), reltol = 1e-12, abstol=1e-12)

    sp_this5 = vector2sp(get_current_state(solver5), 3, 3)
    sp_this8 = vector2sp(get_current_state(solver8), 3, 3)
    sp_de = vector2sp(sol[end], 3, 3)
    for i in 1:length(sp_this5.onsite)
        for j in 1:3
            @test abs(sp_this5.M[i][j] - sp_de.M[i][j]) <= 1
            # @test abs(sp_this8.M[i][j] - sp_de.M[i][j]) <= 1e-2
            # @test abs(sp_this5.M[i][j] - sp_this8.M[i][j]) <= 1
        end
    end
    sp_r = spinglass_mapping(3, 3)
    Vtrans = fill(1.0, length(sp_r.onsite))
    # runge_kutta_integrate!(sp_r, 1e-4, Tmax, Vtrans)
    for i in 1:length(sp_r.onsite)
        @test abs(1 - sum(sp_de.M[i] .* sp_de.M[i])) <= 1e-5
        @test abs(1 - sum(sp_this5.M[i] .* sp_this5.M[i])) <= 1e-5
        @test abs(1 - sum(sp_this8.M[i] .* sp_this8.M[i])) <= 1e-5
        @info "i = $i, sp_de.M = $(sp_de.M[i])"
        @info "i = $i, sp_this5.M = $(sp_this5.M[i])"
        for j in 1:3
            # @test abs(sp_r.M[i][j] - sp_de.M[i][j]) <= 1e-3
        end
    end
end

@testset "simple_magneticfield" begin
    # ss = spinglassmodel(1, 1, Vector{Tuple{Int64,Int64,Float64}}(), [2.0], [renorm([0.0, 0.0, 1.0])])
    # tt = 0
    # # TT = round(2*π, digits=4)
    # TT = 50
    # pnt = Vector{Point{3, Float64}}()
    # while tt < TT
    #     delta_t = min(1e-3, TT-tt)
    #     Mdot = integrator(ss, [[sin((tt / TT) * π / 2), 0.0, cos((tt / TT) * π / 2)]])
    #     # Mdot = integrator(ss, [[0.0, 0.0, 1.0]])
    #     ss.M = renorm.(ss.M .+ (Mdot .* delta_t))
    #     push!(pnt, Point((ss.M[1][1], ss.M[1][2], ss.M[1][3])))
    #     tt += delta_t
    # end
    # println(ss.M)
    # f = Figure(size=(900,900))
    # # lscene = LScene(f[1,1])
    # # axislimits!(lscene, Rect3f(Vec3f(-1, -1, -1), Vec3f(1, 1, 1)))
    # ax1 = Axis3(f[1,1], aspect = (2, 2, 2))
    # ax2 = Axis3(f[1,2], aspect = (2, 2, 2), elevation = (0.5 * π))

    # scatter!(ax1, pnt, color=:blue, fxaa=true)
    # scatter!(ax2, pnt, color=:blue, fxaa=true)
    # # xlims!(-1, 1)
    # # ylims!(f, -1, 1)
    # # zlims!(f, -1, 1)
    # save("1.png", f)
end
