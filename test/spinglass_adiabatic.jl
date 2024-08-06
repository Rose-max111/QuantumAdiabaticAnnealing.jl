using Test
using DormandPrince
using DifferentialEquations
using QuantumAdiabaticAnnealing
using QuantumAdiabaticAnnealing: initialvector, spinglass_mapping, vector2sp, spingls!, runge_kutta_integrate!, euclidean_integrate!, fcn, freeze_input!
using QuantumAdiabaticAnnealing: printsp, sp_energy
using QuantumAdiabaticAnnealing: instantaneous_field, instantaneous_field_autodiff

@testset "integrater_double_check" begin
    Tmax = 1000.0
    init_dp5 = initialvector(Tmax, 3, 4;gradient=1.0)
    init_dp8 = initialvector(Tmax, 3, 4;gradient=1.0)
    solver5 = DP5Solver(fcn, 0.0, init_dp5; atol=1e-12, rtol=1e-12, maximum_allowed_steps=20000000)
    @time integrate!(solver5, Tmax)
    solver8 = DP8Solver(fcn, 0.0, init_dp8; atol=1e-12, rtol=1e-12, maximum_allowed_steps=5000000)
    @time integrate!(solver8, Tmax)

    init_de = initialvector(Tmax, 3, 4;gradient=1.0)
    # cb = ManifoldProjection(gproject)
    tspan = (0.0, Tmax)
    prob = ODEProblem(spingls!, init_de, tspan)
    @time sol = DifferentialEquations.solve(prob, Vern7(), reltol = 1e-12, abstol=1e-12)

    sp_this5 = vector2sp(get_current_state(solver5))
    sp_this8 = vector2sp(get_current_state(solver8))
    sp_de = vector2sp(sol[end])

    # for i in 1:length(sp_this5.onsite)
    #     for j in 1:3
    #         @test abs(sp_this5.M[i][j] - sp_de.M[i][j]) <= 1e-1
    #         @test abs(sp_this8.M[i][j] - sp_de.M[i][j]) <= 1e-1
    #         @test abs(sp_this5.M[i][j] - sp_this8.M[i][j]) <= 1e-1
    #     end
    # end

    @info "dp5 solution"
    printsp(sp_this5)
    println(sp_energy(sp_this5, Tmax, Tmax, fill(1.0, length(sp_this5.onsite))))
    @info "dp8 solution"
    printsp(sp_this8)
    println(sp_energy(sp_this8, Tmax, Tmax, fill(1.0, length(sp_this5.onsite))))
    @info "de solution"
    printsp(sp_de)
    println(sp_energy(sp_de, Tmax, Tmax, fill(1.0, length(sp_this5.onsite))))
end

@testset "runge_kutta_vs_integrator" begin
    Tmax = 100.0
    n=3
    m=3
    gradient = 1.0
    init_dp8 = initialvector(Tmax, n, m;gradient=gradient)
    solver8 = DP8Solver(fcn, 0.0, init_dp8; atol=1e-10, rtol=1e-10, maximum_allowed_steps=5000000)
    @time integrate!(solver8, Tmax)
    sp_dp8 = vector2sp(get_current_state(solver8))

    sp_rk = spinglass_mapping(n, m; gradient=gradient)
    freeze_input!(sp_rk)
    Vtrans = fill(1.0, length(sp_rk.onsite))
    @time runge_kutta_integrate!(sp_rk, 1e-4, Tmax, Vtrans; T_end = Tmax)

    for i in 1:length(sp_rk.onsite)
        for j in 1:3
            @test abs(sp_rk.M[i][j] - sp_dp8.M[i][j]) <= 1e-3
        end
    end
end

@testset "runge_kutta_error" begin
    sp_r1 = spinglass_mapping(4, 2)
    sp_r2 = spinglass_mapping(4, 2)
    Vtrans = fill(1.0, length(sp_r1.onsite))
    runge_kutta_integrate!(sp_r1, 1e-4, 100.0, Vtrans; T_end = 100.0)
    runge_kutta_integrate!(sp_r2, 1e-3, 100.0, Vtrans; T_end = 100.0)

    for i in 1:length(sp_r1.onsite)
        for j in 1:3
            @test abs(sp_r1.M[i][j] - sp_r2.M[i][j]) <= 1e-5
        end
    end
end

@testset "runge_kutta_with_euclidean" begin
    sp_r = spinglass_mapping(3, 2)
    sp_e = spinglass_mapping(3, 2)
    Vtrans = fill(1.0, length(sp_r.onsite))
    runge_kutta_integrate!(sp_r, 1e-2, 10.0, Vtrans)
    euclidean_integrate!(sp_e, 1e-6, 10.0, Vtrans)

    for i in 1:length(sp_r.onsite)
        for j in 1:3
            @test abs(sp_r.M[i][j] - sp_e.M[i][j]) <= 1e-2
        end
    end
end

function runge_kutta_visualize()
    T_end = 1000.0
    n=3
    m=4
    sp = spinglass_mapping(n, m)

    # for i in sp.n+1:length(sp.onsite)
    #     aa = rand()*0.2 - 0.1
    #     bb = rand()*0.2 - 0.1
    #     # @info "aa = $aa, bb=$bb, i=$i"
    #     sp.M[i][3] = bb
    #     sp.M[i][2] = aa
    #     sp.M[i][1] = -sqrt(1 - bb^2 - aa^2)
    #     # @info "i=$i, sp.M[i]= $(sp.M[i])"
    # end

    Vtrans = fill(1.0, length(sp.onsite))
    model_print, field_print = runge_kutta_integrate!(sp, 1e-2, T_end, Vtrans; T_end = T_end)
    @info "energy is $(sp_energy(sp, T_end, T_end, Vtrans))"
    @info "using autodiff energy is $(spinglass_hamiltonian(sp, T_end, T_end, Vtrans))"

    open("rk.txt","w") do io
        for i in 1:length(model_print)
            for j in model_print[i]
                print(io, j, " ")
            end
            println(io, "")
        end
        for i in 1:length(field_print)
            for j in field_print[i]
                print(io, j, " ")
            end
            println(io, "")
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

@testset "basic_spinglass" begin
    n=3
    m=5
    sp = spinglass_mapping(3, 5)
    for edge in sp.edges
        i, j, w = edge[1], edge[2], edge[3]
        @test i <= n*m || j<=n*m # no ancilla connection
    end
    sp = spinglass_mappint(1, 1)
    sp_gs = sp_energy(sp)
end

using Random
@testset "autodiff vs exact" begin
    sp=spinglass_mapping(3, 4)
    for i in 1:length(sp.onsite)
        sp.M[i]=[rand()*2-1, rand()*2-1, rand()*2-1]
    end
    @info "sp.M[1] = $(sp.M[1])"
    Vtrans = fill(1.0, length(sp.onsite))
    Tmax = 10.0
    this_t = 3.0
    H_autodiff = instantaneous_field_autodiff(sp, this_t, Tmax, Vtrans)
    H_exact = instantaneous_field(sp, this_t, Tmax, Vtrans)
    
    for i in 1:length(sp.onsite)
        for j in 1:3
            @test abs(H_exact[i][j] - H_autodiff[i][j]) <= 1e-6
        end
    end
end

@testset "singlemodel" begin
    single_onsite = [1.0, 2.0, 2.0, 2.0, 5.0] # last one is ancilla
    single_edge_weights = [1.0, 1.0, 2.0, 3.0, 2.0, 2.0, 5.0, 2.0, 5.0, 6.0]
    configs = [1.0, 1.0, 1.0, -1.0, -1.0]
    ret = 0.0
    cnt = 0
    for i in 1:5
        for j in i+1:5
            cnt += 1
            ret += single_edge_weights[cnt] * configs[i] * configs[j]
        end
    end
    ret += sum(configs .* single_onsite)
    @test ret == -11
end