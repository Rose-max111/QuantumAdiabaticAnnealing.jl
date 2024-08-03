using Test
using DormandPrince
using DifferentialEquations
using QuantumAdiabaticAnnealing
using QuantumAdiabaticAnnealing: initialvector, spinglass_mapping, vector2sp, spingls!, runge_kutta_integrate!, euclidean_integrate!, fcn, freeze_input!

@testset "integrater_double_check" begin
    Tmax = 2000.0
    init_dp5 = initialvector(Tmax, 3, 3)
    init_dp8 = initialvector(Tmax, 3, 3)
    solver5 = DP5Solver(fcn, 0.0, init_dp5; atol=1e-10, rtol=1e-10, maximum_allowed_steps=5000000)
    @time integrate!(solver5, Tmax)
    solver8 = DP8Solver(fcn, 0.0, init_dp8; atol=1e-10, rtol=1e-10, maximum_allowed_steps=5000000)
    @time integrate!(solver8, Tmax)

    init_de = initialvector(Tmax, 3, 3)
    # cb = ManifoldProjection(gproject)
    tspan = (0.0, Tmax)
    prob = ODEProblem(spingls!, init_de, tspan)
    @time sol = DifferentialEquations.solve(prob, Vern7(), reltol = 1e-10, abstol=1e-10)

    sp_this5 = vector2sp(get_current_state(solver5))
    sp_this8 = vector2sp(get_current_state(solver8))
    sp_de = vector2sp(sol[end])

    for i in 1:length(sp_this5.onsite)
        for j in 1:3
            @test abs(sp_this5.M[i][j] - sp_de.M[i][j]) <= 1e-4
            @test abs(sp_this8.M[i][j] - sp_de.M[i][j]) <= 1e-4
            @test abs(sp_this5.M[i][j] - sp_this8.M[i][j]) <= 1e-4
        end
    end
end

@testset "runge_kutta_vs_integrator" begin
    Tmax = 123.0
    n=4
    m=4
    gradient = 1.1
    init_dp8 = initialvector(Tmax, n, m;gradient=1.1)
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
