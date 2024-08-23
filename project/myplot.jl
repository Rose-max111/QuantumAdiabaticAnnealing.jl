using CairoMakie
using CurveFit

function draw(xdata, ydata)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "loglog(λ)+(1/log(λ))", ylabel = "log(sweep times)",
        title = "Time v.s. Gradient(50% success probability)")
    scatter!(ax, xdata, ydata)

    fit = curve_fit(LinearFit, xdata, ydata)
    lines!(ax, xdata, fit.(xdata))
    # ylims!(ax, low=0)
    f
end

function draw_timevsdepth(xdata, ydata)
    f = Figure()
    ax = Axis(f[1, 1], xlabel = "depth * log^2(depth)", ylabel = "sweep times",
        title = "Time v.s. Depth(50% success probability, Λ=1.3, width=15) -- Energy Gradient")
    scatter!(ax, xdata, ydata)

    # fit = curve_fit(Polynomial, xdata, ydata, 2)
    fit = curve_fit(LinearFit, xdata, ydata)
    lines!(ax, xdata, fit.(xdata))
    # ylims!(ax, low=0)
    f
end

# draw(xx[1:10], yy[1:10])

function __main__(graph_width, graph_depth, λ, GW)
    evaluate_time = Vector{Float64}()
    for d in graph_depth
        filepath = joinpath(@__DIR__, "data_toymodel_pulse/W=$(graph_width)_D=$(d)_GW=$(GW)_E=$(λ).txt")
        # filepath = joinpath(@__DIR__, "data_toymodel_pulse/W=$(graph_width)_D=$(d)_E=$(λ).txt")
        open(filepath, "r") do file
            push!(evaluate_time, parse(Float64, readline(file)))
        end
    end
    graph_depth = Float64.(graph_depth)
    draw_timevsdepth(((graph_depth) .* log.(graph_depth) .* log.(graph_depth))[1:end], (evaluate_time[1:end]))

    # filepath = joinpath(@__DIR__, "test.txt")
    # open(filepath, "w") do file
    #     for i in graph_depth
    #         println(file, i)
    #     end
    #     for i in evaluate_time
    #         println(file, i)
    #     end
    # end
end

function __main__gradient(graph_width, graph_depth, λ)
    evaluate_time = Vector{Float64}()
    for d in graph_depth
        filepath = joinpath(@__DIR__, "data_toymodel/W=$(graph_width)_D=$(d)_E=$(λ).txt")
        # filepath = joinpath(@__DIR__, "data_toymodel_pulse/W=$(graph_width)_D=$(d)_E=$(λ).txt")
        open(filepath, "r") do file
            push!(evaluate_time, parse(Float64, readline(file)))
        end
    end
    graph_depth = Float64.(graph_depth)
    draw_timevsdepth(((graph_depth) .* log.(graph_depth) .* log.(graph_depth))[1:end], (evaluate_time[1:end]))

    # filepath = joinpath(@__DIR__, "test.txt")
    # open(filepath, "w") do file
    #     for i in graph_depth
    #         println(file, i)
    #     end
    #     for i in evaluate_time
    #         println(file, i)
    #     end
    # end
end


# __main__(12, [6, 8, 10, 12, 15, 18, 22, 25, 28, 31, 34, 36, 39, 42, 45, 48, 50, 53, 56, 58, 61, 64, 67, 69, 72, 74, 77, 80, 82, 85, 88, 90, 92], 1.5, 1.0)
# __main__(15, [4, 6, 8, 10, 12, 15, 18, 20, 22, 25, 28, 30, 33, 35, 37], 1.3, 1.0)

__main__gradient(15, [4, 6, 8, 10, 12, 15, 18, 20, 22, 25, 28, 30, 33, 35, 37], 1.3)