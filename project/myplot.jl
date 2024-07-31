using CairoMakie
using CurveFit

graph_width=12
# graph_depth=10

evaluate_time = Vector{Float64}()
# gradient = append!([1.2, 1.22, 1.23, 1.25, 1.27, 1.28, 1.29], Vector(1.3:0.1:3.0))
# for λ in gradient
#     filepath = joinpath(@__DIR__, "data_toymodel/lambda_sweep/W=$(graph_width)_D=$(graph_depth)_E=$(λ).txt")
#     open(filepath, "r") do file
#         push!(evaluate_time, parse(Float32, readline(file)))
#     end
# end

graph_depth = [6, 8, 10, 12, 15, 18, 22, 25, 28, 31, 34, 36]
λ = 1.5
for d in graph_depth
    filepath = joinpath(@__DIR__, "data_toymodel/W=$(graph_width)_D=$(d)_E=$(λ).txt")
    open(filepath, "r") do file
        push!(evaluate_time, parse(Float64, readline(file)))
    end
end
# xx = (log.(log.(gradient)) .+ 1.0 ./ log.(gradient))
# yy = log.(evaluate_time)

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
    ax = Axis(f[1, 1], xlabel = "depth", ylabel = "sweep times",
        title = "Time v.s. Depth(50% success probability, λ=1.5, width=12)")
    scatter!(ax, xdata, ydata)

    # fit = curve_fit(Polynomial, xdata, ydata, 2)
    fit = curve_fit(LinearFit, xdata, ydata)
    lines!(ax, xdata, fit.(xdata))
    # ylims!(ax, low=0)
    f
end

draw_timevsdepth(Float64.(graph_depth)[1:end], (evaluate_time[1:end]))
# draw(xx[1:10], yy[1:10])

filepath = joinpath(@__DIR__, "test.txt")
open(filepath, "w") do file
    for i in graph_depth
        println(file, i)
    end
    for i in evaluate_time
        println(file, i)
    end
end