using GLMakie

function renorm(vec)
    return vec ./ sqrt(sum(vec.^2))
end

function init_data(n, m)
    Tstep=1000 / 1e-2 + 1
    repoint = Vector{Vector{Vector{Float64}}}()
    refield = Vector{Vector{Vector{Float64}}}()
    open("rk.txt","r") do io
        for T in 1:Tstep
            this_step = Vector{Vector{Float64}}()
            for i in 1:(n*m+n*(m-1))
                push!(this_step, [parse(Float64, x) for x in split(readline(io))])
            end
            push!(this_step, [0.0, 0.0, 0.1])
            push!(repoint, this_step)
        end
        for T in 1:Tstep
            this_step = Vector{Vector{Float64}}()
            for i in 1:(n*m+n*(m-1))
                push!(this_step, [parse(Float64, x) for x in split(readline(io))])
                # this_step[end] = renorm(this_step[end])
            end
            if T == 1
                @info this_step
            end
            push!(this_step, [0.0, 0.0, 0.0])
            push!(refield, this_step)
        end
    end
    return repoint, refield
end
# f = Figure(size = (800, 800))
# as = Axis3(f[1, 1], aspect = (6, 6, 3))

n=3
m=4

point_data, field_data = init_data(n, m)


ps = [Point3f(x, y, z) for y in 0:3:3*(m-1) for x in 0:3:3*(n-1) for z in 0:0]
append!(ps, [Point3f(x, y, z) for y in 0:3:3*(m-2) for x in 0:3:3*(n-1) for z in 4:4])
push!(ps,Point3f(3*(n-1)+3, 3*(m-1)+3, 8))
ns_point = Observable([Vec3f(point_data[1][i]...) for i in 1:length(point_data[1])])
ns_field = Observable([Vec3f(field_data[1][i]...) for i in 1:length(field_data[1])])
f = Figure(size = (800, 800))
ax = Axis3(f[1, 1])
arrows!(ax,
    ps, ns_point, fxaa=true, # turn on anti-aliasing
    linecolor = :gray, arrowcolor = :black,
    linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.5),
    align = :center, depth_shift = 1
)
arrows!(ax, 
    ps, ns_field, fxaa=true, # turn on anti-aliasing
    linecolor = :red, arrowcolor = :black,
    linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.5),
    align = :center, depth_shift = 1
)
scatter!(ps[1:n*m], color = :blue, markersize = 10, fxaa=true)
scatter!(ps[n*m+1:n*(2m-1)], color = :red, markersize = 10, fxaa=true)
# scatter!(f[1,2], ps[1:n*m], color = :blue, markersize = 10, fxaa=true)
# scatter!(f[1,2], ps[n*m+1:n*(2m-1)], color = :red, markersize = 10, fxaa=true)
save("my.png", fig)

record(f, "my.mp4", 1:50:100000) do val
    @info "val = $val"
    ns_point[] = [Vec3f(point_data[val][i]...) for i in 1:length(point_data[val])]
    ns_field[] = [Vec3f(field_data[val][i]...) for i in 1:length(field_data[val])]
end