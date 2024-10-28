using QuantumAdiabaticAnnealing
using Random

function ioinfile()
    edge = Vector{Tuple{Int, Int}}()
    weights = Vector{Float64}()
    
    open("simple.txt","r") do file
        for line in eachline(file)
            thisline = split(line)
            println(thisline)
            push!(edge, (parse(Int, thisline[1]), parse(Int, thisline[2])))
            push!(weights, parse(Float64, thisline[3]))
        end
    end
    return edge, weights
end

function initial_randomize!(sp)
    for i in 1:length(sp.onsite)
        aa = rand()*0.2 - 0.1
        bb = rand()*0.2 - 0.1
        # @info "aa = $aa, bb=$bb, i=$i"
        sp.M[i][3] = bb
        sp.M[i][2] = aa
        sp.M[i][1] = -sqrt(1 - bb^2 - aa^2)
        # @info "i=$i, sp.M[i]= $(sp.M[i])"
    end
end

function randominit(n)
    edge = Vector{Tuple{Int, Int}}()
    weights = Vector{Float64}()
    total = n*(n-1) / 2
    edge_set = [(i,j) for i in 1:n for j in i+1:n]

    for howmanyedge in 1:Int(floor(total * 0.15))
        while true
            ed = rand(edge_set)
            if !(ed in edge)
                # @info "ed = $ed, edge = $edge, flag = $(ed in edge)"
                push!(edge, (ed[1], ed[2]))
                push!(weights, rand() > 0.5 ? 1 : - 1)
                break
            end
        end
    end
    edge ,weights 
end

# edge, weights = ioinfile()
# sp = spinglass_random_mapping(108, edge, weights)
edges, weights = randominit(40)
sp = spinglass_random_mapping(40, edges, weights)
initial_randomize!(sp)

# @info "energy is $(sp_energy(sp, 0, 1, fill(1.0, length(sp.onsite))))"
sp_ground_state(sp)

runge_kutta_integrate!(sp, 1e-2, 1000.0, fill(1.0, length(sp.onsite));pin_input = false)
for i in 1:length(sp.onsite)
    sp.M[i][3] = sp.M[i][3] > 0 ? 1 : -1
end
@info "energy is $(sp_energy(sp, 1, 1, fill(1.0, length(sp.onsite))))"