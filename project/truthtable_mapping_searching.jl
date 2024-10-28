using GenericTensorNetworks
using Random

function rule110(p, q, r)
    return (q + r + q*r + p*q*r) % 2
end

function search_weighted(id, length, this_set, set, edge_set)
    if id == length+1
        flag = search_spinglass_symmetry(this_set, edge_set)
        if flag == false
            return false
        end
        return true
    else
        for w in set
            this_set[id] = w
            if id == 2
                @info "now begin 1, 2 is $(this_set[1]), $(this_set[2])"
            end
            if search_weighted(id+1, length, this_set, set, edge_set) == true
                return true
            end
        end
    end
    return false
end

@inline function query_energy(msk, hyperedges, weights_selected)
    sum = 0
    for iid in 1:length(hyperedges)
        x = 1
        for u in hyperedges[iid]
            x *= ((msk >> (u-1)) & 1 == 1) ? -1 : 1
        end
        sum += x*weights_selected[iid]
    end
    return sum
end

function search_spinglass_symmetry(weights_set, edge_set)
    hyperedges = copy(edge_set)
    weights_selected = copy(weights_set)
    push!(hyperedges, [1, 3])
    push!(hyperedges, [3, 4])
    push!(hyperedges, [3])
    append!(weights_selected, [weights_selected[end-2], weights_selected[end-1], weights_selected[end]])

    val0 = query_energy(0, hyperedges, weights_selected)
    minn = val0
    for msk in 1:2^3-1
        sum = query_energy((msk)+(rule110(msk&1, (msk>>1)&1, (msk>>2)&1)<<3), hyperedges, weights_selected)
        if sum != minn
            return false
        end
    end
    @info "weights_set = $weights_set, minn_lst = $(minn), cnt = $cnt"
    for msk in 0:2^3-1
        sum = query_energy((msk)+((rule110(msk&1, (msk>>1)&1, (msk>>2)&1)⊻1)<<3), hyperedges, weights_selected)
        @info "msk = $((msk)+((rule110(msk&1, (msk>>1)&1, (msk>>2)&1)⊻1)<<3)), energy = $sum"
        if sum <= minn
            return false
        end
    end
    return true
    # hyperspinglass = SpinGlass(4, hyperedges, weights_selected)
    # hyperproblem = GenericTensorNetwork(hyperspinglass)
    # counting_min_eneregy = solve(hyperproblem, CountingMin())[]
    # @assert counting_min_eneregy.n == minn
    # @assert counting_min_eneregy.c == cnt
    # @info "counting_min_eneregy = $(counting_min_eneregy)"
    # return false
    @info "weights_set = $weights_set, minn = $(minn), cnt = $cnt"
    for input in inputs
        this_energy = 0
        for iid in 1:length(hyperedges)
            x = 1
            for u in hyperedges[iid]
                x *= input[u] == 1 ? -1 : 1
            end
            this_energy += x*weights_selected[iid]
        end
        @info "input = $input, this_energy = $this_energy"
        if this_energy != minn
            return false
        end
    end
    return true
end

num_vertices = 4
hyperedges = [[1,4], [2,3]]
for i in [1,4]
    push!(hyperedges, [i])
end
append!(hyperedges, [[1,2], [2,4], [2]])

search_weighted(1, length(hyperedges), zeros(length(hyperedges)), Vector(13:-1:-13), hyperedges)

# weights_set = []
# search_weighted!(weights_set, 1, length(hyperedges), zeros(length(hyperedges)), [-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5])

# for id in 1:length(weights_set)
#     if search_spinglass_symmetry(id, weights_set, hyperedges) == true
#         println("weights_set = $(weights_set[id])")
#         break
#     end
# end