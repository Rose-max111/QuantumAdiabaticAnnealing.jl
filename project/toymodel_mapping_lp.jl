using JuMP
import HiGHS
using GenericTensorNetworks

function rule110(p, q, r)
    return (q + r + q*r + p*q*r) % 2
end


# let state1 - state2 = val
function make_constriant(model, state1, state2, val = nothing)
    state1_spin = map(x->2*x-1, state1)
    state2_spin = map(x->2*x-1, state2)
    weights = []
    for i in 1:5
        for j in i+1:5
            push!(weights, state1_spin[i]*state1_spin[j] - state2_spin[i]*state2_spin[j])
        end
    end
    for i in 1:5
        push!(weights, state1_spin[i] - state2_spin[i])
    end
    return weights
end

function find_proper_model()
    answer = []
    for msk in 0:2^8-1
    # for msk in 238:238
        model = Model(HiGHS.Optimizer)
        @variable(model, x[1:15])

        delta = 1
        cnt = 0
        proper_states = Vector{Vector{Int}}()
        for p in [0,1]
            for q in [0,1]
                for r in [0,1]
                    cnt += 1
                    push!(proper_states, [p, q, r, rule110(p, q, r), (msk>>(cnt-1))&1])
                end
            end
        end

        # for i in proper_states
        #     println(i)
        # end

        wrong_states = Vector{Vector{Int}}()
        for msk in 0:2^5-1
            state = [msk>>i&1 for i in 0:4]
            if in(state, proper_states) == false
                push!(wrong_states, state)
            end
        end

        @assert length(wrong_states) + length(proper_states) == 2^5

        for id in 2:length(proper_states)
            weights = make_constriant(model, proper_states[1], proper_states[id], nothing)
            @constraint(model, sum(weights[i] * x[i] for i in 1:length(weights)) == 0)
        end

        for id in 1:length(wrong_states)
            weights = make_constriant(model, proper_states[1], wrong_states[id], delta)
            @constraint(model, sum(weights[i] * x[i] for i in 1:length(weights)) <= -delta)
        end

        @objective(model, Min, x[1])

        optimize!(model)
        if is_solved_and_feasible(model)
            @info "msk = $msk, successfully find a solution"
            # weights = value.(x)
            # push!(answer, [maximum(weights), minimum(weights), maximum(weights) / minimum(weights), msk])
            # println(maximum(weights)," ", minimum(weights)," ", maximum(weights) / minimum(weights))
            return [value(x[i]) for i in 1:15]
        else
            @info "msk = $msk, failed to find a solution"
        end
    end
end

weights = find_proper_model()
weights = round.(weights, digits=2)
weights = weights ./ 0.25
# solution_summary(model)

hyperedges = [[i,j] for i in 1:5 for j in i+1:5]
for i in 1:5
    push!(hyperedges, [i])
end
hyperspinglass = SpinGlass(5, hyperedges, weights)
hyperproblem = GenericTensorNetwork(hyperspinglass)

counting_min_eneregy = GenericTensorNetworks.solve(hyperproblem, CountingMin())[]
msk = 17
cnt = 0
output1 = []
output2 = []
for p in [0,1]
    for q in [0,1] 
        for r in [0,1]
            cnt += 1
            st = StaticBitVector([p, q, r, rule110(p, q, r), (msk>>(cnt-1))&1])
            @info st
            st = map(x->x⊻1, st)
            st = Int.(st)
            push!(output1, "$p$q$r$(rule110(p, q, r))")
            push!(output2, "$(st[1])$(st[2])$(st[3])$(st[4])$(st[5])")
        end
    end
end
for i in output1
    print(i, " | ")
end
# println("")
for i in output2
    print(i, " | ")
end


output_weights = zeros(Int, 5, 5)

for i in 1:10
    output_weights[hyperedges[i][1], hyperedges[i][2]] = weights[i]
end
for i in 11:15
    output_weights[hyperedges[i][1], hyperedges[i][1]] = weights[i]
end
for i in 1:5
    for j in 1:i-1
        output_weights[i, j] = output_weights[j, i]
    end
end

for i in 1:5
    for j in 1:5
        print(output_weights[i, j], " ")
    end
    println("")
end