using JuMP
import HiGHS
using GenericTensorNetworks
using COPT
using Suppressor

function rule110(p, q, r)
    return (q + r + q*r + p*q*r) % 2
end


function generic_logic_grate(p, q, r, N)
    return (N >> (p << 2 | q << 1 | r)) & 1
end

# let state1 - state2 = val
function make_constriant(state1, state2, total_atoms)
    state1_spin = map(x->2*x-1, state1)
    state2_spin = map(x->2*x-1, state2)
    weights = []
    for i in 1:total_atoms
        for j in i+1:total_atoms
            push!(weights, state1_spin[i]*state1_spin[j] - state2_spin[i]*state2_spin[j])
        end
    end
    for i in 1:total_atoms
        push!(weights, state1_spin[i] - state2_spin[i])
    end
    return weights
end

function generate_vaild_mask!(mask_list, now, endpos, now_mask)
    if now == endpos+1
        push!(mask_list, copy(now_mask))
        return 
    end
    for i in 0:255
        now_mask[now] = i
        generate_vaild_mask!(mask_list, now+1, endpos, now_mask)
    end
end

function find_proper_model(ruleid, total_atoms)
    answer = []
    total_mask = Vector{Vector{Int}}()
    generate_vaild_mask!(total_mask, 1, total_atoms - 4, zeros(Int, total_atoms - 4))
    for msk in total_mask 
    # for T in 1:1
        # msk = [1, 15]
        model = @suppress Model(COPT.Optimizer)
        @variable(model, x[1:Int(total_atoms * (total_atoms - 1) / 2 + total_atoms)])

        delta = 0.1
        cnt = 0
        proper_states = Vector{Vector{Int}}()
        for p in [0,1]
            for q in [0,1]
                for r in [0,1]
                    cnt += 1
                    st = [p, q, r, generic_logic_grate(p, q, r, ruleid)]
                    for t in msk
                        push!(st, (t>>(cnt-1))&1)
                    end
                    push!(proper_states, st)
                end
            end
        end

        wrong_states = Vector{Vector{Int}}()
        for thismsk in 0:2^total_atoms-1
            state = [thismsk>>i&1 for i in 0:total_atoms-1]
            if in(state, proper_states) == false
                push!(wrong_states, state)
            end
        end
        @assert length(wrong_states) + length(proper_states) == 2^total_atoms

        for id in 2:length(proper_states)
            weights = make_constriant(proper_states[1], proper_states[id], total_atoms)
            @constraint(model, sum(weights[i] * x[i] for i in 1:length(weights)) == 0)
            # open("test.txt","a") do io
            #     # println(io, weights)
            #     for i in weights
            #         print(io, i, " ")
            #     end
            #     println(io)
            # end
        end

        for id in 1:length(wrong_states)
            weights = make_constriant(proper_states[1], wrong_states[id], total_atoms)
            @constraint(model, sum(weights[i] * x[i] for i in 1:length(weights)) <= -delta)
            # open("test.txt","a") do io
            #     for i in weights
            #         print(io, i, " ")
            #     end
            #     println(io)
            # end
        end

        @objective(model, Min, x[1])

        # @info "ruleid = $ruleid, msk = $msk, proper_states = $(length(proper_states)), wrong_states = $(length(wrong_states))"
        # @info "model = $model"
        @suppress optimize!(model)
        if is_solved_and_feasible(model)
            # @info "msk = $msk, successfully find a solution"
            # weights = value.(x)
            # push!(answer, [maximum(weights), minimum(weights), maximum(weights) / minimum(weights), msk])
            # println(maximum(weights)," ", minimum(weights)," ", maximum(weights) / minimum(weights))
            return msk, [value(x[i]) for i in 1:length(x)]
        else
            @info "msk = $msk, ruleid = $ruleid, failed to find a solution"
        end
    end
    return -1, false
end

function query_model(ruleid, total_atoms)
    msk, weights = find_proper_model(ruleid, total_atoms)
    if msk == -1
        return -1, false
    end
    weights = round.(weights, digits=2)
    weights = weights ./ 0.25
    return msk, weights
    # solution_summary(model)

    # hyperedges = [[i,j] for i in 1:5 for j in i+1:5]
    # for i in 1:5
    #     push!(hyperedges, [i])
    # end
    # hyperspinglass = SpinGlass(5, hyperedges, weights)
    # hyperproblem = GenericTensorNetwork(hyperspinglass)

    # counting_min_eneregy = GenericTensorNetworks.solve(hyperproblem, CountingMin())[]
    # msk = 17
    # cnt = 0
    # output1 = []
    # output2 = []
    # for p in [0,1]
    #     for q in [0,1] 
    #         for r in [0,1]
    #             cnt += 1
    #             st = StaticBitVector([p, q, r, generic_logic_grate(p, q, r, ruleid), (msk>>(cnt-1))&1])
    #             @info st
    #             st = map(x->x⊻1, st)
    #             st = Int.(st)
    #             push!(output1, "$p$q$r$(generic_logic_grate(p, q, r, ruleid))")
    #             push!(output2, "$(st[1])$(st[2])$(st[3])$(st[4])$(st[5])")
    #         end
    #     end
    # end
end

natoms = Vector{Vector{Int}}()
previous = Vector{Int}() 
for total_atoms = 6:1:6
    ok = Vector{Int}()
    for id in 0:255
        if (ok in previous) == false
            msk, weights = query_model(id, total_atoms)
            if msk != -1
                push!(ok, id)
                push!(previous, id)
            end
        end
    end
    push!(natoms, ok)
end



# output_weights = zeros(Int, 5, 5)

# for i in 1:10
#     output_weights[hyperedges[i][1], hyperedges[i][2]] = weights[i]
# end
# for i in 11:15
#     output_weights[hyperedges[i][1], hyperedges[i][1]] = weights[i]
# end
# for i in 1:5
#     for j in 1:i-1
#         output_weights[i, j] = output_weights[j, i]
#     end
# end

# for i in 1:5
#     for j in 1:5
#         print(output_weights[i, j], " ")
#     end
#     println("")
# end

function individualsolver()
    model = @suppress Model(COPT.Optimizer)
    # set_attribute(model, "output_flag", false)
    @variable(model, x[1:15+6])
    open("test.txt","r") do io
        for t in 1:7
            weights = parse.(Int, split(readline(io)))
            @constraint(model, sum(weights[i] * x[i] for i in 1:length(weights)) == 0)
        end
        for t in 8:63
            weights = parse.(Int, split(readline(io)))
            @constraint(model, sum(weights[i] * x[i] for i in 1:length(weights)) <= -0.1)
        end
    end
    @objective(model, Min, x[1])
    # @info "$model"
    @suppress optimize!(model)
    if is_solved_and_feasible(model)
        @info "yes"
    else
        @info "no"
    end
end
individualsolver()