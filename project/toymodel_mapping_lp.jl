using JuMP
using GenericTensorNetworks
using COPT
using HiGHS
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
    # msk = [1, 15]
    # for T in 1:1
        model = @suppress Model(HiGHS.Optimizer)
        set_silent(model)
        @variable(model, x[1:Int(total_atoms * (total_atoms - 1) / 2 + total_atoms)])
        delta = 1
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
        end

        for id in 1:length(wrong_states)
            weights = make_constriant(proper_states[1], wrong_states[id], total_atoms)
            @constraint(model, sum(weights[i] * x[i] for i in 1:length(weights)) <= -delta)
        end

        # @objective(model, Min, x[1])
        optimize!(model)
        if is_solved_and_feasible(model)
            return msk, [value(x[i]) for i in 1:length(x)]
        else
            @info "msk = $msk, ruleid = $ruleid, failed to find a solution"
        end
        for con in all_constraints(model; include_variable_in_set_constraints = false)
            delete(model, con)
        end
    end
    return -1, false
end

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

function check_vaild(total_atoms, weights, ruleid, msk)
    hyperedges = [[i, j] for i in 1:total_atoms for j in i+1:total_atoms]
    for i in 1:total_atoms
        push!(hyperedges, [i])
    end
    hyperweights = weights
    spgls = SpinGlass(total_atoms, hyperedges, hyperweights)
    spproblem = GenericTensorNetwork(spgls)
    gs = GenericTensorNetworks.solve(spproblem, CountingMin())[]
    @assert gs.c == 8.0
    cnt = 0
    for p in [0, 1]
        for q in [0, 1]
            for r in [0, 1]
                cnt += 1
                state = [p, q, r, generic_logic_grate(p, q, r, ruleid)]
                for i in 1:total_atoms-4
                    push!(state, (msk[i]>>(cnt-1))&1)
                end
                state .⊻= 1
                stbit = StaticBitVector(state)
                this_energy = spinglass_energy(spgls, stbit)
                # @info "state = $(state.⊻1), energy = $this_energy"
                @assert this_energy == gs.n
            end
        end
    end
end

function query_model(ruleid, total_atoms)
    msk, weights = find_proper_model(ruleid, total_atoms)
    if msk == -1
        return -1, false
    end
    weights = round.(weights, digits=2)
    weights = weights ./ 0.25
    return msk, weights
end

function write_data(total_atoms, weights, ruleid, msk)
    filename = (pwd()) * "/data/spin_glass_mapping/$(ruleid).txt"
    weights = map(x -> x==-0.0 ? 0.0 : x, weights)
    open(filename, "w") do io
        println(io, "Total Atoms = $total_atoms")
        println(io, "")
        cnt = 0
        for i in 1:total_atoms
            for j in i+1:total_atoms
                cnt = cnt + 1
                println(io,"J($i, $j)=", weights[cnt])
            end
        end
        for i in 1:total_atoms
            println(io, "h($i)=$(weights[cnt+i])")
        end

        println(io, "")
        hyperedges = [[i, j] for i in 1:total_atoms for j in i+1:total_atoms]
        for i in 1:total_atoms
            push!(hyperedges, [i])
        end
        hyperweights = weights
        spgls = SpinGlass(total_atoms, hyperedges, hyperweights)
        spproblem = GenericTensorNetwork(spgls)
        gs = GenericTensorNetworks.solve(spproblem, CountingMin())[]
        println(io, "Ground state energy = $(gs.n), ", "Ground state degenerate = $(gs.c)")

        println(io, "")
        println(io, "The 8 degenerate ground states are:")
        cnt = 0
        for p in [0,1]
            for q in [0,1]
                for r in [0,1]
                    state = [p, q, r, generic_logic_grate(p, q, r, ruleid)]
                    for i in 1:total_atoms-4
                        push!(state, (msk[i]>>(cnt-1))&1)
                    end
                    state = map(x -> 2*x-1, state)
                    println(io, state)
                end
            end
        end
    end
end

function __main__()
    natoms = Vector{Vector{Int}}()
    previous = Vector{Int}()
    for total_atoms = 4:1:5
        ok = Vector{Int}()
        for id in 0:255
            if (id in previous) == false
                msk, weights = query_model(id, total_atoms)
                if msk != -1
                    push!(ok, id)
                    push!(previous, id)
                    @info "now testing ruleid = $id, total_atoms = $total_atoms, weights = $(weights)"
                    check_vaild(total_atoms, weights, id, msk)
                    write_data(total_atoms, weights, id, msk)
                end
            end
        end
        push!(natoms, ok)
    end
end

__main__()




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