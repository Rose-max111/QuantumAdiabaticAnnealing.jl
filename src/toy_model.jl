function rule110(p, q, r)
    return (q + r + q*r + p*q*r) % 2
end

function get_bit(state, l)
    return (state & (2^(l-1))) > 0 ? 1 : 0
end

function query_condition(state, l, mid, r, output)
    p = get_bit(state, l)
    q = get_bit(state, mid)
    r = get_bit(state, r)
    next = get_bit(state, output)
    return rule110(p, q, r) ⊻ next
end

function toy_model_state_energy(state, n; on_site_energy = nothing, period_condition = false)
    ret = 0
    if on_site_energy == nothing
        on_site_energy = fill(0, n)
    end
    for l in 1:n
        ret += on_site_energy[l] * ((state & (2^(l-1))) > 0 ? 1 : -1)
    end
    if period_condition == false
        for l in 1:n-2
            ret += query_condition(state, l, l+1, l+2, l+n)
        end
    else
        for mid in 2:n-1
            ret += query_condition(state, mid-1, mid, mid+1, mid+n)
        end
        ret += query_condition(state, n, 1, 2, 1+n)
        ret += query_condition(state, n-1, n, 1, n+n)
    end
    return ret
end

function toy_model_transition_matrix(n, Temp; on_site_energy = nothing, period_condition = false)
    total_atoms = 2 * n - 2
    if period_condition == true
        total_atoms = 2 * n
    end
    if on_site_energy == nothing
        on_site_energy = fill(0, n)
    end
    state_energy = [toy_model_state_energy(i, n; on_site_energy, period_condition) for i in 0:(2^total_atoms - 1)]

    row = Vector{Int}()
    col = Vector{Int}()
    val = Vector{Float64}()
    @info "finish calculating state energy"
    # P = spzeros(2^total_atoms, 2^total_atoms)
    for i in 0:(2^total_atoms - 1)
        other_prob = 0
        for j in 1:total_atoms
            reverse_mask = i ⊻ (2^(j-1))
            push!(row, reverse_mask + 1)
            push!(col, i+1)
            if state_energy[reverse_mask + 1] < state_energy[i + 1] # accept transition
                push!(val, 1.0 / total_atoms)
                # P[reverse_mask + 1, i + 1] = 1.0 / total_atoms
            else
                push!(val, exp(- (state_energy[reverse_mask + 1] - state_energy[i + 1]) / Temp) * 1.0 / total_atoms)
                # P[reverse_mask + 1, i + 1] = exp(- (state_energy[reverse_mask + 1] - state_energy[i + 1]) / Temp)
            end
            other_prob += val[end]
        end
        push!(row, i+1)
        push!(col, i+1)
        push!(val ,1 - other_prob)
    end
    # @info val
    return sparse(row, col, val)
end