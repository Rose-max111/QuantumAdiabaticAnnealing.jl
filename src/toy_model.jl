abstract type TransitionRule end
struct HeatBath <: TransitionRule end
struct Metropolis <: TransitionRule end
function update(::Metropolis, Temp, ΔE, prior)
    if ΔE < 0
        1.0
    else
        exp(- (ΔE) / Temp)
    end * prior
end

function update(::HeatBath, Temp, ΔE, prior)
    exp(-ΔE / Temp) / (1 + exp(-ΔE / Temp)) * prior
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

function query_condition_with_gradient(state, l, mid, r, output, m, row; energy_gradient = nothing)
    ret = query_condition(state, l, mid, r, output)
    if energy_gradient != nothing
        return ret * energy_gradient^(m - row - 1)
    end
    return ret 
end

function toy_model_state_energy(state, n, m; on_site_energy = nothing, period_condition = false, energy_gradient = nothing)
    ret = 0
    if on_site_energy == nothing
        on_site_energy = fill(0, n)
    end
    for l in 1:n
        ret += on_site_energy[l] * ((state & (2^(l-1))) > 0 ? 1 : -1)
    end
    if period_condition == false
        for l in 1:n-2
            for row in 1:m-1
                actual_l = (row-1) * n + l
                ret += query_condition_with_gradient(state, actual_l, actual_l+1, actual_l+2, actual_l+n, m, row; energy_gradient)
            end
        end
    else
        for row in 1:m-1
            for mid in 2:n-1
                real_mid = (row-1) * n + mid
                ret += query_condition_with_gradient(state, real_mid-1, real_mid, real_mid+1, real_mid+n, m, row; energy_gradient)
            end
            real_mid = (row-1) * n + 1
            ret += query_condition_with_gradient(state, (row-1) * n + n, real_mid, real_mid+1, real_mid+n, m, row; energy_gradient)
            real_mid = (row-1) * n + n
            ret += query_condition_with_gradient(state, real_mid-1, real_mid, (row-1) * n + 1, real_mid+n, m, row; energy_gradient)
        end
    end
    return ret
end


function toy_model_transition_matrix(rule::TransitionRule, n, m, Temp; on_site_energy = nothing, period_condition = false, energy_gradient = nothing)
    total_atoms = n * m - 2
    if period_condition == true
        total_atoms = n * m
    end
    if on_site_energy == nothing
        on_site_energy = fill(0, n)
    end
    state_energy = [toy_model_state_energy(i, n, m; on_site_energy, period_condition, energy_gradient) for i in 0:(2^total_atoms - 1)]

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

            ΔE = state_energy[reverse_mask + 1] - state_energy[i + 1]
            push!(val, update(rule, Temp, ΔE, 1.0 / total_atoms))
            other_prob += val[end]
        end
        push!(row, i+1)
        push!(col, i+1)
        push!(val ,1 - other_prob)
    end
    # @info val
    return sparse(row, col, val)
end