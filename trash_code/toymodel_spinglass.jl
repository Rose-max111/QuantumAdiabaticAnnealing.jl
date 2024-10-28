using QuantumAdiabaticAnnealing

function evaluate_adiabatic_time(n, m, gradient)
    # Vtrans = fill(1.0, n*m + n*(m-1))
    adiabatic_time = 0.0
    max_try = 1.0
    while max_try > 0
        next_try = max_try
        @info "Stage1, max_try = $max_try, next_try = $next_try, checking if vaild"
        sp = spinglass_adiabatic_dp8(n, m, next_try, gradient)
        flag = sp_check_vaild(sp)
        if flag == true
            break
        end
        max_try *= 2
    end
    if max_try == 1
        return 1
    end
    adiabatic_time += max_try / 2
    max_try /= 2
    while max_try != 1
        max_try /= 2
        next_try = adiabatic_time + max_try
        @info "Stage2, max_try = $max_try, next_try = $next_try, checking if vaild"
        sp = spinglass_adiabatic_dp8(n, m, next_try, gradient)
        flag = sp_check_vaild(sp)
        if flag == false
            adiabatic_time = next_try
        end
    end
    return adiabatic_time + 1
end

evaluate_adiabatic_time(3, 4, 1)

sp_ground_state_sa(3, 5)
sp_ground_state(spinglass_mapping(3, 5))
