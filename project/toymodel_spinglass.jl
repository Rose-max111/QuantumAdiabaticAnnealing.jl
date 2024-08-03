using QuantumAdiabaticAnnealing

function evaluate_adiabatic_time(n, m)
    Vtrans = fill(1.0, n*m + n*(m-1))
    adiabatic_time = 0.0
    max_try = 1.0
    while max_try > 0
        next_try = max_try
        sp = spinglass_mapping(n, m)
        freeze_input!(sp)
        @info "Stage1, max_try = $max_try, next_try = $next_try, checking if vaild"
        flag = check_vaild_time(sp, next_try, Vtrans)
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
        sp = spinglass_mapping(n, m)
        freeze_input!(sp)
        @info "Stage2, max_try = $max_try, next_try = $next_try, checking if vaild"
        flag = check_vaild_time(sp, next_try, Vtrans)
        if flag == false
            adiabatic_time = next_try
        end
    end
    return adiabatic_time + 1
end

