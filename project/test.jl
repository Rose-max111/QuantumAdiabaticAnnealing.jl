using DormandPrince

function test(y)
    return [0.85*y[1], 0.75*y[2]]
end

function fcn(x, y, f)
    f = test(y)
    # f[1] = 0.85*y[1]
    # f[2] = 0.75*y[2]
end

solver = DP5Solver(fcn, 0.0, [10.0, 20.0])
integrate!(solver, 5.0)

get_current_state(solver)
