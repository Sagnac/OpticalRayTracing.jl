import Optim

get_constraint(system, fields) = foldl(getfield, fields; init = system)

function optimize(surfaces, system, v, constraints,
                  aberr = fieldnames(Aberration)[1:5], weights = ones(length(aberr)))
    (; a) = system
    h′ = system.chief.y[end]
    prescription = Prescription(copy(surfaces))
    cts = Dict((c[1],) => c[2] for c in constraints)
    ∑weights = sum(weights)
    p = 1e3
    function objective(u)
        prescription[v] .= u
        sys = solve(prescription, a, h′)
        c = [get_constraint(sys, c[1]) - c[2] for c in cts]
        W = aberrations(prescription, sys)
        return abs.(getfield(W, i) for i in aberr)' * weights / ∑weights + p * c' * c
    end
    u0 = prescription[v]
    result = Optim.optimize(objective, u0)
    u = Optim.minimizer(result)
    prescription[v] .= u
    system = solve(prescription, a, h′)
    return prescription, system
end
