import Optim
using LinearAlgebra

get_constraint(system, fields) = foldl(getfield, fields; init = system)

# This employs an Augmented Lagrangian method
# the magic numbers and arbitrary values are somewhat based on heuristics
# scale is not yet taken into account
function optimize(surfaces, system, v, constraints,
                  aberr = fieldnames(Aberration)[1:5], weights = ones(length(aberr)))
    (; a) = system
    h′ = system.chief.y[end]
    prescription = Prescription(copy(surfaces))
    cts = Dict((c[1],) => c[2] for c in constraints)
    ∑wts = sum(weights)
    # Lagrangian multipliers
    λ = zeros(length(cts))
    # penalty
    μ = Ref(10.0) # since closed over and mutated, this avoids Core.Box heap alloc.
    # init
    x0 = x = prescription[v]
    x0[isinf.(x0)] .= 1e7
    # toleration / patience counter for increasing multipliers
    k = 0
    # monitor constraint satisfaction
    cᵢ =[get_constraint(system, c[1]) - c[2] for c in cts]
    local new_system
    function update!(x)
        prescription[v] .= x
        sys = solve(prescription, a, h′)
        c = [get_constraint(sys, c[1]) - c[2] for c in cts]
        return sys, c
    end
    for i = 1:24
        # Lagrangian / objective function
        function L(x)
            sys, c = update!(x)
            W = aberrations(prescription, sys)
            Wᵢ = abs.(getfield(W, i) for i ∈ aberr)
            return Wᵢ ⋅ weights / ∑wts + 0.5 * μ[] * c ⋅ c + λ ⋅ c
        end
        result = Optim.optimize(L, x0)
        new_system, cᵢ′ = update!(x)
        λ .*= μ[] .* cᵢ′
        # update the penalty parameter if residuals stagnate
        # or multipliers start blowing up
        if norm(cᵢ′) > 0.8 * norm(cᵢ) || k == 3
            μ[] *= 10.0
            k = 0
        elseif norm(λ) > 10.0 * μ[]
            k += 1
        end
        cᵢ = cᵢ′
        x = x0 = Optim.minimizer(result)
    end
    return prescription, new_system
end
