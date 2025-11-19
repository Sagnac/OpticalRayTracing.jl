import Optim
using LinearAlgebra

# This employs an Augmented Lagrangian method
# the magic numbers and arbitrary values are somewhat based on heuristics
function optimize(surfaces, system, v, constraints,
                  aberr = fieldnames(Aberration)[1:5], weights = ones(length(aberr)))
    (; a) = system
    h′ = system.chief.y[end]
    prescription = Prescription(copy(surfaces))
    # normalize keys
    cts = Dict((c[1],) => c[2] for c in constraints)
    # scale constraint violations to improve well-conditioning / stability
    function get_constraints(system)
        Float64[foldl(getfield, c[1]; init = system) / c[2] - 1 for c in cts]
    end
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
    cᵢ = get_constraints(system)
    local new_system
    function update!(x)
        prescription[v] .= x
        sys = solve(prescription, a, h′)
        c = get_constraints(sys)
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
