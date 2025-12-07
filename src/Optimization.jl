# This employs an Augmented Lagrangian method
# the magic numbers and arbitrary values are somewhat based on heuristics
function optimize(surfaces, system, v, constraints,
                  aberr = fieldnames(Aberration)[1:5], weights = ones(length(aberr)))
    (; a) = system
    h′ = system.chief.y[end]
    layout = Layout(copy(surfaces))
    # normalize keys
    cts = Dict((c[1],) => c[2] for c in constraints)
    # scale constraint violations to improve well-conditioning / stability
    function get_constraints(system)
        cx = Vector{Float64}(undef, length(cts))
        for (i, c) in enumerate(cts)
            c1 = foldl(getfield, c[1]; init = system)
            c2 = c[2]
            cx[i] = iszero(c2) ? c1 : c1 / c2 - 1
        end
        return cx
    end
    ∑wts = sum(weights)
    # Lagrangian multipliers
    λ = zeros(length(cts))
    # penalty
    μ = Ref{Float64}(10.0) # since closed over and mutated, this avoids Core.Box
    # init
    x0 = layout[v]
    x0[isinf.(x0)] .= 1e7
    # toleration / patience counter for increasing multipliers
    k = 0
    # monitor constraint satisfaction
    c_prev = get_constraints(system)
    function update!(x)
        layout[v] .= x
        local system = solve(layout, a, h′)
        local c = get_constraints(system)
        return system, c
    end
    for i = 1:24
        # Lagrangian / objective function
        function L!(x)
            local system, c = update!(x)
            W = aberrations(layout, system)
            Wᵢ = abs.(getfield(W, i) for i ∈ aberr)
            return Wᵢ ⋅ weights / ∑wts + 0.5 * μ[] * c ⋅ c + λ ⋅ c
        end
        result = Optim.optimize(L!, x0)
        x0 = Optim.minimizer(result)
        system, c = update!(x0)
        λ .*= μ[] .* c
        # update the penalty parameter if residuals stagnate
        # or multipliers start blowing up
        if norm(c) > 0.8 * norm(c_prev) || k == 3
            μ[] *= 10.0
            k = 0
        elseif norm(λ) > 10.0 * μ[]
            k += 1
        end
        c_prev = c
    end
    return system
end

function optimize(
    system::System{Layout}, v::Vector{Int}, constraints::Dict,
    aberr = fieldnames(Aberration)[1:5],
    weights = ones(length(aberr))
)
    return optimize(system.layout, system, v, constraints, aberr, weights)
end
