function raypoints(system::System)
    (; lens) = system
    (; M, n) = lens
    (; marginal, chief, EBFD) = system
    t = map(*, @view(M[:,1]), @view(n[begin:end-1]))
    k = length(t) + 2
    l = sum(t)
    BFD = isfinite(EBFD) ? abs(EBFD) * n[end] : l / 4
    l += BFD
    s = 0.1
    d = -l * s / n[1]
    at_objective = iszero(t[1])
    z = Vector{Float64}(undef, k)
    z[1] = at_objective ? d : 0.0
    z[2] = at_objective ? 0.0 : t[1]
    for i = 2:lastindex(t)
        z[i+1] = z[i] + t[i]
    end
    z[end] = z[end-1] + BFD
    y0 = zeros(k)
    y1 = marginal.y
    y2 = -y1
    nū = chief.nu[1]
    ȳ1 = chief.y[2]
    ȳo1 = ȳ1 + nū * d
    ȳ = [ȳo1; @view(chief.y[begin+1:end])]
    y3 = ȳ + y1
    y4 = ȳ + y2
    return z, y0, y1, y2, ȳ, y3, y4
end

# requires Makie
function rayplot end
function rayplot! end
