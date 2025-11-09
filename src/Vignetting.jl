function _vignetting(system::SystemOrRayBasis, a::AbstractVector, stop::Int)
    (; marginal, chief) = system
    ȳ = abs.(surface_ray(chief.y))
    y = abs.(surface_ray(marginal.y))
    vig = Matrix{Float64}(undef, length(a), 5)
    vig[:,1] .= a
    limited = vig[:,2] .= y
    unvignetted = vig[:,3] .= y .+ ȳ
    half_vignetted = vig[:,4] .= ȳ
    fully_vignetted = vig[:,5] .= ȳ .- y
    half_vignetted[half_vignetted .< y] .= NaN
    fully_vignetted[fully_vignetted .< y] .= NaN
    a_unvig = a .≥ unvignetted .|| .≈(a, unvignetted)
    limit = partial = full = Int[]
    un = all(a_unvig)
    min_un = minimum((a[i] - y[i]) / ȳ[i] for i in eachindex(a) if i != stop)
    min_half = minimum(a ./ ȳ)
    min_full = minimum((a .+ y) ./ ȳ)
    FOV = Matrix{Float64}(undef, 3, 3)
    for (i, s) in enumerate((min_un, min_half, min_full))
        ū = abs(chief.u[1] * s)
        FOV[i,1] = 2atand(ū)
        FOV[i,2] = ū
        FOV[i,3] = abs(chief.y[end] * s)
    end
    limit = findall(a .< limited .&& .!(a .≈ unvignetted))
    full = findall(a .≤ fully_vignetted)
    partial = setdiff(findall(.!a_unvig), full)
    return Vignetting(vig, FOV, un, limit, partial, full)
end

function vignetting(system::System, a::AbstractVector = system.a)
    _vignetting(system, system.a, system.stop)
end

function vignetting(rays::RayBasis, system::System, a::AbstractVector = system.a)
    _vignetting(rays, a, system.stop)
end
