function _vignetting(system::SystemOrRayBasis, a::AbstractVector)
    (; marginal, chief, stop) = system
    ȳ = abs.(@view(chief.y[begin+1:end-1]))
    y = abs.(@view(marginal.y[begin+1:end-1]))
    vig = Matrix{Float64}(undef, length(a), 5)
    vig[:,1] .= a
    limited = vig[:,2] .= y
    unvignetted = vig[:,3] .= y .+ ȳ
    half_vignetted = vig[:,4] .= ȳ
    fully_vignetted = vig[:,5] .= ȳ .- y
    half_vignetted[half_vignetted .< y] .= NaN
    fully_vignetted[fully_vignetted .< y] .= NaN
    a_unvig = a .≥ unvignetted .|| .≈(a, unvignetted)
    bold = true
    limit = partial = full = Int[]
    if all(a_unvig)
        printstyled("Unvignetted.\n\n"; bold, color = :green)
    else
        limit = findall(a .< limited .&& .!(a .≈ unvignetted))
        full = findall(a .≤ fully_vignetted)
        partial = setdiff(findall(.!a_unvig), full)
        if !isempty(partial)
            printstyled("Partially vignetted:\n"; bold, color = :yellow)
            show(partial)
            println('\n')
        end
        if !isempty(full)
            printstyled("Fully vignetted:\n"; bold, color = :magenta)
            show(full)
            println('\n')
        end
        if !isempty(limit)
            printstyled("Stop limited:\n"; bold, color = :red)
            show(limit)
            println('\n')
        end
    end
    min_un = minimum((a[i] - y[i]) / ȳ[i] for i in eachindex(a) if i != stop)
    min_half = minimum(a ./ ȳ)
    min_full = minimum((a .+ y) ./ ȳ)
    FOVs = [(2atand(abs(chief.u[1] * s)), abs(chief.y[end] * s)) for s in
        (min_un, min_half, min_full)]
    println("Maximum supported FOVs:")
    @printf("Unvignetted: %.4f°, h′: %.4f\n", FOVs[1]...)
    @printf("Half-vignetted: %.4f°, h′: %.4f\n", FOVs[2]...)
    @printf("Fully-vignetted: %.4f°, h′: %.4f\n\n", FOVs[3]...)
    s1, s2, s3, s4 = (' ' ^ i for i = (7, 4, 7, 5))
    printstyled("a", s1, "limit", s2, "un", s3, "half", s4, "full\n";
                bold, color = :cyan)
    return vig, FOVs, limit, partial, full
end

vignetting(system::SystemOrRayBasis, a::AbstractVector) = _vignetting(system, a)[1]

function vignetting(system::SystemOrRayBasis, a::AbstractVector, ::Type{FOV})
    r = redirect_stdout(devnull) do
        _vignetting(system, a)[2]
    end
    ū, h′ = eachrow(reinterpret(reshape, Float64, r))
    (tand.(0.5 * ū), copy(h′))
end
