function vignetting(system::SystemOrRayBasis, a::AbstractVector)
    (; marginal, chief) = system
    ȳ = abs.(@view(chief.y[begin+1:end-1]))
    y = abs.(@view(marginal.y[begin+1:end-1]))
    vig = Matrix{Float64}(undef, length(a), 5)
    vig[:,1] .= a
    vig[:,2] .= y
    unvignetted = vig[:,3] .= y .+ ȳ
    half_vignetted = vig[:,4] .= ȳ
    fully_vignetted = vig[:,5] .= ȳ .- y
    half_vignetted[half_vignetted .< y] .= NaN
    fully_vignetted[fully_vignetted .< y] .= NaN
    a_unvig = a .≥ unvignetted
    if all(a_unvig)
        printstyled("\nUnvignetted.\n\n"; bold = true, color = :green)
    else
        full = findall(a .≤ fully_vignetted)
        partial = setdiff(findall(.!a_unvig), full)
        if !isempty(partial)
            printstyled("\nPartially vignetted:\n"; bold = true, color = :yellow)
            show(partial)
            println('\n')
        end
        if !isempty(full)
            printstyled("\nFully vignetted:\n"; bold = true, color = :red)
            show(full)
            println('\n')
        end
    end
    s1, s2, s3, s4 = (' ' ^ i for i = (7, 4, 7, 5))
    printstyled("a", s1, "limit", s2, "un", s3, "half", s4, "full\n";
                bold = true, color = :cyan)
    return vig
end
