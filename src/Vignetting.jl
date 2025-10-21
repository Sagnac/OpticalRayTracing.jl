function vignetting(system::System, a::AbstractVector)
    (; marginal, chief) = system
    k = length(a)
    ȳ = abs.(@view(chief.y[begin+1:end-1]))
    y = abs.(@view(marginal.y[begin+1:end-1]))
    vig = Matrix{Float64}(undef, k, 4)
    vig[:,1] .= a
    unvignetted = vig[:,2] .= y .+ ȳ
    half_vignetted = vig[:,3] .= max.(ȳ, y)
    fully_vignetted = vig[:,4] .= max.(ȳ .- y, y)
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
            printstyled("\nFully vignetted or limit:\n"; bold = true, color = :red)
            show(full)
            println('\n')
        end
    end
    s = ' ' ^ 8
    printstyled("a", s, "un", s, "half", s, "full\n"; bold = true, color = :cyan)
    return vig
end
