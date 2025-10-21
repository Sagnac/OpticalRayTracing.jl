module MakieExtension

using OpticalRayTracing, Makie

using OpticalRayTracing: System

import OpticalRayTracing: rayplot

import Makie: plot!

function rayplot(x...; theme = Attributes(), kwargs...)
    with_theme(() -> _rayplot(x...; kwargs...), theme)
end

@recipe(RayTracePlot) do scene
    Attributes(;
        ray_colors = (:green, :blue),
        surface_color = :black,
        default_theme(Lines)...
    )
end

_rayplot(system::System; kwargs...) = raytraceplot(system; kwargs...)

function _rayplot(lens::Lens, a::AbstractVector, h′ = -0.5; kwargs...)
    _rayplot(solve(lens, a, h′); kwargs...)
end

function _rayplot(surfaces::Matrix{Float64}, a::AbstractVector, h′ = -0.5; kwargs...)
    _rayplot(Lens(surfaces), a, h′; kwargs...)
end

function _rayplot(
    surfaces::Matrix{Float64},
    system::System,
    a::AbstractVector;
    kwargs...
)
    raytraceplot(surfaces, system, a; kwargs...)
end

function plot!(p::RayTracePlot{Tuple{System}})
    system = p[1][]
    attr = p.attributes
    ray_colors = pop!(attr, :ray_colors)[]
    surface_color = pop!(attr, :surface_color)[]
    z, y... = raypoints(system)
    i = 0
    for yi in y
        i += 1
        ray_color = i > 3 ? ray_colors[1] : ray_colors[2]
        scatter!(p, z, yi; color = ray_color)
        lines!(p, attr, z, yi; color = ray_color)
    end
    zf = z[end]
    yf = y[4][end]
    lines!(p, attr, [zf for i = 1:2], [-yf, yf]; color = surface_color)
    DataInspector()
    return p
end

function plot!(p::RayTracePlot{<:Tuple{Matrix{Float64}, System, <:AbstractVector}})
    surfaces, system, a = (getindex(p, i)[] for i = 1:3)
    attr = p.attributes
    surface_color = pop!(attr, :surface_color)[]
    delete!(attr, :ray_colors)
    z, y... = raypoints(system)
    (; n) = system.lens
    R = @view surfaces[2:end,1]
    z_surfaces = @view z[2:end-1]
    z_corners = Float64[]
    y_corners = Float64[]
    for (i, r) in pairs(R)
        z_vertex = z_surfaces[i]
        ai = a[i]
        z_max = ai ^ 2 / 2r + z_vertex
        if isfinite(r)
            zi = range(z_vertex, z_max, 70)
            y_surface = @.(sqrt(abs(2r * (zi - z_vertex))))
        else
            zi = [z_vertex, z_vertex]
            y_surface = [0.0, ai]
        end
        lines!(p, attr, zi, y_surface; color = surface_color)
        lines!(p, attr, zi, -y_surface; color = surface_color)
        push!(z_corners, z_max)
        push!(y_corners, y_surface[end])
    end
    n_1 = n[1]
    for i = 1:length(z_corners)-1
        n_i = n[i+1]
        n_i == n_1 && continue
        z1, z2 = z_edge = @view z_corners[i:i+1]
        y1, y2 = y_edge = @view y_corners[i:i+1]
        z_min = z_edge[argmin(y_edge)]
        y_max = max(y1, y2)
        z_edge = [z1, z_min, z2]
        y_edge = [y1, y_max, y2]
        lines!(p, attr, z_edge, y_edge; color = surface_color)
        lines!(p, attr, z_edge, -y_edge; color = surface_color)
    end
    raytraceplot!(system)
    return p
end

end
