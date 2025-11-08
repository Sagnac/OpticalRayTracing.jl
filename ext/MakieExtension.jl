module MakieExtension

using OpticalRayTracing, Makie, Printf

using OpticalRayTracing: System, RayBasis, SystemOrRayBasis, Aberration

import OpticalRayTracing: rayplot, rayplot!, wavefan, rayfan,
                          field_curves, percent_distortion, spot_size, caustic

import Makie: plot!

const default_plot_theme = Attributes()

function rayplot(x...; theme = default_plot_theme, kwargs...)
    with_theme(() -> _rayplot(x...; kwargs...), theme)
end

function rayplot!(x...; theme = default_plot_theme, kwargs...)
    with_theme(() -> _rayplot!(x...; kwargs...), theme)
end

@recipe RayTracePlot begin
    ray_colors = (:blue, :green)
    surface_color = :black
end

const k = 1024

const k_rays = 22

_rayplot(x...; kwargs...) = raytraceplot(raypoints(x...)...; kwargs...)

_rayplot!(x...; kwargs...) = raytraceplot!(raypoints(x...)...; kwargs...)

function _rayplot(lens::Lens, a::AbstractVector, h′ = -0.5; kwargs...)
    _rayplot(solve(lens, a, h′); kwargs...)
end

function _rayplot!(lens::Lens, a::AbstractVector, h′ = -0.5; kwargs...)
    _rayplot!(solve(lens, a, h′); kwargs...)
end

function _rayplot(
    surfaces::Matrix{Float64},
    a::AbstractVector,
    h′::Float64 = -0.5;
    kwargs...
)
    _rayplot(Lens(surfaces), a, h′; kwargs...)
end

function _rayplot!(
    surfaces::Matrix{Float64},
    a::AbstractVector,
    h′::Float64 = -0.5;
    kwargs...
)
    _rayplot!(Lens(surfaces), a, h′; kwargs...)
end

function _rayplot(
    surfaces::Matrix{Float64},
    a::AbstractVector,
    stop::Int,
    rays::RayBasis;
    kwargs...
)
    raytraceplot(surfaces, a, stop, raypoints(rays)...; kwargs...)
end

function _rayplot!(
    surfaces::Matrix{Float64},
    a::AbstractVector,
    stop::Int,
    rays::RayBasis;
    kwargs...
)
    raytraceplot!(surfaces, a, stop, raypoints(rays)...; kwargs...)
end

function _rayplot(
    surfaces::Matrix{Float64},
    system::System;
    kwargs...
)
    raytraceplot(surfaces, system.a, system.stop, raypoints(system)...; kwargs...)
end

function _rayplot!(
    surfaces::Matrix{Float64},
    system::System;
    kwargs...
)
    raytraceplot!(surfaces, system.a, system.stop, raypoints(system)...; kwargs...)
end

function wavefan(W::Aberration; k = k, kwargs...)
    fig = Figure()
    tangential_axis = Axis(fig[1,1];
        title = "Wave Aberration\nTangential",
        xlabel = L"y_p",
        ylabel = L"OPD"
    )
    sagittal_axis = Axis(fig[1,2];
        title = "Wave Aberration\nSagittal",
        xlabel = L"x_p",
        ylabel = L"OPD"
    )
    ρ = x = H = range(0.0, 1.0, k)
    grid = GridLayout(fig[1,3])
    grid[2,1] = slider = Slider(fig;
        horizontal = false,
        range = H,
        startvalue = 1.0,
        snap = false,
    )
    H = slider.value
    on(H) do _
        reset_limits!(tangential_axis)
        reset_limits!(sagittal_axis)
    end
    grid[1,1] = Label(fig, @lift("H: " * @sprintf("%.3f", $H)))
    y = range(-1.0, 1.0, k)
    Wy = @lift W.(ρ, 0.0, $H)
    Wx = @lift W.(ρ, π/2, $H)
    lines!(tangential_axis, y, Wy; color = :black, kwargs...)
    lines!(sagittal_axis, x, Wx; color = :black, kwargs...)
    DataInspector(fig)
    return fig
end

function rayfan(W::Aberration, s::SystemOrRayBasis, ; k = k, kwargs...)
    fig = Figure()
    tangential_axis = Axis(fig[1,1];
        title = "Transverse Ray Error\nTangential",
        xlabel = L"y_p",
        ylabel = L"\varepsilon_Y"
    )
    sagittal_axis = Axis(fig[1,2];
        title = "Transverse Ray Error\nSagittal",
        xlabel = L"x_p",
        ylabel = L"\varepsilon_X"
    )
    x = H = range(0.0, 1.0, k)
    grid = GridLayout(fig[1,3])
    grid[2,1] = slider = Slider(fig;
        horizontal = false,
        range = H,
        startvalue = 1.0,
        snap = false,
    )
    H = slider.value
    on(H) do _
        reset_limits!(tangential_axis)
        reset_limits!(sagittal_axis)
    end
    grid[1,1] = Label(fig, @lift("H: " * @sprintf("%.3f", $H)))
    y = range(-1.0, 1.0, k)
    ε_y = RayError{Tangential}(W, s)
    ε_x = RayError{Sagittal}(W, s)
    εy = @lift ε_y.(y, $H)
    εx = @lift ε_x.(x, $H)
    lines!(tangential_axis, y, εy; color = :black, kwargs...)
    lines!(sagittal_axis, x, εx; color = :black, kwargs...)
    DataInspector(fig)
    return fig
end

function field_curves(W::Aberration, s::SystemOrRayBasis; k = k, kwargs...)
    (; W220P, W220, W220T, W222, λ) = W
    nu = s.marginal.nu[end]
    u = s.marginal.u[end]
    α = -inv(nu * u) * λ
    fig = Figure()
    axis = Axis(fig[1,1];
        title = "Longitudinal Astigmatic Field Curves",
        xlabel = "z",
        ylabel = "H"
    )
    max_W220P = α * W220P
    max_W220T = α * W220T
    z = range(min(max_W220P, max_W220T), max(max_W220P, max_W220T), k)
    H = range(0.0, 1.0, k)
    P = max_W220P .* H .^ 2
    T = max_W220T .* H .^ 2
    S = α * W220 .* H .^ 2
    lines!(axis, P, H; label = "P", kwargs...)
    lines!(axis, T, H; label = "T", kwargs...)
    lines!(axis, S, H; label = "S", kwargs...)
    axislegend(axis)
    DataInspector(fig)
    return fig
end

function percent_distortion(W::Aberration, s::SystemOrRayBasis; k = k, kwargs...)
    (; W311, λ) = W
    ȳ = s.chief.y[end]
    n′u′ = s.marginal.nu[end]
    fig = Figure()
    axis = Axis(fig[1,1];
        title = "Percent Distortion",
        xlabel = "%",
        ylabel = "H"
    )
    p = range(0.0, 100.0, k)
    H = range(0.0, 1.0, k)
    pd = @. W311 * λ / n′u′ * H ^ 3 / ȳ * 100.0
    lines!(axis, pd, H; kwargs...)
    DataInspector(fig)
    return fig
end

function spot_size(W::Aberration, s::SystemOrRayBasis;
                   k = round(Int, √k) + 1, kwargs...)
    fig = Figure()
    axis = Axis(fig[1,1];
        xlabel = L"\varepsilon_X",
        ylabel = L"\varepsilon_Y"
    )
    H = range(0.0, 1.0, k)
    grid = GridLayout(fig[1,2])
    grid[2,1] = slider = Slider(fig;
        horizontal = false,
        range = H,
        startvalue = 1.0,
        snap = false,
    )
    H = slider.value
    on(_ -> reset_limits!(axis), H)
    grid[1,1] = Label(fig, @lift("H: " * @sprintf("%.3f", $H)))
    x = y = range(-1.0, 1.0, k)
    lattice = [(xᵢ, yᵢ) for xᵢ ∈ x for yᵢ ∈ y if hypot(xᵢ, yᵢ) ≤ 1.0]
    ε = RayError{Skew}(W, s)
    ε_x = RayError{Sagittal}(W, s)
    ε_y = RayError{Tangential}(W, s)
    ε_x_i = @lift ε_x.(x, $H)
    ε_y_i = @lift ε_y.(y, $H)
    εx_εy = @lift ε.(lattice, $H)
    mean_ε_x = @lift sum($ε_x_i) / k
    mean_ε_y = @lift sum($ε_y_i) / k
    var_x = @lift sum(($ε_x_i .- $mean_ε_x) .^ 2) / k
    var_y = @lift sum(($ε_y_i .- $mean_ε_y) .^ 2) / k
    RMS = @lift sqrt($var_x + $var_y)
    on(RMS; update = true) do RMS
        axis.title[] = "RMS Spot Size: " * @sprintf("%.5f", RMS)
    end
    scatter!(axis, εx_εy; kwargs...)
    DataInspector(fig)
    return fig
end

function caustic(surfaces::Matrix{Float64}, system::System,
                 k_rays::Int = k_rays; theme = default_plot_theme, kwargs...)
    with_theme(() -> raytraceplot(surfaces, system, k_rays; kwargs...), theme)
end

function caustic!(surfaces::Matrix{Float64}, system::System,
                  k_rays::Int = k_rays; theme = default_plot_theme, kwargs...)
    with_theme(() -> raytraceplot!(surfaces, system, k_rays; kwargs...), theme)
end

# plots ray bundles
function plot!(p::RayTracePlot{Tuple{T, Vector{T}}}) where T <: Vector{Float64}
    z, y = p.arg1[], p.arg2[]
    attr = p.attributes
    ray_colors = attr.ray_colors[]
    surface_color = attr.surface_color[]
    i = 0
    for yi in y
        i += 1
        ray_color = i > 3 ? ray_colors[2] : ray_colors[1]
        lines!(p, attr, z, yi; color = ray_color)
        scatter!(p, z, yi; color = ray_color)
    end
    zf = z[end]
    yf = y[end][end]
    lines!(p, attr, [zf for i = 1:2], [-yf, yf]; color = surface_color)
    DataInspector(textcolor = :black)
    return p
end

# plots lens surfaces + rays
function plot!(p::RayTracePlot{<:Tuple{Matrix{Float64}, <:AbstractVector, Int,
                               T, Vector{T}}}) where T <: Vector{Float64}
    surfaces, a, stop, z, y = p.arg1[], p.arg2[], p.arg3[], p.arg4[], p.arg5[]
    attr = p.attributes
    raytraceplot!(p, p.attributes, surfaces, a, stop, z)
    raytraceplot!(p, attr, z, y)
    return p
end

# plots lens surfaces
function plot!(p::RayTracePlot{<:Tuple{Matrix{Float64}, <:AbstractVector, Int,
                               Vector{Float64}}})
    surfaces, a, stop, z = p.arg1[], p.arg2[], p.arg3[], p.arg4[]
    attr = p.attributes
    surface_color = attr.surface_color[]
    n = @view surfaces[:,3]
    R = @view surfaces[2:end,1]
    z_surfaces = @view z[2:end-1]
    z_corners = Float64[]
    y_corners = Float64[]
    for (i, r) in pairs(R)
        z_vertex = z_surfaces[i]
        ai = a[i]
        z_max = r - sign(r) * sqrt(r ^ 2 - ai ^ 2) + z_vertex
        if isfinite(r)
            zi = range(z_vertex, z_max, 70)
            y_surface = @. sqrt(abs(r ^ 2 - (abs(zi - z_vertex) - abs(r)) ^ 2))
        else
            z_max = z_vertex
            zi = [z_vertex, z_vertex]
            y_surface = [0.0, ai]
        end
        linestyle = i == stop ? (; linestyle = :dot) : (; )
        lines!(p, attr, zi, y_surface; color = surface_color, linestyle...)
        lines!(p, attr, zi, -y_surface; color = surface_color, linestyle...)
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
    return p
end

# plots the caustic
function plot!(p::RayTracePlot{Tuple{Matrix{Float64}, System, Int}})
    surfaces, system, k_rays = p.arg1[], p.arg2[], p.arg3[]
    attr = p.attributes
    color = attr.ray_colors[][1]
    surface_color = attr.surface_color[]
    (; a, marginal, stop) = system
    (; z) = marginal
    raytraceplot!(p, p.attributes, surfaces, a, stop, z)
    real_marginal_ray = trace_marginal_ray(surfaces, system)
    a_stop = marginal.y[begin+stop]
    d = marginal.y[1] / k_rays
    yi = real_marginal_ray.y[1]
    paraxial_BFD = z[end] - z[end-1]
    marginal_BFD = -real_marginal_ray.y[end] / tan(real_marginal_ray.u[end])
    for i = 1:k_rays
        ray = raytrace(surfaces, yi, 0.0, RealRay)
        if marginal_BFD < paraxial_BFD
            # extend the rays for negative longitudinal spherical aberration
            zf = z[end]
            s = ray.z[end] - ray.z[end-1]
            yf = ray.y[end] + tan(ray.u[end]) * (paraxial_BFD + s)
        else
            # extend out to the marginal focus for positive LSA
            zf = real_marginal_ray.z[end-1] + marginal_BFD
            yf = ray.y[end] + tan(ray.u[end]) * marginal_BFD
        end
        # object space
        lines!(p, attr, [z[1], ray.z[1]], [ray.y[1]; ray.y[1]]; color)
        lines!(p, attr, [z[1], ray.z[1]], [-ray.y[1]; -ray.y[1]]; color)
        # surface trace
        lines!(p, attr, @view(ray.z[1:end-1]), @view(ray.y[2:end]); color)
        lines!(p, attr, @view(ray.z[1:end-1]), -@view(ray.y[2:end]); color)
        # paraxial focus
        lines!(p, attr, [z[end], z[end]], [-a_stop, a_stop]; color = surface_color)
        # image space
        lines!(p, attr, [ray.z[end-1], zf], [ray.y[end], yf]; color)
        lines!(p, attr, [ray.z[end-1], zf], [-ray.y[end], -yf]; color)
        yi -= d
    end
    # optical axis
    lines!(p, attr, [z[1], z[end]], [0.0, 0.0]; color = surface_color)
    DataInspector(textcolor = :black)
    return p
end

end
