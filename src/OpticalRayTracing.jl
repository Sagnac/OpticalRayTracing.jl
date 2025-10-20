module OpticalRayTracing

using Printf

export TransferMatrix, Lens, transfer, reverse_transfer, raytrace,
       trace_marginal_ray, trace_chief_ray, scale!, Ray, Marginal, Chief,
       solve, flatten, raypoints, rayplot, vignetting, aberrations

abstract type ParaxialRay end

abstract type TangentialRay <: ParaxialRay end

struct Marginal <: TangentialRay end

struct Chief <: TangentialRay end

struct Ray{T <: TangentialRay}
    y::Vector{Float64}
    n::Vector{Float64}
    u::Vector{Float64}
    yu::Matrix{Float64}
    nu::Vector{Float64}
    ynu::Matrix{Float64}
    function Ray{T}(ynu, n) where T <: TangentialRay
        y, nu = eachcol(ynu)
        n = [n; n[end]]
        u = map(/, nu, n)
        yu = [y u]
        new(y, n, u, yu, nu, ynu)
    end
end

struct TransferMatrix
    M::Matrix{Float64}
end

struct Lens
    M::Matrix{Float64}
    n::Vector{Float64}
end

struct Aberration
    spherical::Float64
    coma::Float64
    astigmatism::Float64
    petzval::Float64
    distortion::Float64
    axial::Float64
    lateral::Float64
    sagittal::Float64
    medial::Float64
    tangential::Float64
    W040::Vector{Float64}
    W131::Vector{Float64}
    W222::Vector{Float64}
    W220P::Vector{Float64}
    W220::Vector{Float64}
    W311::Vector{Float64}
    W020::Vector{Float64}
    W111::Vector{Float64}
end

function TransferMatrix(lens::Lens)
    (; M) = lens
    τ, ϕ = eachcol(M)
    M = prod([1.0 τ[i]; -ϕ[i] 1.0-τ[i]*ϕ[i]] for i = reverse(axes(M, 1)))
    TransferMatrix(M)
end

struct Pupil
    D::Float64
    t::Float64
end

struct System
    f::Float64
    EBFD::Float64
    EFFD::Float64
    N::Float64
    FOV::Float64
    stop::Int
    EP::Pupil
    XP::Pupil
    marginal::Ray{Marginal}
    chief::Ray{Chief}
    H::Float64
    M::TransferMatrix
    lens::Lens
end

function scale!(M::Matrix{Float64})
    M[:,2] .*= 1E-3
    return M
end

function Lens(surfaces::Matrix{Float64})
    rows = size(surfaces, 1)
    R, t, n = eachcol(surfaces)
    M = Matrix{Float64}(undef, rows, 2)
    t[1] *= isfinite(t[1])
    @. M[:,1] = t / n
    for i = 1:rows-1
        M[i,2] = (n[i+1] - n[i]) / R[i+1]
    end
    if iszero(t[end]) || !isfinite(t[end])
        M = M[1:end-1,:]
    else
        M[end,2] = 0.0
    end
    return Lens(M, n)
end

function transfer(y, ω, τ, ϕ)
    y′ = isfinite(τ) ? ω * τ + y : y
    ω′ = ω - y′ * ϕ
    return y′, ω′
end

extend(M, τ, τ′) = [1.0 τ′; 0.0 1.0] * M * [1.0 τ; 0.0 1.0]

transfer(M::Matrix, v::Vector, τ, τ′) = extend(M, τ, τ′) * v
transfer(system::System, v::Vector, τ, τ′) = transfer(system.M, v, τ, τ′)
transfer(M::TransferMatrix, v::Vector, τ, τ′) = transfer(M.M, v, τ, τ′)

reverse_transfer(M::Matrix, v::Vector, τ′, τ) = extend(M, τ, τ′) \ v

function reverse_transfer(system::System, v::Vector, τ′, τ)
    reverse_transfer(system.M, v, τ′, τ)
end

function reverse_transfer(M::TransferMatrix, v::Vector, τ′, τ)
    reverse_transfer(M.M, v, τ′, τ)
end

flatten(transfer_matrix::TransferMatrix) = flatten(transfer_matrix.M)

function flatten(M::Matrix)
    f = -inv(M[2,1])
    EFFD = -M[2,2] * f
    EBFD = M[1,1] * f
    return (; f, EFFD, EBFD)
end

function rev(objective, stop, n)
    v = @view (reverse ∘ transpose)(view(objective, 1:stop, :))[begin+1:end-1]
    Lens(reshape(v, 2, stop-1)', n)
end

function rev(objective_chief_ray)
    v = transpose(objective_chief_ray)[:]
    push!(v, 0.0)
    reverse!(v)
    push!(v, objective_chief_ray[1,2])
    reshape(v, 2, size(objective_chief_ray, 1) + 1)'
end

function raytrace(lens::Lens, y, ω, a = fill(Inf, size(lens.M, 1)); clip = false)
    (; M) = lens
    τ, ϕ = eachcol(M)
    rt = similar(M, size(M, 1) + 1, size(M, 2))
    rt[1,:] .= y, ω
    for i = axes(M, 1)
        y, ω = transfer(y, ω, τ[i], ϕ[i])
        if clip && y > a[i]
            rt[i+1:end,:] .= NaN
            break
        end
        rt[i+1,1] = y
        rt[i+1,2] = ω
    end
    return rt
end

function raytrace(surfaces::Matrix, y, ω, a = fill(Inf, size(surfaces, 1)))
    raytrace(Lens(surfaces), y, ω, a)
end

function trace_marginal_ray(lens::Lens, a, ω = 0.0)
    marginal_ray = raytrace(lens, 1.0, ω, a)
    y, ω = eachcol(marginal_ray)
    f = -inv(ω[end])
    EBFD = y[end] * f
    sv = a ./ @view(y[begin+1:end])
    s, stop = findmin(sv)
    marginal_ray *= s
    ωf = marginal_ray[end,2]
    yf = iszero(ωf) ? marginal_ray[end,1] : 0.0
    marginal_ray = [marginal_ray; [yf ωf]]
    return marginal_ray, stop, f, EBFD
end

function trace_chief_ray(lens::Lens, stop, EBFD, h′ = -0.5)
    (; M, n) = lens
    M = [M; [EBFD 0.0]]
    rear_lens = Lens(M[stop+1:end,:], n)
    objective = rev(M, stop, n)
    rear_chief_ray = raytrace(rear_lens, 0.0, -1.0)
    ȳ′ = rear_chief_ray[end,1]
    s = h′ / ȳ′
    rear_chief_ray *= s
    objective_chief_ray = rev(raytrace(objective, 0.0, s))
    objective_chief_ray[:,2] .*= -1.0
    chief_ray = [objective_chief_ray; @view(rear_chief_ray[begin+1:end,:])]
    return chief_ray
end

function solve(lens::Lens, a::AbstractVector, h′::Float64 = -0.5)
    marginal_ray, stop, f, EBFD = trace_marginal_ray(lens, a)
    chief_ray = trace_chief_ray(lens, stop, EBFD, h′)
    ȳ = chief_ray[begin+1,1]
    nū = chief_ray[begin,2]
    ȳ′ = chief_ray[end,1]
    n′ū′ = chief_ray[end,2]
    y = marginal_ray[begin,1]
    nu = marginal_ray[begin,2]
    y′ = marginal_ray[end,1]
    n′u′ = marginal_ray[end,2]
    ȳ′b = chief_ray[end-1,1]
    # δ′ = EBFD - f
    δ = (ȳ′ - n′ū′ * f - ȳ) / nū
    EFFD = δ - f
    EP = Pupil(abs(y) * 2, -ȳ / nū)
    H = nū * y
    XP = Pupil(abs(2H / n′ū′), -ȳ′b / n′ū′)
    N = abs(f / EP.D)
    n = lens.n
    FOV = 2atand(abs(nū / n[1]))
    return System(f, EBFD, EFFD, N, FOV, stop, EP, XP,
                  Ray{Marginal}(marginal_ray, n), Ray{Chief}(chief_ray, n), H,
                  TransferMatrix(lens), lens)
end

solve(surfaces, a, h′ = -0.5) = solve(Lens(surfaces), a, h′)

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

# Helium d-line wavelength in millimeters
const λ = 587.5618e-6

Δ(x, y, i) = x[i+1] / y[i+1] - x[i] / y[i]

function aberrations(surfaces::Matrix{Float64}, system::System,
                     δn::Vector{Float64} = zeros(size(surfaces, 1)),
                     λ::Float64 = λ)
    (; marginal, chief, H) = system
    R = @view surfaces[2:end,1]
    n = marginal.n
    y = @view marginal.y[2:end-1]
    ȳ = @view chief.y[2:end-1]
    nu = marginal.nu
    nū = chief.nu
    u = marginal.u
    j = eachindex(R)
    A = [nu[i] + n[i] * y[i] / R[i] for i ∈ j]
    Ā = [(H + A[i] * ȳ[i]) / y[i] for i ∈ j]
    yΔ = [y[i] * Δ(u, n, i) for i ∈ j]
    yδ = [y[i] * Δ(δn, n, i) for i ∈ j]
    P = [(inv(n[i+1]) - inv(n[i])) / R[i] for i ∈ j]
    # per surface contributions
    W040 = @. -A ^ 2 * yΔ / 8λ
    W131 = @. -A * Ā * yΔ / 2λ
    W222 = @. -Ā ^ 2 * yΔ / 2λ
    W220P = @. -H ^ 2 * P / 4λ
    W220 = @. W220P + W222 / 2
    W311 = @. Ā / A * (W222 + 2W220P)
    W020 = @. A * yδ / 2λ
    W111 = @. Ā * yδ / λ
    # system coefficients
    spherical = sum(W040)
    coma = sum(W131)
    astigmatism = sum(W222)
    distortion = sum(W311)
    # field curvatures
    petzval = sum(W220P)
    sagittal = petzval + 0.5 * astigmatism
    medial = petzval + astigmatism
    tangential = petzval + 1.5 * astigmatism
    # chromatic aberrations
    # defocus from longitudinal
    axial = sum(W020)
    # tilt from transverse
    lateral = sum(W111)
    Aberration(spherical, coma, astigmatism, petzval, distortion, axial, lateral,
               sagittal, medial, tangential,
               W040, W131, W222, W220P, W220, W311, W020, W111)
end

function Base.show(io::IO, system::T) where T <: System
    print(io, "f: ")
    show(IOContext(io, :compact => true), system.f)
end

function Base.show(io::IO, ::MIME"text/plain", system::T) where T <: System
    println(T)
    if haskey(io, :typeinfo)
        show(io, system)
        return
    end
    for property in fieldnames(T)
        property === :marginal && break
        value = getproperty(system, property)
        @printf(io, "\n%4s: ", property)
        if typeof(value) === Pupil
            @printf(io, "D = %.4f, t = %.4f", value.D, value.t)
        elseif property === :stop
            @printf(io, "%d", value)
        else
            @printf(io, "%.4f", value)
        end
    end
    return
end

function Base.show(io::IO, m::MIME"text/plain", ray::Ray)
    summary(io, ray)
    println(io, ".yu:")
    show(IOContext(io, :displaysize => displaysize(io) .- (1, 0)), m, ray.yu)
end

function Base.show(io::IO, ::MIME"text/plain", aberr::T) where T <: Aberration
    println(T)
    for property in fieldnames(T)
        property === :lateral && break
        @printf(io, "\n%4s: %.4f", property, getproperty(aberr, property))
    end
    return
end

Base.getindex(M::TransferMatrix) = M.M

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
# TODO: create a recipe
function rayplot(system::System)
    fig, axis, z, y = _rayplot(system)
    _rayplot!(fig, axis, z, y)
end

function _rayplot(system::System)
    fig = Main.Figure()
    axis = Main.Axis(fig[1,1])
    z, y... = raypoints(system)
    return fig, axis, z, y
end

function _rayplot!(fig, axis, z, y)
    i = 0
    for yi in y
        i += 1
        color = i > 3 ? :green : :blue
        Main.scatter!(axis, z, yi; color)
        Main.lines!(axis, z, yi; color)
    end
    zf = z[end]
    yf = y[4][end]
    Main.lines!(axis, [zf for i = 1:2], [-yf, yf]; color = :black)
    Main.DataInspector(fig)
    return fig
end

function rayplot(lens::Lens, a::AbstractVector, h′ = -0.5)
    system = solve(lens, a, h′)
    rayplot(system)
end

function rayplot(surfaces::Matrix{Float64}, a::AbstractVector, h′ = -0.5)
    lens = Lens(surfaces)
    rayplot(lens, a, h′)
end

function rayplot(surfaces::Matrix{Float64}, system::System, a::AbstractVector)
    fig, axis, z, y = _rayplot(system)
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
        Main.lines!(axis, zi, y_surface; color = :black)
        Main.lines!(axis, zi, -y_surface; color = :black)
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
        Main.lines!(axis, z_edge, y_edge; color = :black)
        Main.lines!(axis, z_edge, -y_edge; color = :black)
    end
    _rayplot!(fig, axis, z, y)
    return fig
end

end
