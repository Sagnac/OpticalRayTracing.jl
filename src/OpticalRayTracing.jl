module OpticalRayTracing

using Printf

export TransferMatrix, Lens, construct, transfer, reverse_transfer, raytrace,
       trace_marginal_ray, trace_chief_ray, scale!, Rays, OpticalRays, ParaxialRays,
       solve, flatten, raypoints, rayplot

abstract type Rays end

struct TransferMatrix
    A::Matrix{Float64}
end

struct Lens
    A::Matrix{Float64}
end

function TransferMatrix(lens::Lens)
    (; A) = lens
    τ, ϕ = eachcol(A)
    A = prod([1.0 τ[i]; -ϕ[i] 1.0-τ[i]*ϕ[i]] for i = reverse(axes(A, 1)))
    TransferMatrix(A)
end

struct Pupil
    D::Float64
    t::Float64
end

struct OpticalRays <: Rays
    marginal::Matrix{Float64}
    chief::Matrix{Float64}
    H::Float64
end

struct ParaxialRays <: Rays
    marginal::Matrix{Float64}
    chief::Matrix{Float64}
    H::Float64
    n::Vector{Float64}
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
    rays::OpticalRays
    transfer_matrix::TransferMatrix
end

function ParaxialRays(rays::OpticalRays, n::Vector)
    (; marginal, chief, H) = rays
    marginal = copy(marginal)
    chief = copy(chief)
    @. marginal[:,2] /= n
    @. chief[:,2] /= n
    return ParaxialRays(marginal, chief, H, n)
end

function OpticalRays(rays::ParaxialRays)
    (; marginal, chief, H, n) = rays
    marginal = copy(marginal)
    chief = copy(chief)
    @. marginal[:,2] *= n
    @. chief[:,2] *= n
    return OpticalRays(marginal, chief, H)
end

function scale!(A::Matrix{Float64})
    A[:,2] .*= 1E-3
    return A
end

function construct(surfaces::Matrix{Float64})
    rows = size(surfaces, 1)
    R, t, n = eachcol(surfaces)
    A = Matrix{Float64}(undef, rows, 2)
    t[1] *= isfinite(t[1])
    @. A[:,1] = t / n
    for i = 1:rows-1
        A[i,2] = (n[i+1] - n[i]) / R[i+1]
    end
    if iszero(t[end]) || !isfinite(t[end])
        A = A[1:end-1,:]
    else
        A[end,2] = 0.0
    end
    return Lens(A)
end

function transfer(y, ω, τ, ϕ)
    y′ = isfinite(τ) ? ω * τ + y : y
    ω′ = ω - y′ * ϕ
    return y′, ω′
end

extend(A, τ, τ′) = [1.0 τ′; 0.0 1.0] * A * [1.0 τ; 0.0 1.0]

transfer(A::Matrix, v::Vector, τ, τ′) = extend(A, τ, τ′) * v
transfer(system::System, v::Vector, τ, τ′) = transfer(system.A, v, τ, τ′)
transfer(A::TransferMatrix, v::Vector, τ, τ′) = transfer(A.A, v, τ, τ′)

reverse_transfer(A::Matrix, v::Vector, τ′, τ) = extend(A, τ, τ′) \ v

function reverse_transfer(system::System, v::Vector, τ′, τ)
    reverse_transfer(system.A, v, τ′, τ)
end

function reverse_transfer(A::TransferMatrix, v::Vector, τ′, τ)
    reverse_transfer(A.A, v, τ′, τ)
end

function flatten(transfer_matrix::TransferMatrix)
    (; A) = transfer_matrix
    f = -inv(A[2,1])
    δ = (1.0 - A[2,2]) * f
    δ′ = (A[1,1] - 1.0) * f
    return (; f, δ, δ′)
end

function rev(objective, stop)
    v = @view (reverse ∘ transpose)(view(objective, 1:stop, :))[begin+1:end-1]
    Lens(reshape(v, 2, stop-1)')
end

function rev(objective_chief_ray)
    v = transpose(objective_chief_ray)[:]
    push!(v, 0.0)
    reverse!(v)
    push!(v, objective_chief_ray[1,2])
    reshape(v, 2, size(objective_chief_ray, 1) + 1)'
end

function raytrace(lens::Lens, y, ω, a = fill(Inf, size(lens.A, 1)); clip = false)
    (; A) = lens
    τ, ϕ = eachcol(A)
    rt = similar(A, size(A, 1) + 1, size(A, 2))
    rt[1,:] .= y, ω
    for i = axes(A, 1)
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
    raytrace(construct(surfaces), y, ω, a)
end

function trace_marginal_ray(lens::Lens, a, ω = 0.0)
    marginal_ray = raytrace(lens, 1.0, ω, a)
    y, ω = eachcol(marginal_ray)
    f = -inv(ω[end])
    EBFD = -y[end] / ω[end]
    sv = a ./ @view(y[begin+1:end])
    s, stop = findmin(sv)
    marginal_ray *= s
    marginal_ray = [marginal_ray; [0.0 marginal_ray[end,2]]]
    return marginal_ray, stop, f, EBFD
end

function trace_chief_ray(lens::Lens, stop, EBFD, h′ = -0.5)
    (; A) = lens
    A = [A; [EBFD 0.0]]
    rear_lens = Lens(@view(A[stop+1:end,:]))
    objective = rev(A, stop)
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
    FOV = 2atand(abs(h′ / f))
    return System(f, EBFD, EFFD, N, FOV, stop, EP, XP,
                  OpticalRays(marginal_ray, chief_ray, H),
                  TransferMatrix(lens))
end

solve(surfaces, a, h′ = -0.5) = solve(construct(surfaces), a, h′)

function Base.show(io::IO, system::T) where T <: System
    print(io, "f: ")
    show(IOContext(io, :compact => true), system.f)
end

function Base.show(io::IO, ::MIME"text/plain", system::T) where T <: System
    if haskey(io, :typeinfo)
        show(io, system)
        return
    end
    for property in fieldnames(T)
        property === :rays && break
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

function Base.getproperty(system::System, property::Symbol)
    if property === :A
        getfield(system, :transfer_matrix).A
    else
        getfield(system, property)
    end
end

function raypoints(n::AbstractVector, lens::Lens, system::System)
    (; A) = lens
    (; marginal, chief) = system.rays
    t = map(*, @view(A[:,1]), @view(n[begin:end-1]))
    k = length(t) + 2
    l = sum(t)
    s = 0.1
    d = -l * s / n[1]
    at_objective = iszero(t[1])
    z = Vector{Float64}(undef, k)
    z[1] = at_objective ? d : 0.0
    z[2] = at_objective ? 0.0 : t[1]
    for i = 2:lastindex(t)
        z[i+1] = z[i] + t[i]
    end
    z[end] = z[end-1] + system.EBFD * n[end]
    y0 = zeros(k)
    y1 = @view marginal[:,1]
    y2 = -y1
    nū = chief[1,2]
    ȳ1 = chief[2,1]
    y11 = y1[1]
    y21 = y2[1]
    ȳ2 = ȳ1 + y11
    ȳ3 = ȳ1 + y21
    ȳo1 = ȳ1 + nū * d
    ȳo2 = ȳo1 + y11
    ȳo3 = ȳo1 + y21
    ȳ = [ȳo1; @view(chief[begin+1:end,1])]
    ȳ′ = chief[end,1]
    y3 = [@view(raytrace(lens, ȳ2, nū)[:,1]); ȳ′]
    y4 = [@view(raytrace(lens, ȳ3, nū)[:,1]); ȳ′]
    y3[1] = ȳo2
    y4[1] = ȳo3
    return z, y0, y1, y2, ȳ, y3, y4
end

# requires Makie
# TODO: create a recipe
function rayplot(n::AbstractVector, lens::Lens, system::System)
    fig = Main.Figure()
    axis = Main.Axis(fig[1,1]; xlabel = "z", ylabel = "y")
    z, y... = raypoints(n, lens, system)
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

function rayplot(n::AbstractVector, lens::Lens, a::AbstractVector, h′ = -0.5)
    system = solve(lens, a, h′)
    rayplot(n, lens, system)
end

function rayplot(surfaces::Matrix{Float64}, system::System)
    n = @view surfaces[:,3]
    lens = construct(surfaces)
    rayplot(n, lens, system)
end

function rayplot(surfaces::Matrix{Float64}, a::AbstractVector, h′ = -0.5)
    n = @view surfaces[:,3]
    lens = construct(surfaces)
    system = solve(lens, a, h′)
    rayplot(n, lens, system)
end

end
