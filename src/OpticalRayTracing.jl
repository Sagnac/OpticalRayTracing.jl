module OpticalRayTracing

using Printf

export TransferMatrix, Lens, transfer, reverse_transfer, raytrace,
       trace_marginal_ray, trace_chief_ray, scale!, Ray, Marginal, Chief,
       solve, flatten, raypoints, rayplot, rayplot!, vignetting, aberrations

include("Types.jl")
include("TransferMatrix.jl")
include("Vignetting.jl")
include("SeidelAberrations.jl")
include("BaseMethods.jl")
include("RayPlot.jl")

function scale!(M::Matrix{Float64})
    M[:,2] .*= 1E-3
    return M
end

reduced_thickness(lens::Lens) = lens.M[:,1]

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
    y′ = transfer(y, ω, τ)
    ω′ = refract(y′, ω, ϕ)
    return y′, ω′
end

function transfer(y, ω, τ)
    y′ = isfinite(τ) ? y + ω * τ : y
    return y′
end

function refract(y, ω, ϕ)
    ω′ = ω - y * ϕ
    return ω′
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
    (; M, n) = lens
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
    return Ray{TangentialRay}(rt, τ, n)
end

function raytrace(surfaces::Matrix, y, ω, a = fill(Inf, size(surfaces, 1)))
    raytrace(Lens(surfaces), y, ω, a)
end

function raytrace(system::System, ȳ, s,
                  a = fill(Inf, size(system.lens.M, 1)); clip = false)
    (; lens, EP) = system
    (; M) = lens
    n = lens.n
    t = -EP.t
    EP_O = s + t
    y = EP.D / 2
    nu = -y / EP_O
    nū = ȳ / EP_O
    lens = Lens([[t M[1,2]]; @view(M[2:end,:])], n)
    marginal_ray_rt = raytrace(lens, y, nu, a; clip)
    marginal_ray = extend(marginal_ray_rt.ynu)
    chief_ray_rt = raytrace(lens, 0.0, nū, a; clip)
    ȳ′ = -nū * y / marginal_ray_rt.nu[end]
    chief_ray = chief_ray_rt.ynu
    n′ū′ = chief_ray_rt.nu[end]
    chief_ray = [chief_ray; [ȳ′ n′ū′]]
    τ = reduced_thickness(lens)
    return Ray{Marginal}(marginal_ray, τ, n), Ray{Chief}(chief_ray, τ, n)
end

function extend(marginal_ray)
    ωf = marginal_ray[end,2]
    yf = iszero(ωf) ? marginal_ray[end,1] : 0.0
    return [marginal_ray; [yf ωf]]
end

function trace_marginal_ray(lens::Lens, a, ω = 0.0)
    rt = raytrace(lens, 1.0, ω, a)
    marginal_ray = rt.ynu
    y = rt.y
    ω = rt.nu
    f = -inv(ω[end])
    EBFD = y[end] * f
    sv = a ./ @view(y[begin+1:end])
    s, stop = findmin(sv)
    marginal_ray *= s
    marginal_ray = extend(marginal_ray)
    τ = reduced_thickness(lens)
    return Ray{Marginal}(marginal_ray, τ, lens.n), stop, f, EBFD
end

function trace_chief_ray(lens::Lens, stop, EBFD, h′ = -0.5)
    (; M, n) = lens
    M = [M; [EBFD 0.0]]
    rear_lens = Lens(M[stop+1:end,:], n)
    objective = rev(M, stop, n)
    rt = raytrace(rear_lens, 0.0, -1.0)
    rear_chief_ray = rt.ynu
    ȳ′ = rt.y[end]
    s = h′ / ȳ′
    rear_chief_ray *= s
    objective_chief_ray = rev(raytrace(objective, 0.0, s).ynu)
    objective_chief_ray[:,2] .*= -1.0
    chief_ray = [objective_chief_ray; @view(rear_chief_ray[begin+1:end,:])]
    τ = reduced_thickness(lens)
    return Ray{Chief}(chief_ray, τ, n)
end

function solve(lens::Lens, a::AbstractVector, h′::Float64 = -0.5)
    marginal_ray, stop, f, EBFD = trace_marginal_ray(lens, a)
    chief_ray = trace_chief_ray(lens, stop, EBFD, h′)
    ȳ = chief_ray.y[begin+1]
    nū = chief_ray.nu[begin]
    ȳ′ = chief_ray.y[end]
    n′ū′ = chief_ray.nu[end]
    y = marginal_ray.y[begin]
    nu = marginal_ray.nu[begin]
    y′ = marginal_ray.y[end]
    n′u′ = marginal_ray.nu[end]
    ȳ′b = chief_ray.y[end-1]
    # δ′ = EBFD - f
    δ = (ȳ′ - n′ū′ * f - ȳ) / nū
    EFFD = δ - f
    EP = Pupil(abs(y) * 2, -ȳ / nū)
    H = nū * y
    XP = Pupil(abs(2H / n′ū′), -ȳ′b / n′ū′)
    N = abs(f / EP.D)
    FOV = 2atand(abs(chief_ray.u[1]))
    return System(f, EBFD, EFFD, N, FOV, stop, EP, XP,
                  marginal_ray, chief_ray, H,
                  TransferMatrix(lens), lens)
end

solve(surfaces, a, h′ = -0.5) = solve(Lens(surfaces), a, h′)

end
