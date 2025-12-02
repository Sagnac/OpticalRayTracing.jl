const ϵ = sqrt(eps())

const k_rays = 22

function scale!(M::Lens)
    M[:,2] .*= 1E-3
    return M
end

reduced_thickness(lens::Lens) = lens[:,1]

function compute_surfaces(lens::Lens)
    (; M, n) = lens
    τ, ϕ = eachcol(M)
    k = length(ϕ)
    surfaces = Matrix{Float64}(undef, k + 1, 3)
    surfaces[1,:] = [Inf 0.0 n[1]]
    for i in 2:k
        nᵢ = n[i]
        ϕᵢ = ϕ[i-1]
        R = iszero(ϕᵢ) ? Inf : (nᵢ - n[i-1]) / ϕᵢ
        t = τ[i] * nᵢ
        surfaces[i,:] = [R t nᵢ]
    end
    nf = n[end]
    surfaces[end,:] = [(n[end] - n[end-1]) / ϕ[end] 0.0 nf]
    return surfaces
end

surface_ray(v::AbstractVector) = @view v[begin+1:end-1]

surface_ray(M::AbstractMatrix) = @view M[begin+1:end-1,:]

function Lens(surfaces::AbstractMatrix)
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

sag(y, R, ::Type{ParaxialRay}) = y ^ 2 / 2R

sag(y, R, ::Type{RealRay}) = R - sign(R) * sqrt(R ^ 2 - y ^ 2)

function sag(y, U, R, K, p)
    if isfinite(R)
        β = R - y * tan(U)
        y2 = y ^ 2
        Δ = β ^ 2 - y2 * (sec(U) ^ 2 + K)
        if Δ ≥ 0.0
            return y2 / (β + sign(R) * sqrt(Δ)) + p(y)
        else
            return NaN # misses surface
        end
    else
        return 0.0
    end
end

# z-coordinate end-point is vertex position
sag(ray::RealRay{Tangential}) = ray.z[end-1] - ray.z[end]

function sag(real::RealRay, paraxial::ParaxialRay{Marginal})
    real.z[end-1] - paraxial.z[end-1]
end

# surface normal slope
tilt(y, R, K, p) = sign(R) * y / sqrt(R ^ 2 - y ^ 2 * (1 + K)) + dp_dy(p, y)

# equivalent to the sine of the above tangent angle for spherical surfaces
tilt(y, R) = y / R

dp_dy(p, y) = imag(p(complex(y, ϵ))) / ϵ

surface_to_focus(BFD, x...) = BFD - sag(x...)

function transfer(ray::RealRay{<:FundamentalRay}, t)
    transfer(ray.y[end-1], ray.u[end-1], t, RealRay)
end

function transfer(ray::RealRay{Tangential}, t)
    transfer(ray.y[end], ray.u[end], t, RealRay)
end

transfer(y, u, t, ::Type{RealRay}) = y + tan(u) * t

function stop_loss(surfaces, y, u, stop, a_stop, ::Type{Marginal})
    ray = raytrace(surfaces, y, u, RealRay)
    return ray, ray.y[begin+stop] - a_stop
end

function stop_loss(surfaces, ȳ′, ū′, t, stop, ::Type{Chief})
    ȳ′ = transfer(ȳ′, ū′, t, RealRay)
    ray = raytrace(surfaces, ȳ′, ū′, RealRay)
    return ray, ray.y[begin+stop]
end

function raytrace(lens::Lens, y, ω,
                  a::AbstractVector = fill(Inf, size(lens, 1)); clip = false)
    (; M, n) = lens
    τ, ϕ = eachcol(M)
    rt = similar(M, size(M, 1) + 1, 2)
    rt[1,:] .= y, ω
    for i = axes(M, 1)
        y, ω = transfer(y, ω, τ[i], ϕ[i])
        if clip && abs(y) - a[i] > 1e-13
            rt[i+1:end,:] .= NaN
            break
        end
        rt[i+1,1] = y
        rt[i+1,2] = ω
    end
    return ParaxialRay{Tangential}(rt, τ, n)
end

function raytrace(surfaces::AbstractMatrix, y, U, ::Type{RealRay};
                  K = zeros(size(surfaces, 1)), p = fill(zero, length(K)))
    R, t, n = eachcol(surfaces)
    ts = copy(t)
    rt = similar(surfaces, size(surfaces, 1), 2)
    rt[1,:] .= y, U
    for i = 1:size(surfaces, 1)-1
        y += tan(U) * ts[i]
        Rs = R[i+1]
        Ks = K[i+1]
        ps = p[i+1]
        s = sag(y, U, Rs, Ks, ps)
        # transfer across surface sagitta
        y += s * tan(U)
        # update distances
        ts[i] += s
        ts[i+1] -= s
        θ = iszero(Ks) && ps ≡ zero ? asin(tilt(y, Rs)) : atan(tilt(y, Rs, Ks, ps))
        sin_i′ = n[i] * sin(U + θ) / n[i+1]
        U = abs(sin_i′) ≤ 1.0 ? asin(sin_i′) - θ : NaN # total internal reflection
        rt[i+1,1] = y
        rt[i+1,2] = U
    end
    return RealRay{Tangential}(rt, ts, n)
end

function raytrace(surfaces::Layout{Aspheric}, y, U, ::Type{RealRay})
    raytrace(surfaces, y, U, RealRay; K = surfaces.K, p = surfaces.p)
end

function raytrace(surfaces::AbstractMatrix, y, ω,
                  a::AbstractVector = fill(Inf, size(surfaces, 1)); clip = false)
    raytrace(Lens(surfaces), y, ω, a; clip)
end

function raytrace(system::System, ȳ, s)
    (; EP, chief, marginal, H, lens) = system
    y = marginal.y[1]
    EP_O = (s - EP.t)
    nu = -y / EP_O
    nū = ȳ / EP_O
    α = y * nu / H
    β = y * nū / H
    marginal_ray = marginal.ynu + α * chief.ynu
    chief_ray = β * chief.ynu
    H = nū * y
    marginal_ray[end,1] = 0.0
    chief_ray[end,1] = -H / marginal_ray[end,2]
    τ = reduced_thickness(lens)
    rays = RayBasis(
        ParaxialRay{Marginal}(marginal_ray, τ, lens.n),
        ParaxialRay{Chief}(chief_ray, τ, lens.n),
        H
    )
    return rays
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
    return ParaxialRay{Marginal}(marginal_ray, τ, lens.n), stop, f, EBFD
end

function trace_marginal_ray(surfaces, system::System; atol = sqrt(eps()))
    (; marginal, a, stop) = system
    y = marginal.y[1]
    u = 0.0
    a_stop = system.a[stop]
    ray, Δy_stop = stop_loss(surfaces, y, u, stop, a_stop, Marginal)
    while abs(Δy_stop) > atol
        δy = stop_loss(surfaces, y + ϵ, u, stop, a_stop, Marginal)[2]
        y -= Δy_stop * ϵ / (δy - Δy_stop)
        ray, Δy_stop = stop_loss(surfaces, y, u, stop, a_stop, Marginal)
    end
    ray.z[end] = ray.z[end-1] - ray.y[end] / tan(ray.u[end])
    pushfirst!(ray.z, -(extrema(ray.z)...) * 0.1)
    push!(ray.y, 0.0)
    push!(ray.u, ray.u[end])
    yu = vcat(ray.yu, hcat(0.0, ray.u[end]))
    return RealRay{Marginal}(ray.y, ray.u, yu, ray.n, ray.z)
end

function trace_chief_ray(lens::Lens, stop::Int,
                         marginal::ParaxialRay{Marginal}, h′ = -0.5)
    (; n) = lens
    y = surface_ray(marginal.y)
    ynu = surface_ray(marginal.ynu)
    y_stop = y[stop]
    rt = raytrace(lens, 0.0, 1.0)
    y2 = @view rt.y[2:end]
    ynu2 = @view rt.ynu[2:end,:]
    y2_stop = y2[stop]
    nū = -marginal.nu[end] * h′ / y[1]
    chief_ray = Matrix{Float64}(undef, size(marginal.ynu))
    @. chief_ray[begin+1:end-1,:] = nū * (ynu2 - ynu * y2_stop / y_stop)
    chief_ray[begin,:] .= 0.0, nū
    chief_ray[end,:] .= h′, chief_ray[end-1,2]
    τ = reduced_thickness(lens)
    return ParaxialRay{Chief}(chief_ray, τ, n)
end

function trace_chief_ray(surfaces, system::System; atol = sqrt(eps()))
    rev_R = -surfaces[end:-1:2, 1]
    rev_t = @view surfaces[end:-1:1, 2]
    rev_n = @view surfaces[end:-1:1, 3]
    rev_surfaces = [[Inf; rev_R] rev_t rev_n]
    (; chief, marginal) = system
    stop = length(rev_R) - system.stop + 1
    ȳ′ = chief.y[end]
    ū′ = -chief.u[end]
    BFD = marginal.z[end] - marginal.z[end-1]
    ray, y_stop = stop_loss(rev_surfaces, ȳ′, ū′, BFD, stop, Chief)
    while abs(y_stop) > atol
        δy = stop_loss(rev_surfaces, ȳ′, ū′ + ϵ, BFD, stop, Chief)[2]
        ū′ -= y_stop * ϵ / (δy - y_stop)
        ray, y_stop = stop_loss(rev_surfaces, ȳ′, ū′, BFD, stop, Chief)
    end
    ȳ = [0.0; reverse(ray.y)]
    ȳ[end] = ȳ′
    ū = [-reverse(ray.u); -ray.u[1]]
    ȳū = [ȳ ū]
    n = surfaces[:,3]
    z = ray.z[end] .- reverse(ray.z)
    z[1] = -ȳ[2] / tan(ū[1]) + z[2] # EP distance from vertex
    push!(z, z[end] - ȳ[end-1] / tan(ū[end-1]))
    return RealRay{Chief}(ȳ, ū, ȳū, n, z)
end

function solve(lens::Lens, a::AbstractVector, h′::Float64 = -0.5)
    marginal_ray, stop, f, EBFD = trace_marginal_ray(lens, a)
    chief_ray = trace_chief_ray(lens, stop, marginal_ray, h′)
    ȳ = chief_ray.y[begin+1]
    nū = chief_ray.nu[begin]
    ȳ′ = h′
    n′ū′ = chief_ray.nu[end]
    y = marginal_ray.y[begin]
    ȳ′b = chief_ray.y[end-1]
    δ′ = EBFD - f
    δ = (ȳ′ - n′ū′ * f - ȳ) / nū
    EFFD = δ - f
    PN = (lens.n[end] - lens.n[begin]) * f
    EP = Pupil(abs(y) * 2, -ȳ / nū)
    H = nū * y
    XP = Pupil(abs(2H / n′ū′), -ȳ′b / n′ū′)
    N = abs(f / EP.D)
    FOV = 2atand(abs(chief_ray.u[1]))
    return System(f, EBFD, EFFD, N, FOV, stop, EP, XP,
                  marginal_ray, chief_ray, [marginal_ray.yu chief_ray.yu], H,
                  δ, δ′, PN, TransferMatrix(lens), lens, a)
end

solve(surfaces, a, h′ = -0.5) = solve(Lens(surfaces), a, h′)

# paraxial incidence angles
function incidences(surfaces::AbstractMatrix, system::SystemOrRayBasis)
    R = @view surfaces[2:end,1]
    (; marginal, chief) = system
    (; n) = marginal
    nu, y = marginal.nu, surface_ray(marginal.y)
    nū, ȳ = chief.nu, surface_ray(chief.y)
    ni = map(nu, n, y, R) do nu, n, y, R
             return nu + n * y / R
         end
    nī = map(nū, n, ȳ, R) do nū, n, ȳ, R
             return nū + n * ȳ / R
         end
    i = map(/, ni, n)
    ī = map(/, nī, n)
    return [ni nī i ī]
end
