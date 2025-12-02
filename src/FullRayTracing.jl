function sag(y, x, u, v, R, K, p)
    if isfinite(R)
        β = R - y * u - x * v
        r2 = x ^ 2 + y ^ 2
        Δ = β ^ 2 - r2 * (1 + K + u ^ 2 + v ^ 2)
        if Δ ≥ 0.0
            return r2 / (β + sign(R) * sqrt(Δ)) + p(y)
        else
            return NaN # misses surface
        end
    else
        return 0.0
    end
end

function tilt(y, x, R, K, p)
    Δ = R ^ 2 - (x ^ 2 + y ^ 2) * (1 + K)
    [sign(R) * α / sqrt(Δ) + dp_dy(p, α) for α ∈ (x, y)]
end

function refract!(k::Vector, m::Vector, n1::Float64, n2::Float64)
    η = n1 / n2
    γ = -k ⋅ m
    Δ = 1 - η ^ 2 * (1 - γ ^ 2)
    if Δ ≥ 0.0
        @. k = η * k + (η * γ - sqrt(Δ)) * m
    else
        @. k += 2.0 * γ * m # TIR
    end
    return k
end

function raytrace(surfaces::AbstractMatrix, y, x, U, V, ::Type{Vector{RealRay}};
                  K = zeros(size(surfaces, 1)), p = fill(zero, length(K)))
    R, t, n = eachcol(surfaces)
    ts = copy(t)
    u = tan(U)
    v = tan(V)
    k = [v, u, 1.0]
    normalize!(k) # direction cosines
    for i = 1:size(surfaces, 1)-1
        y += u * ts[i]
        x += v * ts[i]
        Rs = R[i+1]
        Ks = K[i+1]
        ps = p[i+1]
        s = sag(y, x, u, v, Rs, Ks, ps)
        y += s * u
        x += s * v
        ts[i] += s
        ts[i+1] -= s
        m = [tilt(y, x, Rs, Ks, ps); -1.0] # surface normal / gradient
        normalize!(m)
        refract!(k, m, n[i], n[i+1])
        u = k[2] / k[3]
        v = k[1] / k[3]
    end
    return x, y
end

function full_trace(surfaces::Matrix, system::SystemOrRayBasis, H, k_rays = k_rays;
                    K = zeros(size(surfaces, 1) + 1), p = fill(zero, length(K)))
    H = abs(H)
    H ≤ 1.0 || throw(DomainError(H, "Domain: |H| ≤ 1.0"))
    real_chief = trace_chief_ray(surfaces, system)
    real_marginal = trace_marginal_ray(surfaces, system)
    EP_t = real_chief.z[1]
    Ū = real_chief.u[1]
    U = H * Ū
    u = tan(U)
    y_EP = real_marginal.y[1]
    h′ = transfer(system, [0.0, u], -EP_t, system.EBFD)[1]
    surfaces = [@view(surfaces[:,1:3]); [Inf 0.0 1.0]]
    surfaces[end-1,2] = system.EBFD * system.lens.n[end]
    V = 0.0
    k_rays_2 = div(k_rays, 2) + 1
    εy = Matrix{Float64}(undef, k_rays_2, k_rays)
    εx = Matrix{Float64}(undef, k_rays_2, k_rays)
    ρ = y_EP .* sqrt.(range(0.0, 1.0, k_rays)) # uniform sampling
    θ = range(0.0, π, k_rays_2) # take advantage of symmetry
    i = 1
    for ρᵢ ∈ ρ, θᵢ ∈ θ
        y = ρᵢ * cos(θᵢ)
        x = ρᵢ * sin(θᵢ)
        if typeof(system) <: RayBasis
            z0 = system.marginal.z[1]
            ȳ = system.chief.y[2] + system.chief.u[1] * z0
            EP_O = z0 - EP_t
            U = (ȳ - y) / EP_O
            V = -x / EP_O
        end
        y -= tan(U) * EP_t
        x, y = raytrace(surfaces, y, x, U, V, Vector{RealRay}; K, p)
        εy[i] = y - h′
        εx[i] = x
        i += 1
    end
    # take advantage of symmetry
    εy = [εy; @view(εy[end-1:-1:begin+1,:])]
    εx = [εx; -@view(εx[end-1:-1:begin+1,:])]
    nu = system.marginal.nu[end]
    # angle convention for Zernike polynomials
    θ_Z = range(start = π/2, step = step(θ), length = k_rays)
    return RealRayError(εx, εy, nu, ρ, θ_Z, H, σ(εx, εy))
end

function full_trace(surfaces::Layout{Aspheric}, system::System, H, k_rays = k_rays)
    full_trace(surfaces, system, H, k_rays; K = surfaces.K, p = surfaces.p)
end

function wavegrad(ε::RealRayError, λ = λ)
    map(field -> getfield(ε, field) * ε.nu / λ, (:x, :y))
end

function σ(εx, εy)
    n = length(εx)
    μx, μy = (sum(ε) / n for ε ∈ (εx, εy))
    sqrt((sum((εx .- μx) .^ 2) + sum((εy .- μy) .^ 2)) / n)
end
