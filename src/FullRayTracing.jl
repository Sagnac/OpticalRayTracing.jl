function refract!(k::Vector, m::Vector, n1::Float64, n2::Float64)
    η = n1 / n2
    γ = -k ⋅ m
    Δ = 1 - η ^ 2 * (1 - γ ^ 2)
    if Δ ≥ 0.0
        @. k = abs(η) * k + (η * γ - sign(η) * sqrt(Δ)) * m
    else
        @. k += 2.0 * γ * m # TIR
    end
    return k
end

function raytrace(surfaces::AbstractMatrix, y, x, U, V, ::Type{Vector{RealRay}};
                  K = zeros(size(surfaces, 1)), p = fill(zero, length(K)))
    R, t, n = eachcol(surfaces)
    tsy = copy(t)
    tsx = copy(t)
    u = tan(U)
    v = tan(V)
    k = [v, u, 1.0]
    normalize!(k) # direction cosines
    m = [0.0, 0.0, -1.0] # surface normal
    for i = 1:size(surfaces, 1)-1
        y += u * tsy[i]
        x += v * tsx[i]
        Rs = R[i+1]
        Ks = K[i+1]
        ps = p[i+1]
        sy = sag(y, U, Rs, Ks, ps)
        sx = sag(x, V, Rs, Ks, ps)
        y += sy * u
        x += sx * v
        tsy[i] += sy
        tsx[i] += sx
        tsy[i+1] -= sy
        tsx[i+1] -= sx
        ds_dy = tilt(y, Rs, Ks, ps)
        ds_dx = tilt(x, Rs, Ks, ps)
        m .= ds_dx, ds_dy, -1.0
        normalize!(m)
        refract!(k, m, n[i], n[i+1])
        u = k[2] / k[3]
        v = k[1] / k[3]
        U = atan(u)
        V = atan(v)
    end
    return x, y
end

function full_trace(surfaces::Matrix, system::SystemOrRayBasis, H, k_rays = k_rays;
                    K = zeros(size(surfaces, 1) + 1), p = fill(zero, length(K)))
    H = abs(H)
    0.0 ≤ H ≤ 1.0 || throw(DomainError(H, "Domain: 0.0 ≤ |H| ≤ 1.0"))
    real_chief = trace_chief_ray(surfaces, system)
    real_marginal = trace_marginal_ray(surfaces, system)
    U = H * real_chief.u[1]
    y_EP = real_marginal.y[1]
    y_max = iszero(U) ? y_EP : real_chief.y[2] + y_EP
    surfaces = [@view(surfaces[:,1:3]); [Inf 0.0 1.0]]
    surfaces[end-1,2] = system.EBFD * system.lens.n[end]
    V = 0.0
    # k_rays_2 = div(k_rays, 2) + 1
    εy = Matrix{Float64}(undef, k_rays, k_rays)
    εx = Matrix{Float64}(undef, k_rays, k_rays)
    # εy = Matrix{Float64}(undef, k_rays_2, k_rays)
    # εx = Matrix{Float64}(undef, k_rays_2, k_rays)
    ρ = range(0.0, y_max, k_rays)
    θ = range(0.0, 2π, k_rays)
    # θ = range(0.0, π, k_rays_2)
    i = 1
    for ρᵢ ∈ ρ, θᵢ ∈ θ
        y = ρᵢ * cos(θᵢ)
        x = ρᵢ * sin(θᵢ)
        if typeof(system) <: RayBasis
            ȳ = system.chief.y[2] + system.chief.u[1] * system.marginal.z[1]
            EP_O = system.marginal.z[1] - system.EP.t
            U = (ȳ - y) / EP_O
            V = -x / EP_O
        end
        x, y = raytrace(surfaces, y, x, U, V, Vector{RealRay}; K, p)
        εy[i] = y - H * system.chief.y[end]
        εx[i] = x
        i += 1
    end
    # εy = [εy; @view(εy[end-1:-1:1,:])]
    # εx = [εx; -@view(εx[end-1:-1:1,:])]
    return RealRayError(εx, εy, system.marginal.nu[end])
end

function full_trace(surfaces::Layout{Aspheric}, system::System, H, k_rays = k_rays)
    full_trace(surfaces, system, H, k_rays; K = surfaces.K, p = surfaces.p)
end

function wavegrad(ε::RealRayError, λ = λ)
    map(field -> getfield(ε, field) * ε.nu / λ, (:x, :y))
end
