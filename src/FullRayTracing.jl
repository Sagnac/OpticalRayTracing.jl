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
        # @. k += 2.0 * γ * m
        return NaN # TIR
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
    vec_length = size(surfaces, 1)-1
    xv = Vector{Float64}(undef, vec_length)
    yv = Vector{Float64}(undef, vec_length)
    for i = 1:vec_length
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
        xv[i] = x
        yv[i] = y
    end
    return xv, yv
end

function full_trace(surfaces::Matrix, system::System, H, k_rays = 10 * k_rays;
                    K = zeros(size(surfaces, 1) + 1), p = fill(zero, length(K)))
    H = abs(H)
    H ≤ 1.0 || throw(DomainError(H, "Domain: |H| ≤ 1.0"))
    stop = system.stop
    a_stop = abs(system.a[stop])
    real_chief = trace_chief_ray(surfaces, system)
    real_marginal = trace_marginal_ray(surfaces, system)
    EP_t = real_chief.z[1]
    Ū = real_chief.u[1]
    U = H * Ū
    u = tan(U)
    y_EP = abs(real_marginal.y[1])
    y0 = -u * EP_t
    h′ = transfer(system, [0.0, u], -EP_t, system.EBFD)[1]
    # extend the surface matrix to the paraxial image plane
    surfaces = [@view(surfaces[:,1:3]); [Inf 0.0 1.0]]
    surfaces[end-1,2] = system.EBFD * system.lens.n[end]
    V = 0.0
    k_rays_2 = div(k_rays, 2)
    εy = Float64[]
    εx = Float64[]
    r = Float64[]
    θ = Float64[]
    y = range(y0 - y_EP, y0 + y_EP, k_rays)
    # x = range(-y_EP, y_EP, k_rays)
    x = range(0.0, y_EP, k_rays_2)
    for yᵢ ∈ y, xᵢ ∈ x
        # if typeof(system) <: RayBasis
            # z0 = system.marginal.z[1]
            # ȳ = system.chief.y[2] + system.chief.u[1] * z0
            # EP_O = z0 - EP_t
            # U = (ȳ - y) / EP_O
            # V = -x / EP_O
        # end
        xv, yv = raytrace(surfaces, yᵢ, xᵢ, U, V, Vector{RealRay}; K, p)
        xf = xv[end]
        yf = yv[end]
        rᵢ = hypot(xv[stop], yv[stop])
        (rᵢ > a_stop || isnan(xf) || isnan(yf)) && continue
        θᵢ = atan(yv[stop], xv[stop])
        push!(εy, yf - h′)
        push!(εx, xf)
        push!(r, rᵢ)
        push!(θ, θᵢ)
    end
    # take advantage of symmetry
    εy = [εy; εy]
    εx = [εx; -εx]
    ρ = r / maximum(r)
    ρ = [ρ; ρ]
    θ = [θ; π .- θ]
    nu = system.marginal.nu[end]
    return RealRayError(εx, εy, nu, ρ, θ, H, σ(εx, εy))
end

function full_trace(surfaces::Layout{Aspheric}, system::System, H, k_rays = k_rays)
    full_trace(surfaces.M, system, H, k_rays;
               K = [surfaces.K; 0.0], p = [surfaces.p; zero])
end

function wavegrad(ε::RealRayError, λ = λ)
    map(field -> getfield(ε, field) * ε.nu / λ, (:x, :y))
end

function σ(εx, εy)
    n = length(εx)
    μx, μy = (sum(ε) / n for ε ∈ (εx, εy))
    sqrt((sum((εx .- μx) .^ 2) + sum((εy .- μy) .^ 2)) / n)
end
