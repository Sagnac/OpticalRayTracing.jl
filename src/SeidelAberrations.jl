# Helium Fraunhofer d-line wavelength in millimeters
const λ = 587.5618e-6

Δ(x, y, i) = x[i+1] / y[i+1] - x[i] / y[i]

function aberrations(surfaces::AbstractMatrix, system::SystemOrRayBasis,
                     λ::Float64 = λ,
                     δn::Vector{Float64} = zeros(size(surfaces, 1)))
    (; marginal, chief, H) = system
    R = @view surfaces[2:end,1]
    n = marginal.n
    y = surface_ray(marginal.y)
    ȳ = surface_ray(chief.y)
    nu = marginal.nu
    nū = chief.nu
    u = marginal.u
    j = eachindex(R)
    A = [nu[i] + n[i] * y[i] / R[i] for i ∈ j]
    Ā = [(H + A[i] * ȳ[i]) / y[i] for i ∈ j]
    yΔ = [y[i] * Δ(u, n, i) for i ∈ j]
    yδ = [y[i] * Δ(δn, n, i) for i ∈ j]
    Δn2 = [inv(n[i+1]) ^ 2 - inv(n[i]) ^ 2 for i ∈ j]
    P = [(inv(n[i+1]) - inv(n[i])) / R[i] for i ∈ j]
    # per surface contributions
    spherical = @. -A ^ 2 * yΔ / 8λ
    coma = @. -A * Ā * yΔ / 2λ
    astigmatism = @. -Ā ^ 2 * yΔ / 2λ
    petzval = @. -H ^ 2 * P / 4λ
    sagittal = @. petzval + astigmatism / 2
    distortion = @. -Ā * (Ā ^ 2 * y * Δn2 - (H + Ā * y) * ȳ * P) / 2λ
    axial = @. A * yδ / 2λ
    lateral = @. Ā * yδ / λ
    medial = @. petzval + astigmatism
    tangential = @. petzval + 1.5 * astigmatism
    # system coefficients
    W040 = sum(spherical)
    W131 = sum(coma)
    W222 = sum(astigmatism)
    W311 = sum(distortion)
    # field curvatures
    W220P = sum(petzval)
    W220 = W220P + 0.5 * W222
    W220M = W220P + W222
    W220T = W220P + 1.5 * W222
    # chromatic aberrations
    # defocus from longitudinal
    W020 = sum(axial)
    # tilt from transverse
    W111 = sum(lateral)
    Aberration(W040, W131, W222, W220, W311, W020, W111, W220P, W220M, W220T,
               spherical, coma, astigmatism, sagittal, distortion, axial, lateral,
               petzval, medial, tangential, λ, sign(chief.y[end]))
end

function aberrations(system::System{Layout},
                     λ::Float64 = λ,
                     δn::Vector{Float64} = zeros(size(system.layout, 1)))
    aberrations(system.layout, system, λ, δn)
end

function (W::Aberration)(ρ, θ, H)
    H = abs(H)
    H ≤ 1.0 || throw(DomainError(H, "Domain: |H| ≤ 1.0"))
    0.0 ≤ ρ ≤ 1.0 || throw(DomainError(ρ, "Domain: 0.0 ≤ ρ ≤ 1.0"))
    H *= W.field_sign
    (; W040, W131, W222, W220, W311, W020, W111) = W
    w1 = W040 * ρ ^ 4
    w2 = W131 * H * ρ ^ 3 * cos(θ)
    w3 = W222 * H ^ 2 * ρ ^ 2 * cos(θ) ^ 2
    w4 = W220 * H ^ 2 * ρ ^ 2
    w5 = W311 * H ^ 3 * ρ * cos(θ)
    w6 = W020 * ρ ^ 2
    w7 = W111 * H * ρ * cos(θ)
    w = w1 + w2 + w3 + w4 + w5 + w6 + w7
    return w
end

function ray_error(ε, x, y, H)
    hypot(x, y) ≤ 1.0 || throw(DomainError((x, y), "Domain: hypot(x, y) ≤ 1.0"))
    H = abs(H)
    H ≤ 1.0 || throw(DomainError(H, "Domain: |H| ≤ 1.0"))
    H *= ε.field_sign
    (; W040, W131, W222, W220, W311, W020, W111, λ) = ε.W
    (; nu) = ε
    ε1_y = 4 * W040 * (x ^ 2 * y + y ^ 3)
    ε1_x = 4 * W040 * (y ^ 2 * x + x ^ 3)
    ε2_y = W131 * H * (x ^ 2 + 3 * y ^ 2)
    ε2_x = W131 * H * (2 * x * y)
    ε3_y = 2 * W222 * H ^ 2 * y
    ε3_x = 0
    ε4_y = 2 * W220 * H ^ 2 * y
    ε4_x = 2 * W220 * H ^ 2 * x
    ε5_y = W311 * H ^ 3
    ε5_x = 0
    ε6_y = 2 * W020 * y
    ε6_x = 2 * W020 * x
    ε7_y = W111 * H
    ε7_x = 0
    ε_y = (ε1_y + ε2_y + ε3_y + ε4_y + ε5_y + ε6_y + ε7_y) * λ / nu
    ε_x = (ε1_x + ε2_x + ε3_x + ε4_x + ε5_x + ε6_x + ε7_x) * λ / nu
    return ε_x, ε_y
end

(ε_y::RayError{Tangential})(y, H) = ray_error(ε_y, 0, y, H)[2]

(ε_x::RayError{Sagittal})(x, H) = ray_error(ε_x, x, 0, H)[1]

(ε::RayError{Skew})(x, y, H) = ray_error(ε, x, y, H)

(ε::RayError{Skew})((x, y)::NTuple{2, Float64}, H) = ε(x, y, H)

function RayError{T}(W::Aberration, s::SystemOrRayBasis) where T <: AbstractRay
    RayError{T}(W, s.marginal.nu[end], W.field_sign)
end

function TSA(surfaces::AbstractMatrix, system, k_rays::Int = k_rays)
    paraxial_marginal = system.marginal
    real_marginal = trace_marginal_ray(surfaces, system)
    real_chief = trace_chief_ray(surfaces, system)
    XP_t = real_chief.z[end] - real_chief.z[end-1]
    y_EP = range(real_marginal.y[1] / k_rays, real_marginal.y[1], k_rays)
    y_XP = Vector{Float64}(undef, k_rays)
    ε = Vector{Float64}(undef, k_rays)
    paraxial_BFD = paraxial_marginal.z[end] - paraxial_marginal.z[end-1]
    t = surface_to_focus(paraxial_BFD, real_marginal, paraxial_marginal)
    y_XP[end] = transfer(real_marginal.y[end-1], real_marginal.u[end], XP_t, RealRay)
    ε[end] = transfer(real_marginal, t)
    for (i, y) in enumerate(Base.Iterators.take(y_EP, k_rays - 1))
        ray = raytrace(surfaces, y, 0.0, RealRay)
        t = surface_to_focus(paraxial_BFD, ray)
        y_XP[i] = transfer(ray.y[end], ray.u[end], XP_t, RealRay)
        ε[i] = transfer(ray, t)
    end
    return y_XP, ε
end

TSA(system::System, k_rays::Int = k_rays) = TSA(system.layout, system, k_rays)

function SA(y::Vector, ε::Vector, degree::Int)
    if iseven(degree) || degree < 3
        throw(DomainError(degree, "Required: isodd(degree) && degree ≥ 3"))
    end
    yₚ = y / maximum(y)
    A = [yₚ ^ k for yₚ ∈ yₚ, k ∈ (3:2:degree)]
    return A \ ε
end
