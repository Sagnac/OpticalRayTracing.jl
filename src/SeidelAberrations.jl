# Helium Fraunhofer d-line wavelength in millimeters
const λ = 587.5618e-6

Δ(x, y, i) = x[i+1] / y[i+1] - x[i] / y[i]

function aberrations(surfaces::Matrix{Float64}, system::SystemOrRayBasis,
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
    P = [(inv(n[i+1]) - inv(n[i])) / R[i] for i ∈ j]
    # per surface contributions
    spherical = @. -A ^ 2 * yΔ / 8λ
    coma = @. -A * Ā * yΔ / 2λ
    astigmatism = @. -Ā ^ 2 * yΔ / 2λ
    petzval = @. -H ^ 2 * P / 4λ
    sagittal = @. petzval + astigmatism / 2
    distortion = @. Ā / A * (astigmatism + 2 * petzval)
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
               petzval, medial, tangential)
end

function (W::Aberration)(ρ, θ, H = 1)
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
