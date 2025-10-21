# Helium Fraunhofer d-line wavelength in millimeters
const λ = 587.5618e-6

Δ(x, y, i) = x[i+1] / y[i+1] - x[i] / y[i]

function aberrations(surfaces::Matrix{Float64}, system::System,
                     λ::Float64 = λ,
                     δn::Vector{Float64} = zeros(size(surfaces, 1)))
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
