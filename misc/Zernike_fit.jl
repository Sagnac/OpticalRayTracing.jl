# This is an example workflow for fitting spherical aberration
# transverse ray errors from a real ray trace to Zernike polynomials

using OpticalRayTracing, Zernike
using OpticalRayTracing: λ

include("../docs/setup.jl")

y, ε = TSA(surfaces, system)

ρ = y / maximum(y)
θ = range(0.0, 2π, length(ρ))

# convert to the derivative of the wavefront error in waves
∂W = ε * system.marginal.nu[end] / λ

OPDx = @. ∂W' * cos(θ)
OPDy = @. ∂W' * sin(θ)

degree = 9
orders_x = [(1, n) for n ∈ (1:2:degree)]
orders_y = [(-1, n) for n ∈ (1:2:degree)]

# can also just use degree instead of orders, but this should be slightly more
# accurate as it only explicitly includes the spherical aberration terms
∂x = W(ρ, θ, OPDx, orders_x)
∂y = W(ρ, θ, OPDy, orders_y)

ΔW = W(∂x, ∂y)

# you can discard any insignificant terms if you'd like
ΔW = Wavefront(Zernike.sieve(ΔW.v, 1e-7))

# if desired you can extract the Seidel coefficient

# Zernike sequential indices
j = [get_j(0, n) for n = 4:2:10]
# Zernike polynomial coefficients
α = getindex.(Z.(j), 4)
# take into account the normalization factor and sum the coefficients for ρ⁴
W040 = sum(Zernike.N(j) .* α .* ΔW[j])

# this can be compared to the Seidel approximation and the fit from the real trace

aberr = aberrations(surfaces, system)
SA3 = SA(y, ε, degree)[1] * system.marginal.nu[end] / 4λ

# compare with half a wave of difference
# results could differ by a wave or so depending on which
# parameter is used for R / rₚ
@show isapprox(W040, aberr.W040; atol = 0.5)
@show isapprox(W040, SA3; atol = 0.5)
@show W040

nothing
