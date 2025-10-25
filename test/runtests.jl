using Test
using OpticalRayTracing
using OpticalRayTracing: λ

# Verification is done using an example found in the book:
# Modern Optical Engineering, fourth edition, by Warren J. Smith
# Chapter Six, starting on page 119

# from Schott catalog
const air = 1.0
const SK4 = 1.61272
const SF2 = 1.64769
# n_F - n_C
const dn1 = 0.010450 # SK4
const dn2 = 0.019151 # SF2

# Cooke triplet prescription
const surfaces = [
    # R          t       n
    Inf          0.0     air; # object space
    +  37.40     5.90    SK4;
    - 341.48    12.93    air;
    -  42.65     2.50    SF2;
    +  36.40     2.00    air; 
    Inf          9.85    air; # stop
    + 204.52     5.90    SK4;
    -  37.05     0.0     air
]

# clear aperture semi-diameters
const a = [
    14.7,
    14.7,
    10.8,
    10.8,
    10.3,
    11.6,
    11.6
]

# Gaussian image height
const h′ = 21.248

# dispersion vector
const δn = [
    0.0,
    dn1,
    0.0,
    dn2,
    0.0,
    0.0,
    dn1,
    0.0
]

const third_order = [
    -0.709019 -0.070068 -0.006924 -0.330897 -0.033385;
    -0.75536 0.844458 -0.944066 -0.036241 1.095939;
    2.049816 -1.416428 0.978756 0.300215 -0.883772;
    0.493011 0.393733 0.314447 0.351763 0.532055;
    0.0 0.0 0.0 0.0 0.0;
    -0.035845 -0.087854 -0.215325 -0.06051 -0.676055;
    -1.229178 0.317877 -0.082206 -0.334022 0.107641
]

const system = solve(surfaces, a, h′)

const (; N) = system

# The third order aberration data in the book are Seidel S_I-S_V coefficients in the
# form of transverse ray errors in units of mm. Transverse ray aberrations are
# proportional to the negative partial derivative of the wavefront error times
# the scaling constant -1/n′ū′ ≈ 2 * f_number, however this book defines the
# coefficients using -1/2n′ū′ instead so the data is scaled appropriately to compare
# the WIJK coefficients in terms of waves.

# Their results appear to be calculated in 32-bit floating point precision.

const aberr = aberrations(surfaces, system, λ, δn)

@testset "coefficients" begin
    for i in axes(third_order, 1)
        @test aberr.W040[i] ≈ (-@view(third_order[:,1])[i] / (8 * λ * N)) atol = 1.0
    end
end
