using Test
using OpticalRayTracing
using OpticalRayTracing: λ

# Verification is done using an example found in the book:
# Modern Optical Engineering, fourth edition, by Warren J. Smith
# Chapter Six, starting on page 119

const air = 1.0
# from Schott catalog
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
const a = [14.7, 14.7, 10.8, 10.8, 10.3, 11.6, 11.6]

# Gaussian image height
const h′ = 21.248

# dispersion vector
const δn = [0.0, dn1, 0.0, dn2, 0.0, 0.0, dn1, 0.0]

const system = solve(surfaces, a, h′)

const EFL = 101.181
const BFL = 77.405
const NA = 0.1443
const HFOV = 11.86
const PTZ_F = -2.935
const VL = 39.08 # system vertex length

const system_scale = 1e-3

@testset "system properties" begin
    @test system.f ≈ EFL atol = system_scale
    @test system.EBFD ≈ BFL atol = system_scale
    @test system.N ≈ 1/2NA atol = system_scale
    @test system.FOV ≈ 2HFOV atol = system_scale
    @test sum(@view(surfaces[2:end,2])) ≈ VL atol = system_scale
end

const yui = [
    14.6 0.0 0.0;
    14.6 -0.148315 0.390374;
    13.724943 -0.263817 -0.188507;
    10.313791 -0.065055 -0.505641;
    10.151154 0.073436 0.213823;
    10.298026 0.073436 0.073436;
    11.021371 0.025062 0.127325;
    11.169234 -0.144296 -0.276402;
    0.0 -0.144296 -0.144296
]

const ȳūī = [
    0.0 0.21 0.21;
    -6.411174 0.195343 0.038578;
    -5.25865 0.324469 0.210743;
    -1.063264 0.187124 0.349399;
    -0.595454 0.297727 0.170765;
    -3.3307e-15 0.297727 0.297727;
    2.932611 0.179164 0.312066;
    3.989677 0.222961 0.07148;
    21.248022 0.222961 0.222961
]

const y, u = eachcol(yui)
const ȳ, ū = eachcol(ȳūī)

const PI = @view(@view(yui[:,3])[begin+1:end-1])
const PIC = @view(@view(ȳūī[:,3])[begin+1:end-1])

const A, Ā = eachcol(@view(incidences(surfaces, system)[:,3:4]))

const (; marginal, chief) = system

const trace_scale = 1e-2

@testset "raytrace" begin
    for i in eachindex(y)
        @test marginal.y[i] ≈ y[i] atol = trace_scale
        @test marginal.u[i] ≈ u[i] atol = trace_scale
        @test chief.y[i] ≈ ȳ[i] atol = trace_scale
        @test chief.u[i] ≈ ū[i] atol = trace_scale
    end
    for i in eachindex(A)
        @test A[i] ≈ PI[i] atol = trace_scale
        @test Ā[i] ≈ PIC[i] atol = trace_scale
    end
end

# The third order aberration data in the book are Seidel S_I-S_V coefficients in the
# form of transverse ray errors in units of mm. Transverse ray aberrations are
# proportional to the negative partial derivative of the wavefront error times
# the scaling constant -1/n′ū′ ≈ 2 * f_number, however this book defines the
# coefficients using -1/2n′ū′ instead, except for the chromatic contributions which
# use -1/n′ū′ so the data is scaled appropriately to compare the WIJK coefficients
# in terms of waves.

# Their results appear to be calculated in 32-bit floating point precision.

const third_order = [
    -0.709019 -0.070068 -0.006924 -0.330897 -0.033385;
    -0.75536 0.844458 -0.944066 -0.036241 1.095939;
    2.049816 -1.416428 0.978756 0.300215 -0.883772;
    0.493011 0.393733 0.314447 0.351763 0.532055;
    0.0 0.0 0.0 0.0 0.0;
    -0.035845 -0.087854 -0.215325 -0.06051 -0.676055;
    -1.229178 0.317877 -0.082206 -0.334022 0.107641
]

const SI = @view third_order[:,1]
const SII = @view third_order[:,2]
const SIII = @view third_order[:,3]
const SIV = @view third_order[:,4]
const SV = @view third_order[:,5]

# primary chromatic aberrations
const PAC = [-0.25596, -0.187385, 0.419729, 0.287842, 0.0, -0.063021, -0.223595]
const PLC = [-0.025295, 0.209488, -0.290034, 0.229879, 0.0, -0.154461, 0.057824]

const spherical = -0.186575
const coma = -0.018282
const astigmatism = 0.044681
const petzval = -0.109691
const distortion = 0.142422
const axial = -0.022389
const lateral = 0.027401

# Absolute tolerance of a quarter wave which is close to the data minimum
const aberr_scale = 0.25

const (; N) = system

const α = -inv(λ * N)

const aberr = aberrations(surfaces, system, λ, δn)

const longitudinal_petzval = 2N * aberr.petzval * (4 / α) # scaled appropriately

# Petzval radius of curvature
const ρ =  h′ ^ 2 / 2longitudinal_petzval

@testset "aberration coefficients" begin
    for i in axes(third_order, 1)
        @test aberr.W040[i] ≈ (α * SI[i] / 8) atol = aberr_scale
        @test aberr.W131[i] ≈ (α * SII[i] / 2) atol = aberr_scale
        @test aberr.W222[i] ≈ (α * SIII[i] / 2) atol = aberr_scale
        @test aberr.W220P[i] ≈ (α * SIV[i] / 4) atol = aberr_scale
        @test aberr.W311[i] ≈ (α * SV[i] / 2) atol = aberr_scale
        @test aberr.W020[i] ≈ (α * PAC[i] / 4) atol = aberr_scale
        @test aberr.W111[i] ≈ (α * PLC[i] / 2) atol = aberr_scale
    end
    @test aberr.spherical ≈ (α * spherical / 8) atol = aberr_scale
    @test aberr.coma ≈ (α * coma / 2) atol = aberr_scale
    @test aberr.astigmatism ≈ (α * astigmatism / 2) atol = aberr_scale
    @test aberr.petzval ≈ (α * petzval / 4) atol = aberr_scale
    @test aberr.distortion ≈ (α * distortion / 2) atol = aberr_scale
    @test aberr.axial ≈ (α * axial / 4) atol = aberr_scale
    @test aberr.lateral ≈ (α * lateral / 2) atol = aberr_scale
    @test ρ / system.f ≈ PTZ_F atol = system_scale
end
