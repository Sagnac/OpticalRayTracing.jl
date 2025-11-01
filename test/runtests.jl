using Test, Suppressor
using OpticalRayTracing
using OpticalRayTracing: λ, _vignetting, surface_ray

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

const (; f, N, marginal, chief, EFFD, EBFD, EP, XP) = system

const EFL = 101.181
const BFL = 77.405
const NA = 0.1443
const HFOV = 11.86
const PTZ_F = -2.935
const VL = 39.08 # system vertex length

const system_scale = 1e-3

@testset "system properties" begin
    @test f ≈ EFL atol = system_scale
    @test EBFD ≈ BFL atol = system_scale
    @test system.N ≈ 1/2NA atol = system_scale
    @test system.FOV ≈ 2HFOV atol = system_scale
    @test sum(@view(surfaces[2:end,2])) ≈ VL atol = system_scale
    @test system.stop == 5
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

const PI = surface_ray(@view(yui[:,3]))
const PIC = surface_ray(@view(ȳūī[:,3]))

const A, Ā = eachcol(@view(incidences(surfaces, system)[:,3:4]))

const h = 10.0
const s = -f + EFFD
const rays_2f = raytrace(system, h, s)
const rand_h, rand_s = (-1000 * rand() for i = 1:2)
const rand_rays = raytrace(system, rand_h, rand_s)

const trace_scale = 1e-2

@testset "raytrace validation" begin
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

@testset "general system ray tracing" begin
    s′ = f + EBFD
    XP_I = s′ - XP.t
    (; y, u, nu) = rays_2f.marginal
    ȳ, ū, nū = rays_2f.chief.y, rays_2f.chief.u, rays_2f.chief.nu
    @test iszero(y[end])
    @test u[end] ≈ -XP.D / 2XP_I
    @test u[end] ≈ -u[1]
    @test u[end] ≈ rays_2f.H / h
    @test ȳ[end] ≈ -h
    @test ū[end] ≈ -h / XP_I
    H_v = @. nū * y - nu * ȳ
    @test all(≈(rays_2f.H), H_v)
    H_v_rand = @. rand_rays[2].nu * rand_rays[1].y - rand_rays[1].nu * rand_rays[2].y
    @test all(≈(rand_rays.H), H_v_rand)
    EP_O = (rand_s - EP.t)
    u_in = -system.marginal.y[1] / EP_O
    ū_in = rand_h / EP_O
    y_in = -u_in * rand_s
    ȳ_in = -ū_in * EP.t
    yb, ub = system.M * [y_in, u_in]
    ȳb, ūb = transfer(system, [rand_h, ū_in], -rand_s, 0.0)
    rt_marginal = raytrace(system.lens, y_in, u_in)
    rt_chief = raytrace(system.lens, ȳ_in, ū_in)
    @test rand_rays.marginal.y[end-1] ≈ yb
    @test rand_rays.marginal.u[end] ≈ ub
    @test rand_rays.chief.y[end-1] ≈ ȳb
    @test rand_rays.chief.u[end] ≈ ūb
    @test all(@views rand_rays.marginal.y[2:end-1] .≈ rt_marginal.y[2:end])
    @test all(@view(rand_rays.marginal.nu[1:end-1,:]) .≈ rt_marginal.nu)
    @test all(isapprox.(surface_ray(rand_rays.chief.y), # zero comparison at stop
                        @view(rt_chief.y[2:end]); atol = 1e-12))
    @test all(@view(rand_rays.chief.nu[1:end-1,:]) .≈ rt_chief.nu)
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

const W040 = -0.186575
const W131 = -0.018282
const W222 = 0.044681
const W220P = -0.109691
const W311 = 0.142422
const W020 = -0.022389
const W111 = 0.027401

# Absolute tolerance of a quarter wave which is close to the data minimum
const aberr_scale = 0.25

const α = -inv(λ * N)

const aberr = aberrations(surfaces, system, λ, δn)

const longitudinal_petzval = 2N * aberr.W220P * (4 / α) # scaled appropriately

# Petzval radius of curvature
const ρ =  h′ ^ 2 / 2longitudinal_petzval

@testset "aberration coefficients" begin
    for i in axes(third_order, 1)
        @test aberr.spherical[i] ≈ (α * SI[i] / 8) atol = aberr_scale
        @test aberr.coma[i] ≈ (α * SII[i] / 2) atol = aberr_scale
        @test aberr.astigmatism[i] ≈ (α * SIII[i] / 2) atol = aberr_scale
        @test aberr.petzval[i] ≈ (α * SIV[i] / 4) atol = aberr_scale
        @test aberr.distortion[i] ≈ (α * SV[i] / 2) atol = aberr_scale
        @test aberr.axial[i] ≈ (α * PAC[i] / 4) atol = aberr_scale
        @test aberr.lateral[i] ≈ (α * PLC[i] / 2) atol = aberr_scale
    end
    @test aberr.W040 ≈ (α * W040 / 8) atol = aberr_scale
    @test aberr.W131 ≈ (α * W131 / 2) atol = aberr_scale
    @test aberr.W222 ≈ (α * W222 / 2) atol = aberr_scale
    @test aberr.W220P ≈ (α * W220P / 4) atol = aberr_scale
    @test aberr.W311 ≈ (α * W311 / 2) atol = aberr_scale
    @test aberr.W020 ≈ (α * W020 / 4) atol = aberr_scale
    @test aberr.W111 ≈ (α * W111 / 2) atol = aberr_scale
    @test ρ / f ≈ PTZ_F atol = system_scale
end

@testset "transfer matrix" begin
    (; f, EBFD, EFFD, P1, P2) = flatten(system.M)
    @test f ≈ system.f
    @test EBFD ≈ system.EBFD
    @test EFFD ≈ system.EFFD
    @test P1 ≈ system.P1
    @test P2 ≈ system.P2
    @test -/(reverse_transfer(system.M, [1.0, 0.0], 0.0, 0.0)...) ≈ EFFD
end

@testset "vignetting" begin
    partial = @suppress _vignetting(system, a)[4]
    @test partial == [1, 2, 3, 6, 7]
    slopes, FOVs = vignetting(system, a, FOV)
    systems = [solve(surfaces, a, h′) for h′ in FOVs]
    vignetting_matrices = @suppress [vignetting(system, a) for system in systems]
    a_ = (a for (i, a) in enumerate(a) if i != system.stop)
    for i = 1:3
        M = vignetting_matrices[i]
        @test any(a_ .≈ view(M, [1:system.stop-1; system.stop+1:length(a)], i+2))
    end
    ū = slopes[2]
    y = -ū * EP.t
    rt_half = raytrace(system.lens, y, ū, a; clip = true).ynu
    rt_clip = raytrace(system.lens, y - 1e-12, ū, a; clip = true).ynu
    @test !any(isnan, rt_half)
    @test any(isnan, rt_clip)
end
