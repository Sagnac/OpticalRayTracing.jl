using Test
using OpticalRayTracing
using OpticalRayTracing: λ, surface_ray, surface_to_focus

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

# Cooke triplet layout
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

const (; f, N, marginal, chief, EFFD, EBFD, EP, XP, stop) = system

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
    @test u[end] ≈ -XP.D / 2XP_I ≈ -u[1] ≈ rays_2f.H / h
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
# proportional to the negative partial derivative of the wavefront error with respect
# to the normalized pupil coordinate times the ratio of the radius of the reference
# sphere to the radius of the exit pupil which is equal to -1/n′u′; however this book
# defines the wave coefficients using a different sign convention and -1/2n′u′ for
# the scaling constant instead, except for the first order chromatic contributions
# which use -1/n′u′; this is because their coefficients correspond to a transverse
# ray error expansion (B_1-B_5 in the book) which absorb all of the derivative and
# scaling prefactors so the data must be scaled appropriately to compare the WIJK
# coefficients in terms of waves.

const α = 2u[end] / λ

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

const SA5 = 0.125426
const SA7 = 0.048670

# Absolute tolerance of a quarter wave which is close to the data minimum
const wave_scale = 0.25

const aberr = aberrations(surfaces, system, λ, δn)

const longitudinal_petzval = -8N ^ 2 * aberr.W220P * λ

# Petzval radius of curvature
const ρ =  h′ ^ 2 / 2longitudinal_petzval

@testset "aberration coefficients" begin
    for i in axes(third_order, 1)
        @test aberr.spherical[i] ≈ (α * SI[i] / 8) atol = wave_scale
        @test aberr.coma[i] ≈ (α * SII[i] / 2) atol = wave_scale
        @test aberr.astigmatism[i] ≈ (α * SIII[i] / 2) atol = wave_scale
        @test aberr.petzval[i] ≈ (α * SIV[i] / 4) atol = wave_scale
        @test aberr.distortion[i] ≈ (α * SV[i] / 2) atol = wave_scale
        @test aberr.axial[i] ≈ (α * PAC[i] / 4) atol = wave_scale
        @test aberr.lateral[i] ≈ (α * PLC[i] / 2) atol = wave_scale
    end
    @test aberr.W040 ≈ (α * W040 / 8) atol = wave_scale
    @test aberr.W131 ≈ (α * W131 / 2) atol = wave_scale
    @test aberr.W222 ≈ (α * W222 / 2) atol = wave_scale
    @test aberr.W220P ≈ (α * W220P / 4) atol = wave_scale
    @test aberr.W311 ≈ (α * W311 / 2) atol = wave_scale
    @test aberr.W020 ≈ (α * W020 / 4) atol = wave_scale
    @test aberr.W111 ≈ (α * W111 / 2) atol = wave_scale
    @test ρ / f ≈ PTZ_F atol = system_scale
    # Petzval curvature
    Φᵢ = @view system.lens[:,2]
    n′ = @view system.lens.n[begin+1:end]
    n = @view system.lens.n[1:end-1]
    PTZC = -sum(Φᵢ ./ (n′ .* n))
    @test inv(ρ) ≈ PTZC
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
    vig = vignetting(system, a)
    @test vig.partial == [1, 2, 3, 6, 7]
    _, slopes, heights = eachcol(vig.FOV)
    systems = [solve(surfaces, a, h′) for h′ in heights]
    vignetting_matrices = [vignetting(system).M for system in systems]
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

@testset "real raytracing" begin
    yu_par = raytrace(surfaces, 1.0, 0.0).yu
    yu_real = raytrace(surfaces, 1.0, 0.0, RealRay).yu
    R, t = eachcol(surfaces)
    k = length(R) - 1
    ϵ_θ = inv(minimum(abs.(R))) ^ 3 / 6 * k * (k + 1) / 2
    ϵ_y = ϵ_θ * maximum(t)
    δ_θ = sum(@views(abs.(yu_par[2:end,2] .- yu_real[2:end,2])))
    δ_y = sum(@views(abs.(yu_par[2:end,1] .- yu_real[2:end,1])))
    @test δ_θ < ϵ_θ
    @test δ_y < ϵ_y
    atol = sqrt(eps())
    real_marginal = trace_marginal_ray(surfaces, system; atol)
    @test real_marginal.y[begin+stop] ≈ a[stop] atol = atol
    rt = raytrace(surfaces, real_marginal.y[1], real_marginal.u[1], RealRay)
    @test rt.yu[2:end,:] ≈ @view(real_marginal.yu[2:end-1,:]) atol = atol
    # fitting order and relative error are arbitrary at the moment
    B1 = SA(TSA(surfaces, system)..., 9)[1]
    @test abs(B1 / W040 - 1) < 0.05
    real_chief = trace_chief_ray(surfaces, system; atol)
    @test abs(real_chief.y[begin+stop]) < atol
    t = surface_to_focus(EBFD, real_chief, marginal)
    @test transfer(real_chief, t) ≈ chief.y[end]
    y_vertex = real_chief.y[2] - tan(real_chief.u[1]) * real_chief.z[2]
    rt = raytrace(surfaces, y_vertex, real_chief.u[1], RealRay)
    @test rt.yu[2:end,:] ≈ @view(real_chief.yu[2:end-1,:]) atol = atol
end

const ray_scale = 1e-3

@testset "transverse ray errors" begin
    ΔW = aberrations(surfaces, system)
    ε_y = RayError{Tangential}(ΔW, system)
    ε_x = RayError{Sagittal}(ΔW, system)
    @test ε_y(1.0, 1.0) ≈ +(W040, 3W131, 3W222, W220P, W311) atol = ray_scale
    @test ε_x(1.0, 1.0) ≈ +(W040, W222, W220P) atol = ray_scale
    ε = RayError{Skew}(ΔW, system)
    ρ = rand()
    θ = 2π * rand()
    H = rand()
    x = ρ * sin(θ)
    y = ρ * cos(θ)
    εx, εy = ε(x, y, H)
    @test εy ≈ W040 * ρ ^ 3 * cos(θ) +
               W131 * ρ ^ 2 * H * (2 + cos(2θ)) +
               (3W222 + W220P) * ρ * H ^ 2 * cos(θ) +
               W311 * H ^3 atol = ray_scale
    @test εx ≈ W040 * ρ ^ 3 * sin(θ) +
               W131 * ρ ^ 2 * H * sin(2θ) +
               (W222 + W220P) * ρ * H ^ 2 * sin(θ) atol = ray_scale
end

const optim_scale = 0.05

@testset "optimization" begin
    n = rand() * 0.2 + 1.5
    R1 = 100.0 * (n - 1)
    R2 = -R1
    t = 5.0
    a = fill(15.0, 2)
    surfaces = [Inf 0.0 1.0; R1 t n; R2 0.0 1.0]
    system = solve(surfaces, a, 10.0)
    v = 2:3
    aberr = [:W040]
    constraint = Dict(:f => system.f)
    layout, new_system = optimize(surfaces, system, v, constraint, aberr)
    # optimal lens bending shape factor for conjugate parameter M = 1
    X_min = 2 * (n ^ 2 - 1) / (n + 2)
    R1, R2 = layout[2:3]
    X = (R2 + R1) / (R2 - R1)
    @test X ≈ X_min rtol = optim_scale
    @test new_system.f ≈ system.f rtol = optim_scale
end

@testset "aspherics" begin
    aspheric_surface = Layout{Aspheric}([ # parabolic reflector
        Inf    0.0  1.0  0.0
        -100.0 0.0 -1.0 -1.0
    ])
    a = 30.0
    aspheric_system = solve(aspheric_surface, fill(a, 2), 21.0)
    real_marginal = trace_marginal_ray(aspheric_surface, aspheric_system)
    # should exhibit zero spherical aberration
    @test aspheric_system.marginal.z[end] == real_marginal.z[end] == -50.0
end
