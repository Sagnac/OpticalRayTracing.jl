"""
    Lens(surfaces::AbstractMatrix)

Construct a lens from a 3 column surface matrix where the first column corresponds to the radii of curvature, the second column to the thicknesses, and the third to the refractive indices.

The computed result wraps a two column lens matrix where the first column refers to the reduced thickness `τ = t /n` of the distance prior to the surface and the second column refers to the surface power.

```jldoctest
julia> Lens([
           Inf 0.0 1.0 # object space
           50.0 3.0 1.5
           -50.0 0.0 1.0
       ])
2×2 Lens:
 0.0  0.01
 2.0  0.01
```

----

    Lens(M, n)

Construct a lens directly from a two column lens matrix `M` containing the reduced thicknesses and surface powers corresponding to a refractive index vector `n`.

```jldoctest
julia> Lens([0.0 0.01; 2.0 0.01], [1.0, 1.5, 1.0])
2×2 Lens:
 0.0  0.01
 2.0  0.01
```

See also: [`compute_surfaces`](@ref).
"""
Lens

"""
    aberrations(surfaces, system, [λ], [δn])

Compute the third order wavefront error coefficients based on Seidel sums for a given system or ray basis. Results are in waves at the d-line by default. The dispersion vector `δn` allows calculation of chromatic aberration.

The returned `Aberration` type holds the system totals in the `WIJK` fields (e.g. `W040`) while the per surface contributions are held in the long form fields (e.g. `spherical`).

The coefficients assume a positive valued field.
"""
aberrations

"""
    caustic(surfaces, system, [k_rays::Int]; [theme], kwargs...)

Traces real rays and plots the caustic extended to either the paraxial or marginal focus. Requires using a `Makie` backend. `k_rays` controls how many rays are plotted.

See also: [`solve`](@ref), [`raytrace`](@ref), [`trace_marginal_ray`](@ref).
"""
caustic

"""
    compute_surfaces(lens::Lens)

Compute the surface matrix from a given `Lens`.

```jldoctest
julia> compute_surfaces(Lens([0.0 0.01; 2.0 0.01], [1.0, 1.5, 1.0]))
3×3 Matrix{Float64}:
  Inf   0.0  1.0
  50.0  3.0  1.5
 -50.0  0.0  1.0
```

See also: [`Lens`](@ref).
"""
compute_surfaces

"""
    field_curves(W::Aberration, s::SystemOrRayBasis; k = k, kwargs...)

Plot the third order longitudinal astigmatic field curves. Requires using a `Makie` backend. `k` controls the discretization / number of plot elements.

`H` on the plot refers to the normalized field parameter.

`P`, `T`, & `S` on the plot refer to the Petzval field, the Tangential field, and the Sagittal field.

See also: [`percent_distortion`](@ref), [`rayfan`](@ref), [`spot_size`](@ref), [`wavefan`](@ref).
"""
field_curves

"""
    flatten(M::TransferMatrix)

Return the cardinal points for the input vertex matrix.
"""
flatten

"""
    incidences(surfaces, system)

Compute the paraxial incidence angles over the system. Returns a four column matrix with the first two columns referring to the optical angles `ni` & `nī` for the marginal & chief ray, respectively, and the last two columns referring to the physical angles `i`, `ī`.
"""
incidences

"""
    percent_distortion(W::Aberration, s::SystemOrRayBasis; k = k, kwargs...)

Plot the percent distortion for the given system of aberrations. Requires using a `Makie` backend. `k` controls the discretization / number of plot elements.

`H` on the plot refers to the normalized field parameter.

See also: [`field_curves`](@ref), [`rayfan`](@ref), [`spot_size`](@ref), [`wavefan`](@ref).
"""
percent_distortion

"""
    rayfan(W::Aberration, s::SystemOrRayBasis; k = k, kwargs...)

Plot the third order ray intercept curves for the given system of aberrations. Requires using a `Makie` backend. `k` controls the discretization / number of plot elements.

`H` on the plot refers to the normalized field parameter and can be adjusted with a slider.

See also: [`field_curves`](@ref), [`percent_distortion`](@ref), [`spot_size`](@ref), [`wavefan`](@ref).
"""
rayfan

"""
    rayplot(system)
    rayplot(surfaces, a, system)

Plots the given system or ray bundle using `Makie`; requires installation of a `Makie` backend such as `GLMakie`.

See also: [`rayplot!`](@ref).
"""
rayplot

"""
    rayplot!(system)
    rayplot!(surfaces, a, system)

Mutating version of `rayplot` which draws on top of the active figure.

See also: [`rayplot`](@ref).
"""
rayplot!

"""
    raypoints(marginal::ParaxialRay{Marginal}, chief::ParaxialRay{Chief})

Return the coordinates for the given ray basis; useful for manual plotting.
"""
raypoints

"""
    raytrace(surfaces, y, ω, [a]; clip = false)
    raytrace(lens::Lens, y, ω)

Trace a paraxial ray with height `y` and angle / slope `ω = nu` returning a `ParaxialRay{Tangential}` holding the `ynu` matrix.

If a set of aperture semi-diameters `a` is specified and `clip = true` any blocked ray will result in `NaN` in the trace matrix.

----

    raytrace(system::System, ȳ, s)

Trace the paraxial marginal and chief ray through a given system for a given object height and distance from the objective; these are signed quantities.

See also: [`solve`](@ref), [`incidences`](@ref).

----

    raytrace(surfaces::AbstractMatrix, y, U, ::Type{RealRay})

Trace a real ray through the surfaces with initial height `y` and initial angle `U`.

See also: [`trace_marginal_ray`](@ref).
"""
raytrace

"""
    reverse_transfer(M::TransferMatrix, v::Vector, τ′, τ)

Transfer an input `[y, nu]` vector in the reverse direction, extending by image & object space distances `τ′`, `τ`.

See also: [`transfer`](@ref).
"""
reverse_transfer

"""
    SA(y::Vector, ε::Vector, degree::Int)

Fits spherical aberration transverse ray errors to a polynomial of the given `degree` using linear least squares. `y` & `ε` are the XP heights and ray errors returned by [`TSA`](@ref). The returned vector corresponds to the coefficients for a tangential ray aberration expansion of the form: `ε = By³ + Cy⁵ + Dy⁷ …`
"""
SA

"""
    TSA(surfaces, system, [k_rays::Int])

Computes the transverse ray aberration errors for spherical aberration using a real ray trace. `k_rays` controls how many rays are plotted. Returns vectors with the ray heights at the exit pupil and the ray errors at the paraxial plane.

See also: [`SA`](@ref), [`caustic`](@ref), [`raytrace`](@ref).
"""
TSA

"""
    scale!(M::Lens)

Convert the second column powers from dioptres to inverse millimeters.
"""
scale!

"""
    solve(surfaces::Matrix, a::Vector, h′::Float64)
    solve(lens::Lens, a, h′)

Solve a lens, returning a `System` holding the relevant properties. `a` denotes the clear aperture semi-diameters while `h′` indicates the image height.

Fields are:

* `f`: effective focal length;
* `EBFD`: effective back focal distance;
* `EFFD`: effective front focal distance;
* `N`: f-number;
* `FOV`: full field of view in degrees;
* `stop`: surface matrix index denoting the aperture stop;
* `EP`: entrance pupil containing fields: `D`: diameter, `t`: distance from the front vertex;
* `XP`: exit pupil containing fields: `D`: diameter, `t`: distance from the back vertex;
* `marginal`: marginal ray;
* `chief`: chief ray;
* `trace`: the full `yu` ray trace for the marginal & chief rays;
* `H`: Lagrange invariant;
* `P1`: front principal plane location w.r.t the front vertex;
* `P2`: back principal plane location w.r.t the back vertex;
* `PN`: nodal point shift from the principal points;
* `M`: the vertex matrix;
* `lens`: `Lens` of the system.

See also: [`raytrace`](@ref), [`incidences`](@ref).
"""
solve

"""
    spot_size(W::Aberration, s::SystemOrRayBasis; k = k, kwargs...)

Plot the spot diagram over the image plane using third order aberration data. Requires using a `Makie` backend. `k` controls the discretization / number of plot elements.

`H` on the plot refers to the normalized field parameter and can be adjusted with a slider.

See also: [`field_curves`](@ref), [`percent_distortion`](@ref), [`rayfan`](@ref), [`spot_size`](@ref).
"""
spot_size

"""
    trace_chief_ray(lens::Lens, stop::Int, marginal::ParaxialRay{Marginal}, h′)

Find the paraxial chief ray for the given lens.

See also: [`solve`](@ref).

----

    trace_chief_ray(surfaces, system::System; atol = sqrt(eps()))

Find the real chief ray for the given surfaces and system.

See also: [`trace_marginal_ray`](@ref), [`raytrace`](@ref).
"""
trace_chief_ray

"""
    trace_marginal_ray(lens::Lens, a, [ω = 0.0])

Find the paraxial marginal ray for the given lens and aperture sizes.

See also: [`solve`](@ref).

----

    trace_marginal_ray(surfaces, system::System; atol = sqrt(eps()))

Find the real marginal ray for the given surfaces and system.

See also: [`trace_chief_ray`](@ref), [`raytrace`](@ref).
"""
trace_marginal_ray

"""
    transfer(M::TransferMatrix, v::Vector, τ, τ′)

Transfer an input `[y, nu]` vector, extending by object & image space distances `τ`, `τ′`.

See also: [`reverse_transfer`](@ref).
"""
transfer

"""
    vignetting(system, [a = system.a])

Perform a vignetting analysis on the `System` for the given aperture sizes.

Determines the maximum FOVs corresponding to those sizes as well as returning a `Vignetting` data structure holding a semi-diameter matrix in its `M` field with columns corresponding to conditions:
* input;
* stop limited;
* unvignetted;
* half-vignetted;
* fully vignetted.

The fields of view are stored in the `FOV` field with rows corresponding to the unvignetted, half-vignetted, and fully vignetted cases and the columns corresponding to the full fields of view in degrees, paraxial chief ray slopes, and image heights, respectively.

----

    vignetting(rays::RayBasis, system::System, a::AbstractVector = system.a)

Same as above, but for a traced `RayBasis` corresponding to the marginal & chief rays of a non-collimated beam / an object not at infinity.
"""
vignetting

"""
    wavefan(W::Aberration; k = k, kwargs...)

Plot the input aberrations as a a composite wavefront error over the pupil in the tangential and sagittal directions. Requires using a `Makie` backend. `k` controls the discretization / number of plot elements.

`H` on the plot refers to the normalized field parameter and can be adjusted with a slider.

See also: [`field_curves`](@ref), [`percent_distortion`](@ref), [`rayfan`](@ref), [`spot_size`](@ref).
"""
wavefan
