# Pupil Sampling

Real rays can be traced over the entirety of the pupil rather than just the meridional direction in order to compute the RMS spot size in reference to the centroid from a set of geometric ray errors at the image plane.

```@setup pupil_sampling

```

```@example pupil_sampling
# method signature: full_trace(system, H, k_rays, focus)
ε = full_trace(system, 0.0) # on-axis
```

The second argument specifies the normalized field parameter, the third the number of rays to trace, and the fourth the focus position which defaults to the paraxial image plane.

The real [`spot_diagram`](@ref) can then be plotted using Makie.

Because transverse ray errors are proportional to the gradient of the wavefront error the [`wavegrad`](@ref) function can then be used to scale the `full_trace` result and return the wavefront partial derivative vectors in waves. `ε` also holds the corresponding pupil polar coordinates in the `r` and `t` fields so that they can be used along with the gradient to fit the geometric data to Zernike polynomials.
