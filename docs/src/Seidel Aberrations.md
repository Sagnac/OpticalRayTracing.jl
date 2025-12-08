# Seidel Aberrations

The total third order system coefficients can be extracted by querying the relevant `WIJK` fields. In addition, chromatic aberration can be quantified for a given `δn` dispersion vector. The results are optical path differences in waves at the `d` spectral line for a system in `mm`; the wavelength parameter can also be tuned.

```@setup aberrations

```

```@example aberrations
# for an object "at infinity"
# -- collimated rays, maximum paraxial field of view, focusing at the focal point
W = aberrations(system)
```

The return value holds the relevant metadata in an `Aberration` structure which can be evaluated over the pupil if called with polar coordinates and a normalized field value. The long form named fields (`spherical`, `coma`, `astigmatism`, etc.) hold the per surface coefficient vectors. `W220P` / `petzval`, `W220T` / `tangential`, & `W220` / `sagittal` refer to the respective field curvatures while the `W020` & `W111` coefficients refer to defocus & tilt from axial & lateral color, respectively.

The coefficients assume a positive valued field.

Third order coefficients can also be computed for a given object height & axial position by passing a `RayBasis` instead of a `System`:

```@example aberrations
rays = raytrace(system, 10.0, -system.f + system.EFFD)
W_4f = aberrations(surfaces, rays)
```
