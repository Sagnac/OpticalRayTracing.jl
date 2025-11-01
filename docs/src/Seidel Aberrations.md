# Seidel Aberrations

The total system coefficients can be extracted by querying the relevant `WIJK` fields. In addition, chromatic aberration can be quantified for a given `δn` dispersion vector. The results are optical path differences in waves at the `d` spectral line for a system in `mm`; the wavelength parameter can also be tuned.

```@setup aberrations

```

```@example aberrations
aberrations(surfaces, system)
```

The return value holds the relevant metadata in an `Aberration` structure. The long form named fields (`spherical`, `coma`, `astigmatism`, etc.) hold the per surface coefficient vectors. `W220P` / `petzval`, `W220T` / `tangential`, & `W220` / `sagittal` refer to the respective field curvatures while the `W020` & `W111` coefficients refer to defocus & tilt from axial & lateral color, respectively.
