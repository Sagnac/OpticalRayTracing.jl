# Seidel Aberrations

The surface coefficients can be extracted by querying the relevant `WIJK` fields. In addition, chromatic aberration can be quantified for a given `δn` dispersion vector. The results are optical path differences in waves at the `d` spectral line for a system in `mm`; the wavelength parameter can also be tuned.

```@setup aberrations

```

```@example aberrations
aberrations(surfaces, system)
```

The return value holds the relevant metadata in an `Aberration` structure. The named fields (`spherical`, etc.) are the sum totals. `petzval`, `tangential`, & `sagittal` refer to the respective field curvatures, while the `W020` and `W111` coefficient vectors refer to defocus & tilt from axial & lateral color, respectively.
