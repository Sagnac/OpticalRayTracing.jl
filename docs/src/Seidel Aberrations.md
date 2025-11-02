# Seidel Aberrations

The total system coefficients can be extracted by querying the relevant `WIJK` fields. In addition, chromatic aberration can be quantified for a given `δn` dispersion vector. The results are optical path differences in waves at the `d` spectral line for a system in `mm`; the wavelength parameter can also be tuned.

```@setup aberrations

```

```@example aberrations
W = aberrations(surfaces, system)
```

The return value holds the relevant metadata in an `Aberration` structure which can be evaluated over the pupil if called with polar coordinates and a normalized field value. The long form named fields (`spherical`, `coma`, `astigmatism`, etc.) hold the per surface coefficient vectors. `W220P` / `petzval`, `W220T` / `tangential`, & `W220` / `sagittal` refer to the respective field curvatures while the `W020` & `W111` coefficients refer to defocus & tilt from axial & lateral color, respectively.

## Plots

Aberration plots are supported with the installation of any Makie backend. The fan plots are interactive and include a slider which allow dynamic adjustment of the field.

```julia
using GLMakie

wavefan(W)
rayfan(W, system)
field_curves(W, system)
percent_distortion(W, system)
```

![wavefan](assets/images/wavefan.png)
![rayfan](assets/images/rayfan.png)
![field_curves](assets/images/field_curves.png)
![percent_distortion](assets/images/percent_distortion.png)
