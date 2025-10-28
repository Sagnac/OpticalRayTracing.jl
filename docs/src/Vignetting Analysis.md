# Vignetting Analysis

Given a system or ray bundle and a set of apertures this package provides a basic tool for evaluating vignetting criteria at each aperture. Along with identifying, by index, which apertures cause partial or full vignetting, the returned 5 column semi-diameter matrix contains:
* The provided aperture size;
* The marginal ray limit: the size at which the aperture becomes the stop due to limiting the on-axis ray bundle;
* The size necessary for no vignetting to occur;
* The size at which the half-vignetting condition is satisified;
* The size at which the aperture clips the entire off-axis ray bundle; this is the minimum aperture size at which the given field of view is still supported down to approximately zero irradiance.

In addition, the maximum system field of view for the unvignetted, half-vignetted, and fully-vignetted cases is determined for the input aperture sizes.

By passing the exported `FOV` type the second method allows directly extracting the maximum supported paraxial principal ray slopes and the corresponding image heights for all three vignetting conditions.

In essence, the matrix represents the supported aperture sizes for each case given the system field of view and the returned vectors from the second method correspond to the maximal fields for the given aperture sizes.

Any `NaN`s in the matrix signify that there is no possible size which fits the criteria due to the aforementioned limiting factor.

```@setup vignetting

```

```@repl vignetting
vignetting(system, a)
```

```@repl vignetting
vignetting(system, a, FOV)
```
