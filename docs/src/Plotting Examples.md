# Plotting Examples

Plotting is supported with the installation of any [Makie](https://github.com/MakieOrg/Makie.jl) backend.

## Ray Trace Plots

Supported keyword arguments are:
* `theme` (using Makie's theme customizability);
* `ray_colors` (a two-tuple of the marginal & chief ray colors);
* `surface_color` (for the optical element surfaces & image plane);
* anything Makie's `lines` accepts.

Using the mutating `rayplot!` version will draw on top of the current figure.

```julia
using GLMakie

surface_color = :white
ray_colors = (:cyan, :red)

rayplot(system; theme = theme_black(), surface_color, ray_colors)

rays = raytrace(system, -24.0, -1.5 * system.f)
fig = rayplot(rays)
rays = raytrace(system, 24.0, -1.5 * system.f)
rayplot!(rays)

surface = [Inf 0.0 1.0; -100.0 0.0 -1.0]
system = solve(surface, fill(30.0, 2), 21.0)
# traces real rays
fig = caustic(system; theme = theme_black(), surface_color)
```

![rayplot](assets/images/rayplot.png)
![rays](assets/images/rays.png)
![caustic](assets/images/caustic.png)

## Aberration Plots

The fan plots are interactive and include a slider which allows dynamic adjustment of the field.

```julia
using GLMakie

wavefan(W)
rayfan(W, system)
field_curves(W, system)
percent_distortion(W, system)
spot_diagram(W, system)
spot_diagram(full_trace(system, 0.0)) # using real rays
```

|    |    |
|----|----|
|![wavefan](assets/images/wavefan.png)|![rayfan](assets/images/rayfan.png)|
|![field_curves](assets/images/field_curves.png)|![percent_distortion](assets/images/percent_distortion.png)|
|![spot_diagram](assets/images/spot_diagram.png)|![real_spot_diagram](assets/images/real_spot_diagram.png)|
