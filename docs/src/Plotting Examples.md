# Plotting Examples

This package defines several plot recipes for [Makie](https://github.com/MakieOrg/Makie.jl).

Supported keyword arguments are:
* `theme` (using Makie's theme customizability);
* `ray_colors` (a two-tuple of the marginal & chief ray colors);
* `surface_color` (for the optical element surfaces & image plane);
* anything Makie's `lines` accepts.

Using the mutating `rayplot!` version will draw on top of the current figure.

```julia
using GLMakie

theme = theme_black()
surface_color = :white
ray_colors = (:cyan, :red)

rayplot(surfaces, a, system; theme, surface_color, ray_colors)

rays = raytrace(system, -24.0, -1.5 * system.f)
fig = rayplot(rays)
rays = raytrace(system, 24.0, -1.5 * system.f)
rayplot!(rays)
```

![rayplot](assets/images/rayplot.png)
![rays](assets/images/rays.png)
