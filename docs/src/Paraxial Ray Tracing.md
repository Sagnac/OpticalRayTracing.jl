# Paraxial Ray Tracing

The `raytrace` function provides methods for tracing over a given system or lens while the `transfer` function provides convenience methods for working with transfer matrices.

```@setup ray_tracing

```

Lenses can be computed either by passing a surface matrix, as in `Lens(surfaces)`, or directly by using the `Lens(lens_matrix, refractive_indices)` constructor.

```@example ray_tracing
# telephoto objective comprised of two thin lenses in air
lens = Lens([
    # t/n   # power
    0.0      0.02
    25.0    -0.02
], ones(3))

y, nu = 1.0, 0.0

raytrace(lens, y, nu)
```

```@example ray_tracing
# object height
h = 10.0
# signed object distance from vertex
s = -100.0

raytrace(system, h, s)
```

```@example ray_tracing
(; f) = system
FFL, BFL = system.EFFD, system.EBFD
# positive distance from front vertex
τ = f + abs(FFL)
# positive distance from back vertex
τ′ = f + BFL
u = 0.07

# 2f-2f system
transfer(system.M, [0.0, u], τ, τ′)
```

In addition, general transfer matrices can be directly constructed by passing a `Lens` to the `TransferMatrix` constructor.
