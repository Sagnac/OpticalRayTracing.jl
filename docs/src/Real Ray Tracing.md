# Real Ray Tracing

This package provides several basic functions for tracing real rays.

```@setup real_ray_tracing

```

The following traces a real ray through the surfaces with initial height `y` and initial angle `U`.
```@example real_ray_tracing
y = 7.0
U = 0.05
raytrace(surfaces, y, U, RealRay)
```

And this finds the real marginal ray for the given surfaces and system.
```@example real_ray_tracing
trace_marginal_ray(surfaces, system)
```

In addition, the amount of spherical aberration produced from a real ray trace can be quantifed.

The following computes the transverse ray errors and returns a vector with the ray heights at the entrance pupil and the ray errors at the paraxial plane.

```@example real_ray_tracing
y, ε = TSA(surfaces, system)
ε[end]
```

These can then be used to approximate, using linear least squares, the coefficients to a tangential ray aberration expansion of the form:

`ε = By³ + Cy⁵ + Dy⁷ …`

```@example real_ray_tracing
degree = 9
SA(y, ε, degree)
```
