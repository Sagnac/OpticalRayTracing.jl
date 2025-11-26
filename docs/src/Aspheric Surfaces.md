# Aspheric Surfaces

There is support for tracing over aspheric surfaces by specifying the base sphere radius of curvature at the vertex `R`, the conic constant `K`, and any optional, additional surface functions (e.g. polynomials) `p`.

```@setup aspheric

```

This can be done either by passing keyword arguments for the latter two to the real `raytrace` method or by using the `Layout` constructor:

```@example aspheric
# example for a single surface, including object space
R = [Inf, -100.0]
t = [0.0, 0.0]
n = [1.0, -1.0] # reflector
K = [1.0, -0.7] # elliptical profile, prolate spheroid surface
p = [zero, y -> y ^ 4 / R[2] ^ 4] # fourth order surface sagittal correction
# this constructor also supports specifying the full 4 column matrix directly
# instead of keyword arguments
surfaces = Layout(; R, t, n, K, p)
```
