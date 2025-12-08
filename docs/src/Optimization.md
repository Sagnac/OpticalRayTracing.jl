# Optimization

Constraint based optimization can be carried out by specifying which surface parameters to vary and which system properties to aim for as well as which third order aberrations to minimize.

```@setup optimization

```

```@example optimization
# you can specify a vector of integers corresponding to linear indices,
# or a CartesianIndex vector
v = [2:5; [8, 9]] # this varies only the curvatures for the example lens
# maps System fields to desired values
constraints = Dict(:f => system.f)
# the following two are optional,
# if not provided a simple mean over the five Seidel coefficients is taken
aberr = [:W040, :W131, :W222]
weights = [2.0, 2.0, 1.0]
new_system = optimize(system, v, constraints, aberr, weights)
```

```@repl optimization
new_system.layout
w = aberrations(new_system)
```

Note that this employs a rudimentary optimization procedure meant to provide a rough initial numerical solution; as a result the constraint values might not be met exactly.
