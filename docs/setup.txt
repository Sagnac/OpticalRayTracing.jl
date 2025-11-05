using OpticalRayTracing

# Tessar lens, units in mm; prescription from Hecht's Optics
surfaces = [
    # R        t        n
    Inf        0.0      1.0
      16.28    3.57     1.6116
    -275.7     1.89     1.0
     -34.57    0.81     1.6053
      15.82    2.345    1.0
    Inf        0.905    1.0
    Inf        2.17     1.5123
      19.2     3.96     1.6116
     -24.0     0.0      1.0
]

# Clear aperture semi-diameters
a = [9.5, 9.5, 9.0, 9.0, 7.63, 8.5, 8.5, 8.5]

# image height
h′ = -21.5

system = solve(surfaces, a, h′)
