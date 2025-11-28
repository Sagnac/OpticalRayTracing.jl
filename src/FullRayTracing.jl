function refract!(k::Vector, m::Vector, n1::Float64, n2::Float64)
    η = n1 / n2
    γ = -k ⋅ m
    Δ = 1 - η ^ 2 * (1 - γ ^ 2)
    if Δ ≥ 0.0
        @. k = abs(η) * k + (η * γ - sign(η) * sqrt(Δ)) * m
    else
        @. k += 2.0 * γ * m # TIR
    end
    return k
end

function raytrace(surfaces::AbstractMatrix, y, x, U, V, ::Type{Vector{RealRay}};
                  K = zeros(size(surfaces, 1)), p = fill(zero, length(K)))
    R, t, n = eachcol(surfaces)
    tsy = copy(t)
    tsx = copy(t)
    u = tan(U)
    v = tan(V)
    k = [v, u, 1.0]
    normalize!(k) # direction cosines
    m = [0.0, 0.0, -1.0] # surface normal
    for i = 1:size(surfaces, 1)-1
        y += u * tsy[i]
        x += v * tsx[i]
        Rs = R[i+1]
        Ks = K[i+1]
        ps = p[i+1]
        sy = sag(y, U, Rs, Ks, ps)
        sx = sag(x, V, Rs, Ks, ps)
        y += sy * u
        x += sx * v
        tsy[i] += sy
        tsx[i] += sx
        tsy[i+1] -= sy
        tsx[i+1] -= sx
        ds_dy = tilt(y, Rs, Ks, ps)
        ds_dx = tilt(x, Rs, Ks, ps)
        m .= ds_dx, ds_dy, -1.0
        normalize!(m)
        refract!(k, m, n[i], n[i+1])
        u = k[2] / k[3]
        v = k[1] / k[3]
        U = atan(u)
        V = atan(v)
    end
    return x, y
end
