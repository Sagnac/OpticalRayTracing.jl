abstract type Paraxial end

abstract type Fundamental <: Paraxial end

struct Tangential <: Paraxial end

struct Sagittal <: Paraxial end

struct Skew <: Paraxial end

struct Marginal <: Fundamental end

struct Chief <: Fundamental end

struct Ray{T <: Paraxial}
    y::Vector{Float64}
    n::Vector{Float64}
    u::Vector{Float64}
    yu::Matrix{Float64}
    nu::Vector{Float64}
    ynu::Matrix{Float64}
    z::Vector{Float64}
    function Ray{T}(ynu, τ, n) where T <: Paraxial
        y, nu = eachcol(ynu)
        n = [n; n[end]]
        u = map(/, nu, n)
        yu = [y u]
        t = map(*, τ, @view(n[1:end-2]))
        if T <: Fundamental
            t = [t; -y[end-1] / u[end-1]]
        end
        z = cumsum(t)
        z0 = iszero(u[1]) ? -(extrema(z)...) * 0.1 : -y[2] / u[1]
        z = [z0; z]
        new(y, n, u, yu, nu, ynu, z)
    end
end

struct RealRay
    y::Vector{Float64}
    u::Vector{Float64}
    yu::Matrix{Float64}
    n::Vector{Float64}
    z::Vector{Float64}
    RealRay(yu, t, n) = new(eachcol(yu)..., yu, n, cumsum(t))
end

struct RayBasis
    marginal::Ray{Marginal}
    chief::Ray{Chief}
    H::Float64
end

struct TransferMatrix <: AbstractArray{Float64, 2}
    M::Matrix{Float64}
end

struct Lens <: AbstractArray{Float64, 2}
    M::Matrix{Float64}
    n::Vector{Float64}
end

struct Pupil
    D::Float64
    t::Float64
end

struct System
    f::Float64
    EBFD::Float64
    EFFD::Float64
    N::Float64
    FOV::Float64
    stop::Int
    EP::Pupil
    XP::Pupil
    marginal::Ray{Marginal}
    chief::Ray{Chief}
    trace::Matrix{Float64}
    H::Float64
    P1::Float64
    P2::Float64
    PN::Float64
    M::TransferMatrix
    lens::Lens
    a::Vector{Float64}
end

struct Aberration
    W040::Float64
    W131::Float64
    W222::Float64
    W220::Float64 # W220S
    W311::Float64
    W020::Float64 # axial color
    W111::Float64 # lateral color
    W220P::Float64
    W220M::Float64
    W220T::Float64
    spherical::Vector{Float64}
    coma::Vector{Float64}
    astigmatism::Vector{Float64}
    sagittal::Vector{Float64}
    distortion::Vector{Float64}
    axial::Vector{Float64}
    lateral::Vector{Float64}
    petzval::Vector{Float64}
    medial::Vector{Float64}
    tangential::Vector{Float64}
    λ::Float64
    field_sign::Int
end

struct RayError{T <: Paraxial}
    W::Aberration
    nu::Float64
    field_sign::Int
end

const SystemOrRayBasis = Union{System, RayBasis}

const LensOrTransferMatrix = Union{Lens, TransferMatrix}
