abstract type ParaxialRay end

struct Tangential <: ParaxialRay end

struct Sagittal <: ParaxialRay end

struct Marginal <: ParaxialRay end

struct Chief <: ParaxialRay end

struct Ray{T <: ParaxialRay}
    y::Vector{Float64}
    n::Vector{Float64}
    u::Vector{Float64}
    yu::Matrix{Float64}
    nu::Vector{Float64}
    ynu::Matrix{Float64}
    z::Vector{Float64}
    function Ray{T}(ynu, τ, n) where T <: ParaxialRay
        y, nu = eachcol(ynu)
        n = [n; n[end]]
        u = map(/, nu, n)
        yu = [y u]
        t = map(*, τ, @view(n[1:end-2]))
        z = [-y[2] / u[1]; cumsum([t; -y[end-1] / u[end-1]])]
        new(y, n, u, yu, nu, ynu, z)
    end
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
end

const OrthogonalRay = Union{Tangential, Sagittal}

struct RayError{T <: OrthogonalRay}
    W::Aberration
    nu::Float64
end

const SystemOrRayBasis = Union{System, RayBasis}

const LensOrTransferMatrix = Union{Lens, TransferMatrix}
