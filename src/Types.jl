abstract type ParaxialRay end

struct Tangential <: ParaxialRay end

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
    M::TransferMatrix
    lens::Lens
end

struct Aberration
    spherical::Float64
    coma::Float64
    astigmatism::Float64
    petzval::Float64
    distortion::Float64
    axial::Float64
    lateral::Float64
    sagittal::Float64
    medial::Float64
    tangential::Float64
    W040::Vector{Float64}
    W131::Vector{Float64}
    W222::Vector{Float64}
    W220P::Vector{Float64}
    W220::Vector{Float64}
    W311::Vector{Float64}
    W020::Vector{Float64}
    W111::Vector{Float64}
end

const SystemOrRayBasis = Union{System, RayBasis}

const LensOrTransferMatrix = Union{Lens, TransferMatrix}

const MatrixOrTransferMatrix = Union{Matrix, TransferMatrix}
