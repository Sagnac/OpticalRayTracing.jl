abstract type ParaxialRay end

abstract type TangentialRay <: ParaxialRay end

struct Marginal <: TangentialRay end

struct Chief <: TangentialRay end

struct Ray{T <: TangentialRay}
    y::Vector{Float64}
    n::Vector{Float64}
    u::Vector{Float64}
    yu::Matrix{Float64}
    nu::Vector{Float64}
    ynu::Matrix{Float64}
    function Ray{T}(ynu, n) where T <: TangentialRay
        y, nu = eachcol(ynu)
        n = [n; n[end]]
        u = map(/, nu, n)
        yu = [y u]
        new(y, n, u, yu, nu, ynu)
    end
end

struct TransferMatrix
    M::Matrix{Float64}
end

struct Lens
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
