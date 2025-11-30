abstract type AbstractRay end

abstract type FundamentalRay <: AbstractRay end

abstract type Profile end

struct Tangential <: AbstractRay end

struct Sagittal <: AbstractRay end

struct Skew <: AbstractRay end

struct Marginal <: FundamentalRay end

struct Chief <: FundamentalRay end

struct Spherical <: Profile end

struct Aspheric <: Profile end

struct Polynomial{F <: Function}
    f::F
end

(p::Polynomial)(x) = p.f(x)

fill_poly(x) = fill(Polynomial(zero), x)

struct ParaxialRay{T <: AbstractRay}
    y::Vector{Float64}
    n::Vector{Float64}
    u::Vector{Float64}
    yu::Matrix{Float64}
    nu::Vector{Float64}
    ynu::Matrix{Float64}
    z::Vector{Float64}
    function ParaxialRay{T}(ynu, τ, n) where T <: AbstractRay
        y, nu = eachcol(ynu)
        n = [n; n[end]]
        u = map(/, nu, n)
        yu = [y u]
        t = map(*, τ, @view(n[1:end-2]))
        if T <: FundamentalRay
            t = [t; -y[end-1] / u[end-1]]
        end
        z = cumsum(t)
        z0 = iszero(u[1]) ? -(extrema(z)...) * 0.1 : -y[2] / u[1]
        z = [z0; z]
        new(y, n, u, yu, nu, ynu, z)
    end
end

struct RealRay{T <: AbstractRay}
    y::Vector{Float64}
    u::Vector{Float64}
    yu::Matrix{Float64}
    n::Vector{Float64}
    z::Vector{Float64}
end

function RealRay{T}(yu, t, n) where T <: AbstractRay
    RealRay{T}(eachcol(yu)..., yu, n, cumsum(t))
end

struct RayBasis
    marginal::ParaxialRay{Marginal}
    chief::ParaxialRay{Chief}
    H::Float64
end

struct TransferMatrix <: AbstractArray{Float64, 2}
    M::Matrix{Float64}
end

struct Lens <: AbstractArray{Float64, 2}
    M::Matrix{Float64}
    n::Vector{Float64}
end

@kwdef struct Layout{T <: Profile} <: AbstractArray{Float64, 2}
    M::Matrix{Float64} = Matrix{Float64}(undef, 0, 0)
    R::Vector{Float64}
    t::Vector{Float64}
    n::Vector{Float64}
    K::Vector{Float64} = Float64[]
    p::Vector{Polynomial} = Polynomial{typeof(zero)}[]
    function Layout{T}(M, R, t, n, K, p) where T <: Profile
        isempty(K) && (K = zeros(length(R)))
        isempty(p) && (p = fill_poly(length(R)))
        typeof(p) <: Vector{<:Polynomial} || (p = Polynomial.(p))
        new{T}([R t n K], R, t, n, K, p)
    end
end

function Layout(M)
    x = size(M, 1)
    K = zeros(x)
    p = fill_poly(x)
    Layout{Spherical}(M, eachcol(M)..., K, p)
end

Layout{Aspheric}(M) = Layout{Aspheric}(M, eachcol(M)..., fill_poly(size(M, 1)))

function Layout(R, t, n)
    x = length(R)
    Layout{Spherical}([R t n], R, t, n, zeros(x), fill_poly(x))
end

function Layout(M, R, t, n, K, p)
    Layout{iszero(K) && all(==(zero), p) ? Spherical : Aspheric}(M, R, t, n, K, p)
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
    marginal::ParaxialRay{Marginal}
    chief::ParaxialRay{Chief}
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

struct Vignetting
    M::Matrix{Float64}
    FOV::Matrix{Float64}
    un::Bool
    limit::Vector{Int}
    partial::Vector{Int}
    full::Vector{Int}
end

struct RayError{T <: AbstractRay}
    W::Aberration
    nu::Float64
    field_sign::Int
end

struct RealRayError
    x::Matrix{Float64}
    y::Matrix{Float64}
    nu::Float64
    r::Vector{Float64}
    t::Vector{Float64}
end

const SystemOrRayBasis = Union{System, RayBasis}

const OpticalMatrix = Union{Lens, Layout, TransferMatrix}
