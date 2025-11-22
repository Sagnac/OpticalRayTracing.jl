import Base: show, getindex, setindex!, size, IndexStyle, iterate

function show(io::IO, system::T) where T <: System
    print(io, "f: ")
    show(IOContext(io, :compact => true), system.f)
end

function show(io::IO, ::MIME"text/plain", system::T) where T <: System
    println(io, T)
    if haskey(io, :typeinfo)
        show(io, system)
        return
    end
    for property in fieldnames(T)
        property === :marginal && break
        value = getproperty(system, property)
        @printf(io, "\n%4s: ", property)
        if typeof(value) === Pupil
            @printf(io, "D = %.4f, t = %.4f", value.D, value.t)
        elseif property === :stop
            @printf(io, "%d", value)
        else
            @printf(io, "%.4f", value)
        end
    end
    return
end

show(io::IO, ray::Union{ParaxialRay, RealRay}) = summary(io, ray)

function show(io::IO, m::MIME"text/plain", ray::Union{ParaxialRay, RealRay})
    show(io, ray)
    haskey(io, :typeinfo) && return
    println(io, ".yu:")
    show(IOContext(io, :displaysize => displaysize(io) .- (1, 0)), m, ray.yu)
end

show(io::IO, ray_basis::RayBasis) = summary(io, ray_basis)

function show(io::IO, m::MIME"text/plain", ray_basis::T) where T <: RayBasis
    show(io, ray_basis)
    haskey(io, :typeinfo) && return
    print(io, '\n', fieldtypes(T)[1:2])
end

show(io::IO, aberr::Aberration) = summary(io, aberr)

function show(io::IO, ::MIME"text/plain", aberr::T) where T <: Aberration
    show(io, aberr)
    haskey(io, :typeinfo) && return
    println(io)
    for property in fieldnames(T)
        property === :W220P && break
        @printf(io, "\n%4s: %.4f", property, getproperty(aberr, property))
    end
    return
end

show(io::IO, vig::Vignetting) = summary(io, vig)

function show(io::IO, mime::MIME"text/plain", vig::Vignetting)
    (; M, un, partial, full, limit) = vig
    FOVs, _, heights = eachcol(vig.FOV)
    bold = true
    if un
        printstyled(io, "Unvignetted.\n\n"; bold, color = :green)
    else
        if !isempty(partial)
            printstyled(io, "Partially vignetted:\n"; bold, color = :yellow)
            show(io, partial)
            println(io, '\n')
        end
        if !isempty(full)
            printstyled(io, "Fully vignetted:\n"; bold, color = :magenta)
            show(io, full)
            println(io, '\n')
        end
        if !isempty(limit)
            printstyled(io, "Stop limited:\n"; bold, color = :red)
            show(io, limit)
            println(io, '\n')
        end
    end
    println(io, "Maximum supported FOVs:")
    @printf(io, "Unvignetted: %.4f°, h′: %.4f\n", FOVs[1], heights[1])
    @printf(io, "Half-vignetted: %.4f°, h′: %.4f\n", FOVs[2], heights[2])
    @printf(io, "Fully-vignetted: %.4f°, h′: %.4f\n\n", FOVs[3], heights[3])
    s1, s2, s3, s4 = (' ' ^ i for i = (7, 4, 7, 5))
    printstyled(io, "a", s1, "limit", s2, "un", s3, "half", s4, "full\n";
                bold, color = :cyan)
    show(io, mime, M)
end

getindex(M::TransferMatrix) = M.M

getindex(M::OpticalMatrix, i::Int) = M.M[i]

setindex!(M::Union{Lens, Layout}, v, i::Int) = M.M[i] = v

size(M::OpticalMatrix) = size(M.M)

IndexStyle(::Type{<:OpticalMatrix}) = IndexLinear()

function getindex(rays::T, i::Int) where T <: RayBasis
    i > 2 ? throw(BoundsError(rays, i)) : getfield(rays, i)
end

function iterate(rays::T, i::Int = 1) where T <: RayBasis
    i > 2 ? nothing : (getfield(rays, i), i + 1)
end
