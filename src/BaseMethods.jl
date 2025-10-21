function Base.show(io::IO, system::T) where T <: System
    print(io, "f: ")
    show(IOContext(io, :compact => true), system.f)
end

function Base.show(io::IO, ::MIME"text/plain", system::T) where T <: System
    println(T)
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

function Base.show(io::IO, m::MIME"text/plain", ray::Ray)
    summary(io, ray)
    println(io, ".yu:")
    show(IOContext(io, :displaysize => displaysize(io) .- (1, 0)), m, ray.yu)
end

function Base.show(io::IO, ::MIME"text/plain", aberr::T) where T <: Aberration
    println(T)
    for property in fieldnames(T)
        property === :sagittal && break
        @printf(io, "\n%4s: %.4f", property, getproperty(aberr, property))
    end
    return
end

Base.getindex(M::TransferMatrix) = M.M
