# the actual plotting functions are in ext/MakieExtension.jl
# these return the plot points

raypoints(system::SystemOrRayBasis) = raypoints(system.marginal, system.chief)

function raypoints(marginal::Ray{Marginal}, chief::Ray{Chief})
    (; z, n) = marginal
    l = abs(z[end-1] - z[begin+1])
    BFD = isfinite(z[end]) ? abs(z[end]) : l / 4
    l += BFD
    if isinf(z[1])
        s = 0.1
        d = -l * s / n[1]
        z = [d; @view(z[2:end])]
    else
        d = z[1]
    end
    y0 = zero(z)
    if iszero(marginal.u[1])
        y1 = marginal.y
    else
        y1 = copy(marginal.y)
        y1[1] = 0.0
    end
    y2 = -y1
    nū = chief.nu[1]
    ȳ1 = chief.y[2]
    ȳo1 = ȳ1 + nū * d
    ȳ = [ȳo1; @view(chief.y[begin+1:end])]
    y3 = ȳ + y1
    y4 = ȳ + y2
    return z, [y0, y1, y2, ȳ, y3, y4]
end

# requires Makie
function rayplot end
function rayplot! end
function wavefan end
function rayfan end
function field_curves end
function percent_distortion end
function spot_size end
