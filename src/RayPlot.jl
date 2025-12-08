# the actual plotting functions are in ext/MakieExtension.jl
# these return the plot points

raypoints(system::SystemOrRayBasis) = raypoints(system.marginal, system.chief)

function raypoints(marginal::ParaxialRay{Marginal}, chief::ParaxialRay{Chief})
    (; z, n) = marginal
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
    ȳo1 = ȳ1 + nū * z[1]
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
function spot_diagram end
function caustic end
