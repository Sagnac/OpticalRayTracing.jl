module OpticalRayTracing

using Printf, LinearAlgebra

export TransferMatrix,
       Lens,
       Marginal,
       Chief,
       ParaxialRay,
       Tangential,
       FOV,
       transfer,
       refract,
       refract!,
       reverse_transfer,
       raytrace,
       trace_marginal_ray,
       trace_chief_ray,
       scale!,
       solve,
       flatten,
       raypoints,
       rayplot,
       rayplot!,
       vignetting,
       aberrations,
       compute_surfaces,
       incidences,
       RayError,
       Sagittal,
       Skew,
       wavefan,
       rayfan,
       field_curves,
       percent_distortion,
       spot_size,
       RealRay,
       caustic,
       Layout,
       Vignetting,
       TSA,
       SA,
       sag,
       Spherical,
       Aspheric,
       RayBasis,
       full_trace,
       wavegrad,
       optimize

include("Types.jl")
include("RayTracing.jl")
include("TransferMatrix.jl")
include("Vignetting.jl")
include("SeidelAberrations.jl")
include("Optimization.jl")
include("FullRayTracing.jl")
include("BaseMethods.jl")
include("RayPlot.jl")
include("API.jl")

end
