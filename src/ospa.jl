__precompile__()


module ospa

using Hungarian

export ospa_kmeans, ospa_barycenter, ospa_dist, Pointcloud

include("types.jl")
include("distance.jl")
include("barycenter.jl")
include("kmeans.jl")

end # module
