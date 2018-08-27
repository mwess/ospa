__precompile__()


module ospa

using Hungarian
using LinearAlgebra
using Statistics
using Random

export ospa_kmeans, ospa_barycenter,
       ospa_dist, Pointcloud,
       generateData
       #mean_initialization, sample_from_mean_initialization,
       #ospa_barycenter_sol2

include("types.jl")
include("distance.jl")
include("barycenter.jl")
include("mean_barycenter_initialization.jl")
include("recursive_barycenter_initialization.jl")
include("cluster_initialization.jl")
include("kmeans.jl")
include("testdata.jl")

end # module
