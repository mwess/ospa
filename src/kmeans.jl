# kmeans implementation for the ospa distance.


function ospa_kmeans(measurements::Vector{Pointcloud}, 
                     c::S
                     ; 
                     nclusters::Integer=2, 
                     niter::Integer=100, 
                     eps::Float64=0.05, 
                     verbose=false
                    ) where {S <: Real}
    clusters = init_barycenter_clusters(measurements, nclusters)
    costs = fill(Inf, nclusters)
    max_barycenter_rounds = 1
    old_indices = fill(Inf, length(measurements))
    for i=1:niter
        if verbose
            println("Iteration: $(i)")
        end
        dmat = ospa_dist(clusters, measurements, c)
        indices = map(x -> indmin(dmat[:,x]), 1:size(dmat)[2])
        for j =1:nclusters
            inds = indices.==j
            weights = fill(1/sum(inds), sum(inds))
            clusters[j], costs[j] = ospa_barycenter(measurements[inds], weights, c, max_barycenter_rounds, eps)
            if verbose
                println("Cost $(j): $(costs[j])")
            end
        end
        # Mechanism for ending prematurely.
        if old_indices == indices
            break
        end
        old_indices = indices
    end
    clusters, costs
end

# This is only for an example. Do something besser here.
function init_barycenter_clusters(measurements::Vector{Pointcloud}, 
                                  nclusters::Integer
                                 )
    init_clusters = Vector{Pointcloud}(nclusters)
    init_clusters[1] = measurements[1]
    init_clusters[end] = measurements[end]
    init_clusters
end
