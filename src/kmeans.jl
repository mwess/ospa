# kmeans implementation for the ospa distance.


function ospa_kmeans(measurements::Vector{Pointcloud},
                     c::S
                     ;
                     nclusters::Integer=2,
                     niter::Integer=100,
                     eps::Float64=0.05,
                     verbose=false
                    ) where {S <: Real}
    # This has to go into the for loop
    clusters = init_barycenter_clusters(measurements, nclusters)
    costs = fill(Inf, nclusters)
    max_barycenter_rounds = 1
    old_indices = fill(Inf, length(measurements))
    for i=1:niter
        if verbose
            println("Iteration: $(i)")
        end
        dmat = ospa_dist(clusters, measurements, c)
        indices = map(x -> findmin(dmat[:,x])[2], 1:size(dmat)[2])
        for j =1:nclusters
            inds = indices.==j
            # this should be provided by the user or automatically derived
            weights = fill(1/sum(inds), sum(inds))
            clusters[j], costs[j] = ospa_barycenter(measurements[inds], c, weights, :binned)
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

function init_barycenter_clusters(datapoints::Vector{Pointcloud}, nclusters::Integer)
    datapoints[[1,end]]
end
