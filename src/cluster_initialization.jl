# Implementation of the kmeans++ initialization for the initial barycenters
function init_barycenter_pp(measurements::Vector{Pointcloud}, 
                            weights::Vector{R}, 
                            c::S, 
                            nclusters::Integer,
                            kwargs...
                           ) where{R <: Real, S <: Real}
    clusters = Vector{Pointcloud}([rand(measurements)])
    for _=2:nclusters
        dmat = ospa_dist(clusters, measurements, c)
        min_distances = vcat(minimum(dmat, 1).^2...)
        min_distances /= sum(min_distances)
        next_inds = draw_from_weighted_distribution(min_distances, 1)
        for ind in next_inds
            push!(clusters, measurements[ind])
        end
    end
    clusters
end

function draw_from_weighted_distribution(normalized_dist::Vector{S},
                                         npoints::Integer,
                                        ) where {S <: Real}
    cdist = cumsum(normalized_dist)
    map(x -> find_group(cdist, x), rand(npoints))
end

function draw_from_weighted_distribution_without_replacement(normalized_dist::Vector{S},
                                                             npoints::Integer,
                                                            ) where {S <: Real}
    if npoints > length(normalized_dist)
        error("Can't draw more points than are available!")
    end
    cdist = cumsum(normalized_dist)
    res = Vector{Integer}()
    while length(res) < npoints
        point = find_group(cdist, rand())
        if !(point in res)
            push!(res, point)
        end
    end
    res
end

function find_group(cum_dist::Vector{S}, point::S) where {S <: Real}
    if 0 <= point < cum_dist[1]
        return 1
    end
    for i=2:length(cum_dist)
        if cum_dist[i-1] <= point < cum_dist[i]
            return i
        end
    end
end

# Not needed anymore :(
function get_cumulative_distribution(vec::Vector{S}) where {S <: Real}
    cdist = copy(vec)
    for i=2:length(vec)
        cdist[i] = vec[i] + cdist[i-1]
    end
    cdist
end


# This is only for the example test data.
function init_for_example(measurements::Vector{Pointcloud}, 
                                  nclusters::Integer
                                 )
    init_clusters = Vector{Pointcloud}(nclusters)
    init_clusters[1] = measurements[1]
    init_clusters[end] = measurements[end]
    init_clusters
end


# Just random initialisation
function init_randomly(measurements::Vector{Pointcloud},
                       nclusters::Integer)
    indices = randperm(length(measurements))[1:nlusters]
    measurements[indices]
end

