# Ospa dist

function ospa_dist(pca1::Vector{Pointcloud}, 
                   pca2::Vector{Pointcloud}, 
                   c::S
                  ) where {S <: Real}
    dmat = Matrix{Float64}(length(pca1), length(pca2))
    for i=1:length(pca1)
        for j=1:length(pca2)
            dmat[i,j] = ospa_dist(pca1[i],pca2[j],c)
        end
    end
    dmat
end

function ospa_dist(pc1::Pointcloud, 
                   pc2::Pointcloud, 
                   c::S
                  ) where {S <: Real}
    if size(pc1)[1] > size(pc2)[1]
        return ospa_dist(pc2, pc1, c)
    end
    dmat = p2dist(pc1, pc2)
    assignments = hungarian(dmat)[1]
    cost = sum([min(dmat[i, assignments[i]], c) for i=1:size(pc1)[1] if assignments[i] != 0])
    cost + c*(size(pc2)[1] - size(pc1)[1])
end

function optimal_assignments(barycenter::Pointcloud, 
                             measurements::Vector{Pointcloud}
                            )
    map(x -> hungarian(p2dist(barycenter, measurements[x]))[1], 1:length(measurements))
end

function p2dist(x,y)
    [sqrt.(sum((x[i,:] .- y[j,:]).^2)) for i=1:size(x)[1], j=1:size(y)[1]]
end

function p2dist(x)
    p2dist(x,x)
end

