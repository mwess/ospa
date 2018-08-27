import Base: length

mutable struct Group
    points
    weight
    Group(x) = new(x, 1.)
end

length(gr::Group) = size(gr.points)[1]

function get_data(gr::Group)
    mean(gr.points, dims=1)
end

function init_groups(pc::Pointcloud)
    Groups([Group(pc[i:i,:]) for i = 1:size(pc)[1]])
end

Groups = Vector{Group}

function get_matrix(groups::Groups)
    points = map(get_data, groups)
    [points[i][j] for i=1:length(groups), j=1:size(points[1])[2]]
end

function add_new_group!(groups::Groups, point)
    new_group = Group(point)
    push!(groups, new_group)
end

function add_points!(group::Group, point)
    group.points = [group.points ; point]
end

function mean_initialization(measurements::Vector{Pointcloud}, c::S) where {S <: Real}
    perminds = randperm(length(measurements))
    groups = init_groups(measurements[perminds[1]])
    addpointcounter = 0
    addgroupcounter = length(groups)
    for i in perminds[2:end]
        current_points = get_matrix(groups)
        assignments, _ = hungarian(p2dist(current_points, measurements[i]))
        second_greater_first = size(current_points) < size(measurements[i])
        for j=1:length(assignments)
            ind = assignments[j]
            if assignments[j] != 0
                if p2dist(current_points[j:j,:], measurements[i][ind:ind, :])[1] < c
                    add_points!(groups[j], measurements[i][ind:ind, :])
                    addpointcounter += 1
                else
                    add_new_group!(groups, measurements[i][ind:ind, :])
                end
            end
        end
        if second_greater_first
            remaining_inds = filter(x -> !(x in assignments), 1:size(measurements[i])[1])
            for rind in remaining_inds
                add_new_group!(groups, measurements[i][rind:rind,:])
                addgroupcounter += 1
            end
        end
    end
    #return groups
    npoints = map(length, groups) |> sum
    [(get_data(g),length(g)/npoints) for g in groups]
end

function sample_from_mean_initialization(sample_space::Array{Tuple{Array{Float64,2},Float64},1},
                                         npoints::Integer,
                                         nsamples::Integer
                                        )
    sample_space_mat, weights = zip(sample_space...) |> collect
    sample_space_mat = vcat(sample_space_mat...)
    weights = [x for x = weights]
    samples = []
    for i=1:nsamples
        idxs = draw_from_weighted_distribution_without_replacement(weights, npoints)
        push!(samples, sample_space_mat[idxs,:])
    end
    samples
end

function binned_barycenter_initialization(measurements::Vector{Pointcloud}
                                        ,weights::Vector{R}
                                        ,c::S
                                        ,min_card::Integer
                                        ,max_card::Integer
                                        ,eps::T
                                        ,nsamples::Integer
                                        ) where {R <: Real, S <: Real, T <: Real}
    sample_space = mean_initialization(measurements, c)
    if max_card == 0
        max_card = length(sample_space)
    end
    max_card = min(max_card, length(sample_space))
    Barycenters = Vector{Pointcloud}()
    costs = Vector{Float64}()
    for i = min_card:max_card
        best_cost = Inf
        best_bc = missing
        bc_inits = sample_from_mean_initialization(sample_space, i, nsamples)
        for _=1:nsamples
            bc_init = sample_from_mean_initialization(sample_space, i, 1)
            cur_bc, cur_cost = ospa_barycenter(measurements, weights, c, bc_init[1], eps)
            if cur_cost < best_cost
                best_cost = cur_cost
                best_bc = cur_bc
            end
        end
#        for sample in bc_inits
#            cur_bc, cur_cost = ospa_barycenter(measurements, weights, c, sample, eps)
#            if cur_cost < best_cost
#                best_cost = cur_cost
#                best_bc = cur_bc
#            end
#        end
        push!(Barycenters, best_bc)
        push!(costs, best_cost)
    end
    (_, ind) = findmin(costs)
    Barycenters[ind], costs[ind]
end
