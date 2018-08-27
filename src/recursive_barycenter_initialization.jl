
function iterative_barycenter_initialization(measurements::Vector{Pointcloud}
                                            ,weights::Vector{R}
                                            ,c::S
                                            ,eps::Float64
                                            ) where {R <: Real, S <: Real}
    #ospa_iterative_sampling_space(measurements, weights, c, eps)
    ospa_barycenter_sol2(measurements, weights, c, eps)
end



function ospa_barycenter_sol2(measurements::Vector{Pointcloud},
                              weights::Vector{R},
                              c::S,
                              eps::Float64
                             ) where {R <: Real, S <: Real}
    #assignments = Vector(map(x -> zeros(Integer, size(x)[1]), measurements))
    assignments = Vector(map(x -> fill(0, size(x)[1]), measurements))
    first_guess = init_by_point(measurements, assignments)
    ospa_barycenter_sol2(measurements, assignments, weights, c, first_guess, eps)
end

#function ospa_barycenter_sol2(measurements::Vector{Pointcloud},
#                              assignments::Vector{Vector{Integer}},
#                              weights::Vector{R},
#                              c::S,
#                              init_first_guess::Function,
#                              ;
#                              do_postprocessing=false,
#                              confidence_threshold=0.05,
#                             ) where {R <: Real, S <: Real}
#    # Start with first guess
#    ospa_barycenter_sol2(measurements, assignments, weights, c, init_fun, bc)
#end


function ospa_iterative_sampling_space(measurements::Vector{Pointcloud}
                                       ,weights::Vector{R}
                                       ,c::S
                                       ,eps::Float64
                                      ) where {R <: Real, S <: Real}
    cost_old = Inf
    #eps = 0.05
    assignments = []
    converged = false
    while !converged
        bc, cost, assignments, subcluster_cost = update_barycenter_sol2(bc, measurements, weights, c)
        if cost_old - cost < eps
            converged = true
        end
        cost_old = cost
    end

    bc, cost, assignments, subcluster_cost = add_remaining_points(bc, measurements, assignments, weights, c, subcluster_cost)
    sorting = sortperm(subcluster_cost)
    bc = bc[sorting, :]
    subcluster_cost = subcluster_cost[sorting, :]
    subcluster_cost = subcluster_cost / sum(subcluster_cost)
    bc, subcluster_cost
end

function ospa_barycenter_sol2(measurements::Vector{Pointcloud},
                              assignments::Vector{Vector{I}},
                              weights::Vector{R},
                              c::S,
                              bc::Pointcloud,
                              eps::Float64
                             ) where {R <: Real, S <: Real, I <: Integer}
    cost_old = Inf
    #eps = 0.05
    assignments = []
    converged = false
    while !converged
        bc, cost, assignments, subcluster_cost = update_barycenter_sol2(bc, measurements, weights, c)
        if cost_old - cost < eps
            converged = true
        end
        cost_old = cost
    end

    # Old
    #bc, cost, assignments, subcluster_cost = update_barycenter_sol2(bc, measurements, weights, c)
    #println("first bc")
    #println(bc)
    bc, cost, assignments, subcluster_cost = add_remaining_points(bc, measurements, assignments, weights, c, subcluster_cost)
    #if do_postprocessing
    #    println(subcluster_cost)
    #    bc, cost, subcluster_cost = confidence_based_postprocessing(bc, measurements, assignments, cost, confidence_threshold, c, subcluster_cost)
    #end
    #now select optimal barycenter.
    #println("Subcluster cost: $(subcluster_cost)")
    sorting = sortperm(subcluster_cost)
    bc = bc[sorting, :]
    subcluster_cost = subcluster_cost[sorting, :]
    cost_rising = false
    first_ind = 1
    last_ind = 1
    previous_cost = Inf
    current_cost = Inf
    best_bc = bc[first_ind:last_ind, :]
    best_cost = mospa(best_bc, measurements, weights, c)
    #println("Now select best barycenter.")
    while !cost_rising && last_ind <= size(bc)[1]
        current_cost = mospa(bc[first_ind:last_ind, :], measurements, weights, c)
#        println("Current_cost: $(current_cost).")
 #       println(bc[first_ind:last_ind, :])
        #println("Current cost for size $(last_ind): $(current_cost).\n Previous cost: $(previous_cost).")
        if current_cost > previous_cost || last_ind >= size(bc)[1]
        #if last_ind >= size(bc)[1]
            cost_rising = true
        else
            last_ind += 1
            previous_cost = current_cost
        end
    end
    best_bc = bc[first_ind:last_ind, :]
    best_bc, current_cost
end

function add_remaining_points(bc::Pointcloud,
                              measurements::Vector{Pointcloud},
                              assignments::Vector{Array{Integer,1}},
                              weights::Vector{R},
                              c::S,
                              subcluster_cost::Vector{Float64},
                             ) where {R <: Real, S <: Real}
    if map(x -> any(x .== 0), assignments) |> any
        exhausted = false
    else
        exhausted = true
    end
    idx_next_sub_cluster = size(bc)[1]
    loop_counter = 1
    while !exhausted
        #println("loop_counter: $(loop_counter).")
        loop_counter += 1
        idx_next_sub_cluster += 1
        next_point, outer_ind, inner_ind = draw_from_remaining_pool(measurements, assignments)
        #println("Next point drawn: $(next_point), outer_ind: $(outer_ind), inner_ind: $(inner_ind).")
        #next_sub_cluster = zeros(next_point)
        #old_next_sub_cluster = zeros(next_point)
        next_sub_cluster = fill(0, size(next_point))
        old_next_sub_cluster = fill(0, size(next_point))
        tmp_assignment_log = []
        tmp_sub_cluster_cost = 0
        converged = false
        while !converged
            opt_ass = selective_optimal_assignments(next_point, measurements, assignments)
            #println("assignment: $(opt_ass).")
            new_weight = 0
            tmp_assignment_log = []
            next_sub_cluster = fill(0, size(next_point))
            tmp_sub_cluster_cost = 0
            for i=1:length(opt_ass)
                idx = opt_ass[i][1]
                if idx != 0
                    #println("Next point: $(next_point), meas: $(measurements[i][idx:idx, :]), dist: $(p2dist(next_point, measurements[i][idx:idx, :])[1]).")
                    if p2dist(next_point, measurements[i][idx:idx, :])[1] <= c
                        next_sub_cluster += measurements[i][idx:idx, :]*weights[i]/max(length(measurements), length(assignments[i]))
                        new_weight += weights[i]/max(length(measurements), length(assignments[i]))
                        push!(tmp_assignment_log, (i, idx, idx_next_sub_cluster))
                        tmp_sub_cluster_cost += 1
                    end
                end
            end
            if new_weight > 0
                next_sub_cluster /= new_weight
            end
            if p2dist(next_sub_cluster, old_next_sub_cluster)[1] == 0
                converged = true
            else
                old_next_sub_cluster = next_sub_cluster
            end
        end
        bc = vcat(bc, next_sub_cluster)
        #println("Next bc:")
        #println(bc)
        push!(subcluster_cost, tmp_sub_cluster_cost)
        #println("New assignments: $(tmp_assignment_log).")
        copy_new_assignments!(assignments, tmp_assignment_log)
        #println("Size: $(size(bc)).")
        exhausted = check_for_remaining_points(assignments)
    end
    bc, mospa(bc, measurements, weights, c), assignments, subcluster_cost
end

function confidence_based_postprocessing(bc::Pointcloud,
                                         measurements::Vector{Pointcloud},
                                         assignments::Vector{Vector{Integer}},
                                         cost::R,
                                         confidence_threshold::Real,
                                         c::S,
                                         subcluster_cost::Vector{Float64},
                                         ;
                                         niter=20
                                        ) where {S <: Real, R <: Real}
    #Maximum sized barcyenter computed. Calc the confidence parameters.
    confidences, points_in_radius = compute_confidence(bc, measurements, c)
    clusters_to_shuffle = find(confidences .< confidence_threshold)
    # Find all clusters that need to be reshuffled now.
    clusters_to_shuffle_final = copy(clusters_to_shuffle)
    for i=1:length(points_in_radius)
        for j=1:length(points_in_radius[i])
            if any(k in points_in_radius[i][j] for k in clusters_to_shuffle)
                clusters_to_shuffle_final = unique([points_in_radius[i][j]; clusters_to_shuffle_final])
            end
        end
    end
    #unassign all points that need to be reshuffled.
    for cluster in clusters_to_shuffle_final
        subcluster_cost[cluster] = 0
    end
    for i =1:length(assignments)
        for j=1:length(assignments[i])
            if assignments[i][j] in clusters_to_shuffle_final
                assignments[i][j] = 0
            end
        end
    end
    #Now recompute the barycenter
    bc_best = copy(bc)
    cost_best = copy(cost)
    for _=1:niter
        bc_new, cost_new, _, subcluster_cost = add_remaining_points(bc, measurements, assignments, weights, c, subcluster_cost)
        if cost_new < cost_best
            bc_best = bc_new
            cost_best = cost_new
        end
    end
    bc_best, cost_best, subcluster_cost
end

function compute_confidence(bc::Pointcloud,
                            measurements::Vector{Pointcloud},
                            c::S
                           ) where {S <: Real}
    #Compute for eacht point to how many sub clusters it can belong
    counters =  map(x -> fill(Vector{Integer}(), size(x)[1]), measurements)
    #point_confidences = [ones(size(counters[i])) for i=1:length(counters)]
    #point_counter_for_clusters = [zeros(size(counters[i])) for i=1:length(counters)]
    cluster_confidence = zeros(size(bc)[1])
    subclusters_cardinality = zeros(size(bc)[1])
    for i=1:length(measurements)
        dmat = p2dist(measurements[i], bc)
        for j=1:size(measurements[i])[1]
            points_in_radius = find(dmat[j,:] .< c)
            counters[i][j] = [counters[i][j] ; points_in_radius]
            for k=1:length(points_in_radius)
                cluster_confidence[points_in_radius[k]] += 1/length(counters[i][j])
                #point_counter_for_clusters[points_in_radius[k]] += 1
                subclusters_cardinality[points_in_radius[k]] += 1
            end
        end
    end
    #println(point_counter_for_clusters)
    #println(cluster_confidence)
    for i=1:length(cluster_confidence)
        #cluster_confidence[i] /= point_counter_for_clusters[i]
        cluster_confidence[i] /= subclusters_cardinality[i]
    end
    cluster_confidence, counters
end

function copy_new_assignments!(assignments::Vector{Vector{Integer}}, new_assignments)
    for triple in new_assignments
        assignments[triple[1]][triple[2]] = triple[3]
    end
end


function update_barycenter_sol2(cur_bc::Pointcloud,
                                measurements::Vector{Pointcloud},
                                weights::Vector{R},
                                c::S,
                               ) where {R <: Real, S <: Real}
    #updated_bc = Matrix{eltype(cur_bc)}(0,size(cur_bc)[2])
    updated_bc = Matrix{eltype(cur_bc)}(fill(0, 0, size(cur_bc)[2]))
    subcluster_weight = Vector{Float64}()
    assignments = optimal_assignments(cur_bc, measurements)
    assigned_clusters = map(x -> zeros(Integer, size(x)[1]), measurements)
    #not sure if weights are handled correctly here
    for idxp = 1:size(cur_bc)[1]
        assigned_points = zeros(eltype(cur_bc), 1,size(cur_bc)[2])
        new_weight = 0
        tmp_subcluster_weight = 0
        for i=1:length(assignments)
            idxa = assignments[i][idxp]
            if idxa != 0 && p2dist(cur_bc[idxp:idxp,:], measurements[i][idxa:idxa,:])[1] <= c
                # Is the correct weight index idxa or i???
                tmp_subcluster_weight += 1
                #is it correct to change that line?
                #assigned_points += measurements[i][idxa:idxa,:]*weights[idxa]/max(length(measurements), length(assignments[i]))
                assigned_points += measurements[i][idxa:idxa,:]*weights[i]/max(length(measurements), length(assignments[i]))
                new_weight += weights[i]/max(length(measurements), length(assignments[i]))
                assigned_clusters[i][idxa] = idxp
            end
        end
        if new_weight > 0
            push!(subcluster_weight, tmp_subcluster_weight)
            updated_bc = vcat(updated_bc, assigned_points /new_weight)
        end
    end
    if size(updated_bc)[1] == 0
        println(cur_bc)
        error("update_barycenter: This should not happen!")
    end
    updated_bc, mospa(updated_bc, measurements, weights, c), assigned_clusters, subcluster_weight
end

function init_by_measurement(measurements::Vector{Pointcloud},
                             rest...
                            )
    measurements[rand(1:length(measurements))]
end

function init_by_point(measurements::Vector{Pointcloud},
                       assignments::Vector{Vector{I}},
                      ) where {I <: Integer}
    outer_ind = map(x -> any(x .== 0), assignments) |> findall |> rand
    inner_ind = (assignments[outer_ind] .== 0) |> findall |> rand
    measurements[outer_ind][inner_ind:inner_ind, :]
end


function draw_random_points(base_set, lower, upper, rand_scaling)
    n = rand(lower:upper)
    idxs = randperm(size(base_set)[1])
    sample = base_set[idxs[1:n], :]
    sample + rand(size(sample))*rand_scaling
end


function draw_from_remaining_pool(measurements::Vector{Pointcloud}, assignments::Vector{Vector{Integer}})
    outer_ind = map(x-> any(x.==0), assignments) |> findall |> rand
    inner_ind = (assignments[outer_ind] .== 0) |> findall |> rand
    measurements[outer_ind][inner_ind:inner_ind, :], outer_ind, inner_ind
end


function selective_optimal_assignments(pc::Pointcloud, measurements::Vector{Pointcloud}, assignments::Vector{Vector{Integer}})
    opt_assignments = []
    for i=1:length(measurements)
        # println("Index: $(i), pc: $(pc), meas: $(measurements[i]), ass: $(assignments[i]).")
        dmat = restricted_p2dist(pc, measurements[i], assignments[i])
        if all(dmat .== Inf) || isempty(dmat) # The last part is ok?
            push!(opt_assignments, [0])
        else
            push!(opt_assignments, hungarian(dmat)[1])
        end
    end
    opt_assignments
end


function restricted_p2dist(x::Pointcloud, y::Pointcloud, assignments::Vector{Integer})
    dmat = p2dist(x,y)
    dmat[:, assignments .!= 0] .= Inf
    dmat
end


function check_for_remaining_points(not_assigned_points::Vector{Vector{Integer}})
    map(x -> all(x .!= 0), not_assigned_points) |> all
end
