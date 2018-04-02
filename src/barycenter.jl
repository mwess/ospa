function ospa_barycenter(measurements::Vector{Pointcloud}, 
                         weights::Vector{R}, 
                         c::S
                         ;
                         niter::Integer=100, 
                         eps::Float64=0.05, 
                         lower_bound::Integer=1, 
                         upper_bound::Integer=0
                        ) where {R <: Real, S <: Real}
    if length(measurements) == 0
        error("No measurements given!")
    elseif length(measurements) != length(weights)
        error("Number of measurements $(length(measurements)) and weights $(length(weights)) does not match!")
    elseif c <= 0
        error("Parameter c ($(c)) must be greater than 0.")
    elseif upper_bound < 0
        error("Upper bound $(upper_bound) needs to be positive (or zero for dynamic computation)!")
    end
    Barycenters = Vector{Pointcloud}()
    costs = Vector{Float64}()
    #What is the upper bound on the size of the optimal barycenter
    if upper_bound == 0
        upper_bound = maximum(map(x -> size(x)[1], measurements))
    end
    for i=lower_bound:upper_bound
        barycenter, cost = compute_fixed_size_barycenter(measurements, weights, c, i, niter, eps)
        push!(Barycenters, barycenter)
        push!(costs, cost)
    end
    index = indmin(costs)
    Barycenters[index], costs[index]
end

function compute_fixed_size_barycenter(measurements::Vector{Pointcloud}, 
                                       weights::Vector{R}, 
                                       c::S, 
                                       sz::Integer, 
                                       niter::Integer
                                       ;
                                       eps::Float64=0.05
                                      ) where {R <: Real, S <: Real}
    lowest_cost = Inf
    lowest_barycenter = []
    for i=1:niter
        old_cost = Inf
        converged = false
        current_barycenter = init_bc(measurements,sz)
        while !converged
            current_barycenter, current_cost = update_barycenter(current_barycenter, measurements, weights, c, sz)
            if old_cost - current_cost < eps
                converged = true
                if current_cost < lowest_cost
                    lowest_cost = current_cost
                    lowest_barycenter = current_barycenter
                end
            end
            old_cost = current_cost
        end
    end
    lowest_barycenter, lowest_cost
end

function init_bc(measurements::Vector{Pointcloud}, 
                 sz::Integer)
    #For now just randomly select measurements
    possible_inds = find(x -> size(x)[1] >= sz, measurements)
    bc = measurements[rand(possible_inds)]
    bc[randperm(size(bc)[1]),:][1:sz,:]
end

function update_barycenter(cur_bc::Pointcloud, 
                           measurements::Vector{Pointcloud}, 
                           weights::Vector{R}, 
                           c::S, 
                           sz::Integer
                          ) where {R <: Real, S <: Real}
    updated_bc = Matrix{eltype(cur_bc)}(0,size(cur_bc)[2])
    assignments = optimal_assignments(cur_bc, measurements)
    for idxp = 1:size(cur_bc)[1]
        assigned_points = zeros(eltype(cur_bc), 1,size(cur_bc)[2])
        new_weight = 0
        for i=1:length(assignments)
            idxa = assignments[i][idxp]
            if idxa != 0 && p2dist(cur_bc[idxp:idxp,:], measurements[i][idxa:idxa,:])[1] <= c
                assigned_points += measurements[i][idxa:idxa,:]*weights[idxa]/max(sz, length(assignments[i])) 
                new_weight += weights[idxa]/max(sz, length(assignments[i]))
            end
        end
        if new_weight > 0
            updated_bc = vcat(updated_bc, assigned_points /new_weight)
        end
    end
    if size(updated_bc)[1] == 0
        println(cur_bc)
        error("update_barycenter: This should not happen!")
    end
    updated_bc, mospa(updated_bc, measurements, weights, c)
end

function mospa(bc::Pointcloud, 
               measurements::Vector{Pointcloud}, 
               weights::Vector{R}, 
               c::S
              ) where {S <: Real, R <: Real}
    sum(map(x -> weights[x]*ospa_dist(bc, measurements[x], c), 1:length(measurements)))
end

