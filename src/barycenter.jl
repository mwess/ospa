
#=
Mean initialization:
    - Compute bins, then sample from bins and compute barycenter from samples.
    - Necessary parameters:
        Pointcloud
        c
        n_samples
        lower_bound
        upper_bound
=#


function ospa_barycenter(measurements::Vector{Pointcloud},
                         c::S,
                         weights::Union{Vector{R}, Missing}=missing,
                         barycenter_init::Symbol=:binned,
                         niter::Integer=100,
                         eps::Float64=0.05,
                         min_card::Integer=1,
                         max_card::Integer=0,
                         n_samples::Integer=10
                        ) where {R <: Real, S <: Real}
    if length(measurements) == 0
        error("No measurements given!")
    elseif c <= 0
        error("Parameter c ($(c)) must be greater than 0.")
    elseif max_card < 0
        error("Maximum cardinality $(max_card) needs to be positive.")
    elseif min_card < 0
        error("Minimum cardinality $(min_card) needs to be positive.")
    end
    if ismissing(weights)
        weights = fill(1/length(measurements), length(measurements))
    end
    if barycenter_init==:binned
        return binned_barycenter_initialization(measurements, weights, c, min_card, max_card, eps, n_samples )
    elseif barycenter_init==:iterative
        return recursive_barycenter_initialization(measurements, weights, c, eps)
    elseif barycenter_init==:standard
        return standard_barycenter_initilization(measurements, weights, c, niter, eps, min_card, max_card)
    else
        error("Invalid initialization method!")
    end
end

function standard_barycenter_initilization(measurements, weights, c, niter, eps, min_card, max_card)
    Barycenters = Vector{Pointcloud}()
    costs = Vector{Float64}()
    #What is the upper bound on the size of the optimal barycenter
    if max_card == 0
        max_card = maximum(map(x -> size(x)[1], measurements))
    end
    max_card = min(max_card, maximum(map(x -> size(x)[1], measurements)))
    for i=min_card:max_card
        barycenter, cost = compute_fixed_size_barycenter(measurements, weights, c, i, niter, eps)
        push!(Barycenters, barycenter)
        push!(costs, cost)
    end
    (_, index) = findmin(costs)
    Barycenters[index], costs[index]
end


function ospa_barycenter(measurements::Vector{Pointcloud},
                         weights::Vector{R},
                         c::S,
                         current_barycenter::Pointcloud,
                         eps::Float64=0.05
                        ) where {R <: Real, S <: Real}
   converged = false
   old_cost = Inf
   current_cost = Inf
   sz = size(current_barycenter)[1]
   while !converged
        current_barycenter, current_cost = update_barycenter(current_barycenter, measurements, weights, c, sz)
        if old_cost - current_cost < eps
            converged = true
        end
        old_cost = current_cost
   end
   current_barycenter, current_cost
end

function compute_fixed_size_barycenter(measurements::Vector{Pointcloud},
                                       weights::Vector{R},
                                       c::S,
                                       sz::Integer,
                                       niter::Integer,
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
    possible_inds = findall(x -> size(x)[1] >= sz, measurements)
    bc = measurements[rand(possible_inds)]
    bc[randperm(size(bc)[1]),:][1:sz,:]
end

function update_barycenter(cur_bc::Pointcloud,
                           measurements::Vector{Pointcloud},
                           weights::Vector{R},
                           c::S,
                           sz::Integer
                          ) where {R <: Real, S <: Real}
    #updated_bc = Matrix{eltype(cur_bc)}(0,size(cur_bc)[2])
    updated_bc = Matrix{eltype(cur_bc)}(fill(0, 0, size(cur_bc)[2]))
    assignments = optimal_assignments(cur_bc, measurements)
    for idxp = 1:size(cur_bc)[1]
        assigned_points = zeros(eltype(cur_bc), 1,size(cur_bc)[2])
        new_weight = 0
        for i=1:length(assignments)
            idxa = assignments[i][idxp]
            if idxa != 0 && p2dist(cur_bc[idxp:idxp,:], measurements[i][idxa:idxa,:])[1] <= c
                assigned_points += measurements[i][idxa:idxa,:]*weights[i]/max(sz, length(assignments[i]))
                new_weight += weights[i]/max(sz, length(assignments[i]))
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
