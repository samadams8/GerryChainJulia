"""
    AbstractConstraint

Abstract supertype for all constraint definitions in GerryChainJulia.
Concrete subtypes must implement `satisfies_constraint(c::AbstractConstraint, graph, partition) -> Bool`.
"""
abstract type AbstractConstraint end

"""
    satisfies_constraint(constraint::AbstractConstraint,
                         graph::AbstractGraph,
                         partition::AbstractPartition) -> Bool

Return `true` iff `partition` satisfies `constraint` on `graph`.
Must be implemented by all concrete subtypes of `AbstractConstraint`.
"""
function satisfies_constraint end

"""
    within_population_bounds(partition::AbstractPartition, min_pop::Int, max_pop::Int) -> Bool

Returns `true` if all district populations in `partition` are between `min_pop` and
`max_pop` (inclusive).
"""
function within_population_bounds(partition::AbstractPartition, min_pop::Int, max_pop::Int)::Bool
    pops = dist_populations(partition)
    for p in pops
        (p < min_pop || p > max_pop) && return false
    end
    return true
end

"""
    within_percent_of_ideal_population(graph::AbstractGraph, partition::AbstractPartition, tolerance::Float64=0.01) -> Bool

Returns `true` if all district populations in `partition` are within `tolerance`
percentage of ideal population.
"""
function within_percent_of_ideal_population(
    graph::AbstractGraph, partition::AbstractPartition, tolerance::Float64=0.01
)::Bool
    ideal_pop = total_pop(graph) / num_dists(partition)
    min_pop = Int(ceil((1 - tolerance) * ideal_pop))
    max_pop = Int(floor((1 + tolerance) * ideal_pop))
    return within_population_bounds(partition, min_pop, max_pop)
end

"""
    PopulationConstraint <: AbstractConstraint

Constraint requiring all district populations in a partition to be within `tolerance`
percentage of ideal population.
"""
struct PopulationConstraint <: AbstractConstraint
    tolerance::Float64
end

PopulationConstraint() = PopulationConstraint(0.01)

function satisfies_constraint(
    c::PopulationConstraint, graph::AbstractGraph, partition::AbstractPartition
)::Bool
    return within_percent_of_ideal_population(graph, partition, c.tolerance)
end

"""
    ContiguityConstraint <: AbstractConstraint

Constraint requiring all districts in a partition to be contiguous subgraphs.
"""
struct ContiguityConstraint <: AbstractConstraint end

function _is_district_contiguous(graph::AbstractGraph, partition::AbstractPartition, d::Int)::Bool
    d_nodes = dist_nodes(partition)[d]
    isempty(d_nodes) && return true

    start_node = first(d_nodes)
    nbrs = neighbors(graph)
    asg = assignments(partition)

    visited_count = 0
    queue = [start_node]
    visited = Set{Int}([start_node])

    while !isempty(queue)
        curr = popfirst!(queue)
        visited_count += 1
        for nbr in nbrs[curr]
            if asg[nbr] == d && !(nbr in visited)
                push!(visited, nbr)
                push!(queue, nbr)
            end
        end
    end

    return visited_count == length(d_nodes)
end

function satisfies_constraint(
    c::ContiguityConstraint, graph::AbstractGraph, partition::AbstractPartition
)::Bool
    for d in 1:num_dists(partition)
        !_is_district_contiguous(graph, partition, d) && return false
    end
    return true
end
