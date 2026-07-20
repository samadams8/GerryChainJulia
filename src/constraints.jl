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
    population_constraint(tolerance::Float64=0.01) -> Function

Returns a validator function `(graph, partition) -> Bool` for population balance
within `tolerance`.
"""
function population_constraint(tolerance::Float64=0.01)
    return (graph, partition) -> within_percent_of_ideal_population(graph, partition, tolerance)
end

"""
    is_contiguous_flip(graph, partition, flip, [visited, [queue]]) -> Bool

Returns `true` if flipping `flip.node` from district `flip.D₁` to `flip.D₂`
leaves `flip.D₁` contiguous. Uses BFS from each remaining neighbor of `flip.node`
in `D₁` to verify they can reach each other without passing through `flip.node`.

`visited` and `queue` are optional pre-allocated scratch buffers; providing them
avoids repeated allocation when calling in a tight loop.
"""
function is_contiguous_flip(
    graph::AbstractGraph,
    partition::AbstractPartition,
    flip::FlipProposal,
    visited::BitVector=BitVector(),
    queue::Vector{Int}=sizehint!(Int[], 64),
)
    nbrs = neighbors(graph)
    asg = assignments(partition)
    node_neighbors = [n for n in nbrs[flip.node] if asg[n] == flip.D₁]
    if isempty(node_neighbors)
        return false
    end
    source_node = pop!(node_neighbors)

    n = num_nodes(graph)
    if length(visited) != n
        resize!(visited, n)
    end

    @inbounds for target_node in node_neighbors
        fill!(visited, false)
        empty!(queue)
        push!(queue, target_node)
        visited[target_node] = true
        found = false
        head = 1
        while head <= length(queue)
            curr_node = queue[head]
            head += 1
            if curr_node == source_node
                found = true
                break
            end
            for neighbor in nbrs[curr_node]
                if (
                    !visited[neighbor] &&
                    asg[neighbor] == flip.D₁ &&
                    neighbor != flip.node
                )
                    visited[neighbor] = true
                    push!(queue, neighbor)
                end
            end
        end
        if !found
            return false
        end
    end
    return true
end
