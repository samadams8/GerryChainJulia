"""
    _is_contiguous_flip(graph, partition, flip, [visited, [queue]]) -> Bool

Returns `true` if flipping `flip.node` from district `flip.D₁` to `flip.D₂`
leaves `flip.D₁` contiguous. Uses BFS from each remaining neighbor of `flip.node`
in `D₁` to verify they can reach each other without passing through `flip.node`.
"""
function _is_contiguous_flip(
    graph::AbstractGraph,
    partition::AbstractPartition,
    flip::FlipPayload,
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

"""
    propose_random_flip(graph::AbstractGraph,
                        partition::AbstractPartition,
                        rng::AbstractRNG=Random.default_rng()) -> FlipPayload

Proposes a random boundary flip: picks a random cut edge, randomly assigns
one endpoint to the other's district, and returns the resulting `FlipPayload`.
"""
function propose_random_flip(
    graph::AbstractGraph,
    partition::AbstractPartition,
    rng::AbstractRNG=Random.default_rng(),
)
    if num_cut_edges(partition) == 0
        throw(ArgumentError("No cut edges in the districting plan"))
    end
    # Select a random cut edge.
    cut_edge_idx = rand(rng, 1:num_cut_edges(partition))
    cut_edge_tracker = 0
    edge_idx = 0
    cut = cut_edges(partition)
    # Iterate through array of bools indicating cut edge, stop at the
    # randomly chosen index-th edge.
    for i in 1:num_edges(graph)
        cut_edge_tracker += cut[i]
        if cut_edge_tracker == cut_edge_idx
            edge_idx = i
            break
        end
    end
    # Randomly choose which of the nodes from the edge get flipped.
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    edge = (srcs[edge_idx], dsts[edge_idx])
    index = rand(rng, (0, 1))
    flipped_node, other_node = edge[index + 1], edge[2 - index]
    node_pop = populations(graph)[flipped_node]
    asg = assignments(partition)
    pops = dist_populations(partition)
    nodes = dist_nodes(partition)
    # Old district.
    D₁ = asg[flipped_node]
    D₁_pop = pops[D₁] - node_pop
    D₁_n = setdiff(nodes[D₁], flipped_node)
    # New district.
    D₂ = asg[other_node]
    D₂_pop = pops[D₂] + node_pop
    D₂_n = union(nodes[D₂], flipped_node)
    return FlipPayload(flipped_node, D₁, D₂, D₁_pop, D₂_pop, D₁_n, D₂_n)
end

"""
    is_valid(graph, partition, min_pop, max_pop, proposal, [visited, queue]) -> Bool

Returns `true` iff the `FlipPayload` satisfies both population balance
(`min_pop ≤ each district ≤ max_pop`) and contiguity of the source district.
Reusable `visited` and `queue` scratch buffers may be provided to avoid allocations.
"""
function is_valid(
    graph::AbstractGraph,
    partition::AbstractPartition,
    min_pop::Int,
    max_pop::Int,
    proposal::FlipPayload,
    visited::BitVector=BitVector(),
    queue::Vector{Int}=sizehint!(Int[], 64),
)
    pop_ok = proposal.D₁_pop >= min_pop && proposal.D₁_pop <= max_pop &&
             proposal.D₂_pop >= min_pop && proposal.D₂_pop <= max_pop
    return pop_ok && _is_contiguous_flip(graph, partition, proposal, visited, queue)
end

"""
    get_valid_proposal(graph, partition, min_pop, max_pop, [rng]) -> FlipPayload

Samples random flips until one satisfies `[min_pop, max_pop]` population bounds
and district contiguity. Scratch buffers are allocated once and reused internally.
"""
function get_valid_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition,
    min_pop::Int,
    max_pop::Int,
    rng::AbstractRNG=Random.default_rng(),
)
    visited = BitVector()
    queue = sizehint!(Int[], 64)
    proposal = propose_random_flip(graph, partition, rng)
    while !is_valid(graph, partition, min_pop, max_pop, proposal, visited, queue)
        proposal = propose_random_flip(graph, partition, rng)
    end
    return proposal
end

"""
    update_partition!(partition::Partition,
                      graph::AbstractGraph,
                      proposal::FlipPayload,
                      copy_parent::Bool=false) -> Partition

Applies a `FlipPayload` to `partition` in-place. When `copy_parent` is `true`,
a shallow snapshot of the current partition is stored in `partition.parent` before
updating.
"""
function update_partition!(
    partition::Partition,
    graph::AbstractGraph,
    proposal::FlipPayload,
    copy_parent::Bool=false,
)
    if copy_parent
        partition.parent = _copy_partition_fields(partition; parent=nothing)
    end

    # Update district population counts.
    partition.dist_populations[proposal.D₁] = proposal.D₁_pop
    partition.dist_populations[proposal.D₂] = proposal.D₂_pop

    # Relabel node with new district.
    partition.assignments[proposal.node] = proposal.D₂

    # Replace (do not mutate shared CoW BitSets).
    partition.dist_nodes[proposal.D₁] = proposal.D₁_nodes
    partition.dist_nodes[proposal.D₂] = proposal.D₂_nodes

    return update_partition_adjacency(partition, graph)
end

"""
    PopulationFlipConfiguration{RNG<:AbstractRNG} <: AbstractProposalConfiguration

Configuration for running a population-balancing boundary node flip proposal.

# Fields
- `ideal_pop::Float64`: Target population for each district.
- `pop_key::String`: Population updater/attribute key (for compatibility).
- `rng::RNG`: Random number generator.
"""
struct PopulationFlipConfiguration{RNG<:AbstractRNG} <: AbstractProposalConfiguration
    ideal_pop::Float64
    pop_key::String
    rng::RNG
end

function PopulationFlipConfiguration(
    ideal_pop::Float64,
    pop_key::String = "population";
    rng::AbstractRNG = Random.default_rng()
)
    return PopulationFlipConfiguration(ideal_pop, pop_key, rng)
end

function propose(
    graph::AbstractGraph,
    partition::Partition,
    config::PopulationFlipConfiguration,
)
    cut = cut_edges(partition)
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    asg = assignments(partition)
    pops = dist_populations(partition)

    candidate_edges = Tuple{Int,Int}[]
    for i in 1:num_edges(graph)
        if cut[i] == 1
            u = srcs[i]
            v = dsts[i]
            D_u = asg[u]
            D_v = asg[v]
            pop_u = pops[D_u]
            pop_v = pops[D_v]

            balanced_u = abs(pop_u - config.ideal_pop) < 1.0
            balanced_v = abs(pop_v - config.ideal_pop) < 1.0

            if !(balanced_u && balanced_v)
                push!(candidate_edges, (u, v))
            end
        end
    end

    if isempty(candidate_edges)
        return partition
    end

    edge = rand(config.rng, candidate_edges)
    index = rand(config.rng, (0, 1))
    flipped_node, other_node = edge[index + 1], edge[2 - index]

    node_pop = populations(graph)[flipped_node]
    D₁ = asg[flipped_node]
    D₁_pop = pops[D₁] - node_pop
    D₁_n = setdiff(dist_nodes(partition)[D₁], flipped_node)

    D₂ = asg[other_node]
    D₂_pop = pops[D₂] + node_pop
    D₂_n = union(dist_nodes(partition)[D₂], flipped_node)

    prop = FlipPayload(flipped_node, D₁, D₂, D₁_pop, D₂_pop, D₁_n, D₂_n)
    p_next = clone_for_update(partition)
    return update_partition!(p_next, graph, prop)
end
