"""
    propose_random_flip(graph::AbstractGraph,
                        partition::AbstractPartition,
                        rng::AbstractRNG=Random.default_rng()) -> FlipProposal

Proposes a random boundary flip: picks a random cut edge, randomly assigns
one endpoint to the other's district, and returns the resulting `FlipProposal`.
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
    return FlipProposal(flipped_node, D₁, D₂, D₁_pop, D₂_pop, D₁_n, D₂_n)
end

"""
    is_valid(graph, partition, min_pop, max_pop, proposal, [visited, queue]) -> Bool

Returns `true` iff the `FlipProposal` satisfies both population balance
(`min_pop ≤ each district ≤ max_pop`) and contiguity of the source district.
Reusable `visited` and `queue` scratch buffers may be provided to avoid allocations.
"""
function is_valid(
    graph::AbstractGraph,
    partition::AbstractPartition,
    min_pop::Int,
    max_pop::Int,
    proposal::FlipProposal,
    visited::BitVector=BitVector(),
    queue::Vector{Int}=sizehint!(Int[], 64),
)
    pop_ok = proposal.D₁_pop >= min_pop && proposal.D₁_pop <= max_pop &&
             proposal.D₂_pop >= min_pop && proposal.D₂_pop <= max_pop
    return pop_ok && is_contiguous_flip(graph, partition, proposal, visited, queue)
end

"""
    get_valid_proposal(graph, partition, min_pop, max_pop, [rng]) -> FlipProposal

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
    get_valid_flip_proposal(graph, partition, [rng]; tolerance=0.01) -> FlipProposal

Computes ideal population bounds from `tolerance` and returns a valid `FlipProposal`.
`tolerance=0.01` means each district must be within ±1% of the ideal population.
"""
function get_valid_flip_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition,
    rng::AbstractRNG=Random.default_rng();
    tolerance::Float64=0.01,
)
    ideal_pop = total_pop(graph) / num_dists(partition)
    min_pop = Int(ceil((1 - tolerance) * ideal_pop))
    max_pop = Int(floor((1 + tolerance) * ideal_pop))
    return get_valid_proposal(graph, partition, min_pop, max_pop, rng)
end

"""
    flip_proposal(graph, partition; tolerance=0.01, rng=Random.default_rng()) -> Partition

High-level flip proposal function. Generates a population-balanced, contiguous
`FlipProposal` and returns a new `Partition` with the flip applied.
Intended for use as the `proposal` argument to `MarkovChain`.
"""
function flip_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition;
    tolerance::Float64=0.01,
    rng::AbstractRNG=Random.default_rng(),
)
    prop = get_valid_flip_proposal(graph, partition, rng; tolerance=tolerance)
    p_next = clone_for_update(partition)
    return update_partition!(p_next, graph, prop)
end

"""
    update_partition!(partition::Partition,
                      graph::AbstractGraph,
                      proposal::FlipProposal,
                      copy_parent::Bool=false) -> Partition

Applies a `FlipProposal` to `partition` in-place. When `copy_parent` is `true`,
a shallow snapshot of the current partition is stored in `partition.parent` before
updating (useful for acceptance functions that need to compare to the previous state).
"""
function update_partition!(
    partition::Partition,
    graph::AbstractGraph,
    proposal::FlipProposal,
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
