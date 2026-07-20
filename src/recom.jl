"""
    sample_subgraph(graph::AbstractGraph,
                    partition::AbstractPartition,
                    rng::AbstractRNG)

Randomly sample two adjacent districts D₁ and D₂ and return a tuple
(D₁, D₂, edges, nodes) where D₁ and D₂ are Ints, `edges` and `nodes` are Sets
containing the Int edges and Int nodes of the induced subgraph.
"""
function sample_subgraph(
    graph::AbstractGraph, partition::AbstractPartition, rng::AbstractRNG
)
    D₁, D₂ = sample_adjacent_districts_randomly(partition, rng)

    # Take all their nodes.
    nodes = union(dist_nodes(partition)[D₁], dist_nodes(partition)[D₂])

    # Get a subgraph of these two districts.
    edges = induced_subgraph_edges(graph, collect(nodes))

    return D₁, D₂, edges, BitSet(nodes)
end

"""
    build_mst(graph::AbstractGraph,
              nodes::BitSet,
              edges::BitSet)::Dict{Int, Array{Int, 1}}

Builds a graph as an adjacency list from the `mst_nodes` and `mst_edges`.
"""
function build_mst(
    graph::AbstractGraph, nodes::BitSet, edges::BitSet
)::Dict{Int,Array{Int,1}}
    mst = Dict{Int,Array{Int,1}}()
    for node in nodes
        mst[node] = Array{Int,1}()
    end
    for edge in edges
        add_edge_to_mst!(graph, mst, edge)
    end
    return mst
end

"""
    remove_edge_from_mst!(graph::AbstractGraph,
                          mst::Dict{Int, Array{Int,1}},
                          edge::Int)

Removes an edge from the graph built by `build_mst()`.
"""
function remove_edge_from_mst!(graph::AbstractGraph, mst::Dict{Int,Array{Int,1}}, edge::Int)
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    filter!(e -> e != dsts[edge], mst[srcs[edge]])
    return filter!(e -> e != srcs[edge], mst[dsts[edge]])
end

"""
    add_edge_to_mst!(graph::AbstractGraph,
                     mst::Dict{Int, Array{Int,1}},
                     edge::Int)

    Adds an edge to the graph built by `build_mst()`.
"""
function add_edge_to_mst!(graph::AbstractGraph, mst::Dict{Int,Array{Int,1}}, edge::Int)
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    push!(mst[srcs[edge]], dsts[edge])
    return push!(mst[dsts[edge]], srcs[edge])
end

"""
    traverse_mst(mst::Dict{Int, Array{Int, 1}},
                 start_node::Int,
                 avoid_node::Int,
                 stack::Stack{Int},
                 traversed_nodes::BitSet)::BitSet

Returns the component of the MST `mst` that contains the vertex
`start_node`.

*Arguments:*
    - mst:        mst to traverse
    - start_node: the node to start traversing from
    - avoid_node: the node to avoid and which separates the MST into
                  two components
    - stack:      an empty Stack
    - traversed_nodes: an empty BitSet that is to be populated.

`stack` and `traversed_nodes` are pre-allocated and passed in to
reduce the number of memory allocations and consequently, time taken.
In the course of calling this function multiple times, it is intended that
we pass in the same (empty) objects repeatedly.
"""
function traverse_mst(
    mst::Dict{Int,Array{Int,1}},
    start_node::Int,
    avoid_node::Int,
    stack::Stack{Int},
    traversed_nodes::BitSet,
)::BitSet
    @assert isempty(stack)
    empty!(traversed_nodes)

    push!(stack, start_node)

    while !isempty(stack)
        new_node = pop!(stack)
        push!(traversed_nodes, new_node)

        for neighbor in mst[new_node]
            if !(neighbor in traversed_nodes) && neighbor != avoid_node
                push!(stack, neighbor)
            end
        end
    end
    return traversed_nodes
end

"""
    get_balanced_proposal(graph, mst_edges, mst_nodes, partition, min_pop, max_pop, D₁, D₂)
        -> Union{RecomProposal, DummyProposal}

Attempts a balanced cut on the spanning tree subgraph formed by merging districts
`D₁` and `D₂`. Each MST edge is tried as a cut; the first that produces two
components whose populations both lie in `[min_pop, max_pop]` is returned as
a `RecomProposal`. Returns `DummyProposal` if no balanced cut exists.

This is the `:edge_scan` method — simpler but O(|V|·|E_mst|) per call.
Prefer `get_balanced_proposal_subtree_population` for performance.
"""
function get_balanced_proposal(
    graph::AbstractGraph,
    mst_edges::BitSet,
    mst_nodes::BitSet,
    partition::AbstractPartition,
    min_pop::Int,
    max_pop::Int,
    D₁::Int,
    D₂::Int,
)
    mst = build_mst(graph, mst_nodes, mst_edges)
    pops = dist_populations(partition)
    subgraph_pop = pops[D₁] + pops[D₂]
    srcs = edge_src(graph)
    dsts = edge_dst(graph)

    # Pre-allocated reusable data structures to reduce the number of memory allocations.
    stack = Stack{Int}()
    component_container = BitSet([])

    for edge in mst_edges
        component₁ = traverse_mst(mst, srcs[edge], dsts[edge], stack, component_container)

        population₁ = get_subgraph_population(graph, component₁)
        population₂ = subgraph_pop - population₁

        if population₁ >= min_pop && population₁ <= max_pop && population₂ >= min_pop && population₂ <= max_pop
            component₂ = setdiff(mst_nodes, component₁)
            proposal = RecomProposal(
                D₁, D₂, population₁, population₂, component₁, component₂
            )
            return proposal
        end
    end
    return DummyProposal("Could not find balanced cut.")
end

"""
    get_balanced_proposal_subtree_population(graph, mst_edges, mst_nodes,
`mst_nodes` such that the population is balanced according to `pop_constraint`.

Note: This is the `:subtree_population` method. It is the default cut method.
It evaluates candidate cuts in O(|V|) total time by rooting the MST at an arbitrary
node and computing subtree populations in a single post-order traversal. This is
significantly faster than `get_balanced_proposal`, and avoids per-call allocations
when a reusable `SubtreeCutScratch` is passed to the `scratch` argument.

Iterates `mst_edges` in the same order as `get_balanced_proposal` and scores the
same side (component containing `edge_src`, avoiding `edge_dst`) so the first
accepted cut matches `:edge_scan`.
"""
function get_balanced_proposal_subtree_population(
    graph::AbstractGraph,
    mst_edges::BitSet,
    mst_nodes::BitSet,
    partition::AbstractPartition,
    min_pop::Int,
    max_pop::Int,
    D₁::Int,
    D₂::Int;
    scratch::Union{SubtreeCutScratch,Nothing}=nothing,
)
    isempty(mst_nodes) && return DummyProposal("Could not find balanced cut.")

    pops = dist_populations(partition)
    subgraph_pop = pops[D₁] + pops[D₂]
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    node_pops = populations(graph)

    max_node = maximum(mst_nodes)
    if scratch === nothing
        scratch = SubtreeCutScratch(max_node)
    else
        _ensure_subtree_cut_scratch!(scratch, max_node)
    end

    adj = scratch.adj
    parent = scratch.parent
    subpop = scratch.subpop
    order = scratch.order
    stack = scratch.stack

    @inbounds for e in mst_edges
        u, v = srcs[e], dsts[e]
        push!(adj[u], v)
        push!(adj[v], u)
    end

    root = first(mst_nodes)
    empty!(stack)
    push!(stack, root)
    parent[root] = root  # Mark root as visited with self-parent.
    while !isempty(stack)
        u = pop!(stack)
        push!(order, u)
        @inbounds for v in adj[u]
            if parent[v] == 0 && v != root
                parent[v] = u
                push!(stack, v)
            end
        end
    end
    parent[root] = 0

    @inbounds for i in length(order):-1:1
        u = order[i]
        subpop[u] = node_pops[u]
        for v in adj[u]
            if parent[v] == u
                subpop[u] += subpop[v]
            end
        end
    end

    for edge in mst_edges
        u = srcs[edge]
        v = dsts[edge]
        # Population of the component containing u when edge (u,v) is cut —
        # same side as get_balanced_proposal / traverse_mst(u, avoid=v).
        population₁ = if parent[u] == v
            subpop[u]
        elseif parent[v] == u
            subgraph_pop - subpop[v]
        else
            continue
        end
        population₂ = subgraph_pop - population₁

        if population₁ >= min_pop && population₁ <= max_pop && population₂ >= min_pop && population₂ <= max_pop
            component₁ = _collect_component_dense!(scratch, u, v)
            component₂ = setdiff(mst_nodes, component₁)
            return RecomProposal(D₁, D₂, population₁, population₂, component₁, component₂)
        end
    end
    return DummyProposal("Could not find balanced cut.")
end

"""Collect nodes reachable from `start` without crossing to `avoid`, using scratch buffers."""
function _collect_component_dense!(
    scratch::SubtreeCutScratch, start::Int, avoid::Int
)::BitSet
    seen = scratch.seen
    fill!(seen, false)
    stack = scratch.stack
    empty!(stack)
    push!(stack, start)
    seen[start] = true
    nodes = BitSet()
    adj = scratch.adj
    while !isempty(stack)
        u = pop!(stack)
        push!(nodes, u)
        @inbounds for v in adj[u]
            if v != avoid && !seen[v]
                seen[v] = true
                push!(stack, v)
            end
        end
    end
    return nodes
end

"""
    _spanning_tree(graph, edges, nodes, rng; tree_method, region_surcharges)

Build a spanning tree on the induced subgraph. `:wilson` draws a uniform
spanning tree (penalties / region surcharges are ignored). `:kruskal` uses
random or weighted Kruskal.
"""
function _spanning_tree(
    graph::AbstractGraph,
    edges::Vector{Int},
    nodes::Vector{Int},
    rng::AbstractRNG;
    tree_method::Symbol=:kruskal,
    scratch::Union{MSTScratch,Nothing}=nothing,
)::BitSet
    if tree_method === :wilson
        return wilson_ust(graph, edges, nodes, rng)
    elseif tree_method === :kruskal
        return _kruskal_mst(graph, edges, nodes, rng; scratch=scratch)
    else
        throw(ArgumentError("tree_method must be :kruskal or :wilson, got $(tree_method)"))
    end
end

"""
    _RecomInternalOptions

Internal options struct to consolidate keyword arguments for proposal generation.
"""
struct _RecomInternalOptions
    num_tries::Int
    tree_method::Symbol
    cut_method::Symbol
    n_parallel::Int
end

"""
    _try_valid_proposal(graph, partition, pop_constraint, rng, opts;
                        scratch, cut_scratch)

Attempt one subgraph sample with up to `opts.num_tries` spanning trees.
Returns a `RecomProposal` or `nothing`.
"""
function _try_valid_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition,
    min_pop::Int,
    max_pop::Int,
    rng::AbstractRNG,
    opts::_RecomInternalOptions;
    scratch::Union{MSTScratch,Nothing}=nothing,
    cut_scratch::Union{SubtreeCutScratch,Nothing}=nothing,
)
    D₁, D₂, sg_edges, sg_nodes = sample_subgraph(graph, partition, rng)
    sg_node_list = collect(sg_nodes)

    for _ in 1:opts.num_tries
        mst_edges = _spanning_tree(
            graph,
            sg_edges,
            sg_node_list,
            rng;
            tree_method=opts.tree_method,
            scratch=scratch,
        )
        proposal = if opts.cut_method === :subtree_population
            get_balanced_proposal_subtree_population(
                graph,
                mst_edges,
                sg_nodes,
                partition,
                min_pop,
                max_pop,
                D₁,
                D₂;
                scratch=cut_scratch,
            )
        elseif opts.cut_method === :edge_scan
            get_balanced_proposal(graph, mst_edges, sg_nodes, partition, min_pop, max_pop, D₁, D₂)
        else
            throw(
                ArgumentError(
                    "cut_method must be :subtree_population or :edge_scan, got $(opts.cut_method)",
                ),
            )
        end
        if proposal isa RecomProposal
            return proposal
        end
    end
    return nothing
end

"""
    _first_valid_proposal(results)

Walk the list of task results and return the first one that is a `RecomProposal` (valid proposal).
Because parallel tasks are sorted by index, this resolves tie-breakers in favor of the lowest index task.
"""
function _first_valid_proposal(results)
    for res in results
        res isa RecomProposal && return res
    end
    return nothing
end

function get_valid_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition,
    min_pop::Int,
    max_pop::Int,
    rng::AbstractRNG,
    opts::_RecomInternalOptions,
)
    n_parallel = opts.n_parallel
    n_parallel < 1 && throw(ArgumentError("n_parallel must be ≥ 1"))
    if opts.cut_method !== :subtree_population && opts.cut_method !== :edge_scan
        throw(
            ArgumentError(
                "cut_method must be :subtree_population or :edge_scan, got $(opts.cut_method)",
            ),
        )
    end

    if n_parallel == 1
        mst_scratch = MSTScratch()
        cut_scratch = SubtreeCutScratch()
        while true
            proposal = _try_valid_proposal(
                graph,
                partition,
                min_pop,
                max_pop,
                rng,
                opts;
                scratch=mst_scratch,
                cut_scratch=cut_scratch,
            )
            proposal !== nothing && return proposal
        end
    end

    while true
        seeds = [rand(rng, UInt64) for _ in 1:n_parallel]
        tasks = map(1:n_parallel) do i
            Threads.@spawn begin
                task_rng = Random.MersenneTwister(seeds[i])
                scratch = MSTScratch()
                cut_scratch = SubtreeCutScratch()
                _try_valid_proposal(
                    graph,
                    partition,
                    min_pop,
                    max_pop,
                    task_rng,
                    opts;
                    scratch=scratch,
                    cut_scratch=cut_scratch,
                )
            end
        end
        results = fetch.(tasks)
        best = _first_valid_proposal(results)
        best !== nothing && return best
    end
end

"""
    get_valid_recom_proposal(graph, partition, [rng], [num_tries]; tolerance=0.01,
                             tree_method=:kruskal, n_parallel=1,
                             cut_method=:subtree_population) -> RecomProposal

Returns a population-balanced ReCom proposal within `tolerance` (default 0.01,
meaning each new district must be within ±1% of ideal population).

- `num_tries`: MST samples per subgraph before resampling a new subgraph.
- `tree_method`: `:kruskal` (default) or `:wilson` for uniform spanning trees.
- `n_parallel`: number of concurrent proposal attempts (`1` = serial).
- `cut_method`: `:subtree_population` (default, O(N)) or `:edge_scan`.
"""
function get_valid_recom_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition,
    rng::AbstractRNG=Random.default_rng(),
    num_tries::Int=3;
    tolerance::Float64=0.01,
    tree_method::Symbol=:kruskal,
    n_parallel::Int=1,
    cut_method::Symbol=:subtree_population,
)
    ideal_pop = total_pop(graph) / num_dists(partition)
    min_pop = Int(ceil((1 - tolerance) * ideal_pop))
    max_pop = Int(floor((1 + tolerance) * ideal_pop))
    opts = _RecomInternalOptions(num_tries, tree_method, cut_method, n_parallel)
    return get_valid_proposal(graph, partition, min_pop, max_pop, rng, opts)
end

"""
    recom_proposal(graph, partition; tolerance=0.01, rng=Random.default_rng()) -> Partition

High-level ReCom proposal function. Generates a population-balanced `RecomProposal`
and returns a new `Partition` with the proposal applied.
Intended for use as the `proposal` argument to `MarkovChain`.
"""
function recom_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition;
    tolerance::Float64=0.01,
    rng::AbstractRNG=Random.default_rng(),
)
    prop = get_valid_recom_proposal(graph, partition, rng; tolerance=tolerance)
    p_next = clone_for_update(partition)
    return update_partition!(p_next, graph, prop)
end

"""
    update_partition!(partition::Partition,
                      graph::AbstractGraph,
                      proposal::RecomProposal,
                      copy_parent::Bool=false)

Updates the `Partition` with the `RecomProposal`.

When `copy_parent` is true, a field-wise snapshot of the current partition
is stored in `partition.parent` (no recursive `deepcopy`). Prefer
`clone_for_update` + `update_partition!(..., false)` when returning a new state.
"""
function update_partition!(
    partition::Partition,
    graph::AbstractGraph,
    proposal::RecomProposal,
    copy_parent::Bool=false,
)
    if copy_parent
        partition.parent = _copy_partition_fields(partition; parent=nothing)
    end

    partition.dist_populations[proposal.D₁] = proposal.D₁_pop
    partition.dist_populations[proposal.D₂] = proposal.D₂_pop

    for node in proposal.D₁_nodes
        partition.assignments[node] = proposal.D₁
    end
    for node in proposal.D₂_nodes
        partition.assignments[node] = proposal.D₂
    end

    # Replace (do not mutate shared BitSets from a CoW clone).
    partition.dist_nodes[proposal.D₁] = proposal.D₁_nodes
    partition.dist_nodes[proposal.D₂] = proposal.D₂_nodes

    return update_partition_adjacency(partition, graph)
end
