"""
    sample_subgraph(graph::AbstractGraph,
                    partition::AbstractPartition,
                    rng::AbstractRNG)

Randomly sample two adjacent districts D₁ and D₂ and return a tuple
(D₁, D₂, edges, nodes) where D₁ and D₂ are Ints, `edges` and `nodes` are Sets
containing the Int edges and Int nodes of the induced subgraph.
"""
function sample_subgraph(
    graph::AbstractGraph,
    partition::AbstractPartition,
    rng::AbstractRNG,
)
    D₁, D₂ = sample_adjacent_districts_randomly(partition, rng)

    # take all their nodes
    nodes = union(dist_nodes(partition)[D₁], dist_nodes(partition)[D₂])

    # get a subgraph of these two districts
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
    graph::AbstractGraph,
    nodes::BitSet,
    edges::BitSet,
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
function remove_edge_from_mst!(
    graph::AbstractGraph,
    mst::Dict{Int,Array{Int,1}},
    edge::Int,
)
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    filter!(e -> e != dsts[edge], mst[srcs[edge]])
    filter!(e -> e != srcs[edge], mst[dsts[edge]])
end

"""
    add_edge_to_mst!(graph::AbstractGraph,
                     mst::Dict{Int, Array{Int,1}},
                     edge::Int)

    Adds an edge to the graph built by `build_mst()`.
"""
function add_edge_to_mst!(
    graph::AbstractGraph,
    mst::Dict{Int,Array{Int,1}},
    edge::Int,
)
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    push!(mst[srcs[edge]], dsts[edge])
    push!(mst[dsts[edge]], srcs[edge])
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
    - avoid_node: the node to avoid adn which seperates the mst into
                  two components
    - stack:      an empty Stack
    - traversed_nodes: an empty BitSet that is to be populated.

`stack` and `traversed_nodes` are are pre-allocated and passed in to
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
    get_balanced_proposal(graph::AbstractGraph,
                          mst_edges::BitSet,
                          mst_nodes::BitSet,
                          partition::AbstractPartition,
                          pop_constraint::PopulationConstraint,
                          D₁::Int,
                          D₂::Int)

Tries to find a balanced cut on the subgraph induced by `mst_edges` and
`mst_nodes` such that the population is balanced according to
`pop_constraint`.
This subgraph was formed by the combination of districts `D₁` and `D₂`.
"""
function get_balanced_proposal(
    graph::AbstractGraph,
    mst_edges::BitSet,
    mst_nodes::BitSet,
    partition::AbstractPartition,
    pop_constraint::PopulationConstraint,
    D₁::Int,
    D₂::Int,
)
    mst = build_mst(graph, mst_nodes, mst_edges)
    pops = dist_populations(partition)
    subgraph_pop = pops[D₁] + pops[D₂]
    srcs = edge_src(graph)
    dsts = edge_dst(graph)

    # pre-allocated reusable data structures to reduce number of memory allocations
    stack = Stack{Int}()
    component_container = BitSet([])

    for edge in mst_edges
        component₁ = traverse_mst(
            mst,
            srcs[edge],
            dsts[edge],
            stack,
            component_container,
        )

        population₁ = get_subgraph_population(graph, component₁)
        population₂ = subgraph_pop - population₁

        if satisfy_constraint(pop_constraint, population₁, population₂)
            component₂ = setdiff(mst_nodes, component₁)
            proposal =
                RecomProposal(D₁, D₂, population₁, population₂, component₁, component₂)
            return proposal
        end
    end
    return DummyProposal("Could not find balanced cut.")
end

"""
    get_balanced_proposal_subtree_population(graph, mst_edges, mst_nodes,
        partition, pop_constraint, D₁, D₂)

Same contract as `get_balanced_proposal`, but evaluates candidate cuts using
rooted-tree subtree populations in O(|V|) instead of walking each side per edge.

Iterates `mst_edges` in the same order as `get_balanced_proposal` and scores the
same side (component containing `edge_src`, avoiding `edge_dst`) so the first
accepted cut matches `:edge_scan`.
"""
function get_balanced_proposal_subtree_population(
    graph::AbstractGraph,
    mst_edges::BitSet,
    mst_nodes::BitSet,
    partition::AbstractPartition,
    pop_constraint::PopulationConstraint,
    D₁::Int,
    D₂::Int,
)
    isempty(mst_nodes) && return DummyProposal("Could not find balanced cut.")

    pops = dist_populations(partition)
    subgraph_pop = pops[D₁] + pops[D₂]
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    node_pops = populations(graph)

    max_node = maximum(mst_nodes)
    adj = [Int[] for _ = 1:max_node]
    @inbounds for e in mst_edges
        u, v = srcs[e], dsts[e]
        push!(adj[u], v)
        push!(adj[v], u)
    end

    parent = zeros(Int, max_node)
    subpop = zeros(Int, max_node)
    order = Int[]
    sizehint!(order, length(mst_nodes))

    root = first(mst_nodes)
    stack = Int[root]
    parent[root] = root  # mark root as visited with self-parent
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

    @inbounds for i = length(order):-1:1
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
            # Should not happen on a tree; fall back to walking
            continue
        end
        population₂ = subgraph_pop - population₁

        if satisfy_constraint(pop_constraint, population₁, population₂)
            component₁ = _collect_component_dense(adj, u, v, max_node)
            component₂ = setdiff(mst_nodes, component₁)
            return RecomProposal(D₁, D₂, population₁, population₂, component₁, component₂)
        end
    end
    return DummyProposal("Could not find balanced cut.")
end

"""Collect nodes reachable from `start` without crossing to `avoid`, using dense adj."""
function _collect_component_dense(
    adj::Vector{Vector{Int}},
    start::Int,
    avoid::Int,
    max_node::Int,
)::BitSet
    seen = falses(max_node)
    stack = Int[start]
    seen[start] = true
    nodes = BitSet()
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
    edges::Array{Int,1},
    nodes::Array{Int,1},
    rng::AbstractRNG;
    tree_method::Symbol = :kruskal,
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
    scratch::Union{MSTScratch,Nothing} = nothing,
)::BitSet
    if tree_method === :wilson
        return wilson_ust(graph, edges, nodes, rng)
    elseif tree_method !== :kruskal
        throw(ArgumentError("tree_method must be :kruskal or :wilson, got $(tree_method)"))
    end

    use_weighted =
        !isempty(region_surcharges) || any(!iszero, edge_penalties(graph))
    if scratch === nothing
        return if use_weighted
            weighted_kruskal_mst(
                graph,
                edges,
                nodes,
                rng;
                region_surcharges = region_surcharges,
            )
        else
            random_kruskal_mst(graph, edges, nodes, rng)
        end
    end

    if use_weighted
        build_mst_weights!(
            scratch,
            graph,
            edges,
            rng;
            region_surcharges = region_surcharges,
        )
    else
        _ensure_mst_scratch!(
            scratch,
            length(edges),
            isempty(nodes) ? 0 : maximum(nodes),
        )
        @inbounds for i = 1:length(edges)
            scratch.weights[i] = rand(rng)
        end
    end
    return kruskal_mst!(scratch, graph, edges, nodes, scratch.weights)
end

"""
    _try_valid_proposal(graph, partition, pop_constraint, rng, num_tries;
                        tree_method, region_surcharges, scratch, cut_method)

Attempt one subgraph sample with up to `num_tries` spanning trees.
Returns a `RecomProposal` or `nothing`.
"""
function _try_valid_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition,
    pop_constraint::PopulationConstraint,
    rng::AbstractRNG,
    num_tries::Int;
    tree_method::Symbol = :kruskal,
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
    scratch::Union{MSTScratch,Nothing} = nothing,
    cut_method::Symbol = :subtree_population,
)
    D₁, D₂, sg_edges, sg_nodes = sample_subgraph(graph, partition, rng)
    sg_node_list = collect(sg_nodes)

    for _ = 1:num_tries
        mst_edges = _spanning_tree(
            graph,
            sg_edges,
            sg_node_list,
            rng;
            tree_method = tree_method,
            region_surcharges = region_surcharges,
            scratch = scratch,
        )
        proposal = if cut_method === :subtree_population
            get_balanced_proposal_subtree_population(
                graph,
                mst_edges,
                sg_nodes,
                partition,
                pop_constraint,
                D₁,
                D₂,
            )
        elseif cut_method === :edge_scan
            get_balanced_proposal(
                graph,
                mst_edges,
                sg_nodes,
                partition,
                pop_constraint,
                D₁,
                D₂,
            )
        else
            throw(ArgumentError(
                "cut_method must be :subtree_population or :edge_scan, got $(cut_method)",
            ))
        end
        if proposal isa RecomProposal
            return proposal
        end
    end
    return nothing
end

"""
    get_valid_proposal(graph::AbstractGraph,
                       partition::AbstractPartition,
                       pop_constraint::PopulationConstraint,
                       rng::AbstractRNG,
                       num_tries::Int=3;
                       region_surcharges=Dict{String,Float64}(),
                       tree_method::Symbol=:kruskal,
                       n_parallel::Int=1,
                       cut_method::Symbol=:subtree_population)

*Returns* a population balanced proposal.

*Arguments:*
    - graph:          AbstractGraph
    - partition:      AbstractPartition
    - pop_constraint: PopulationConstraint to adhere to
    - num_tries:      num times to try getting a balanced cut from a subgraph
                      before giving up
    - rng:            A random number generator that implements the
                      [AbstractRNG type](https://docs.julialang.org/en/v1/stdlib/Random/#Random.AbstractRNG)
                      (e.g. `Random.default_rng()` or `MersenneTwister(1234)`)
    - region_surcharges: Optional per-region-column surcharges added to MST
                      weights when an edge crosses a region boundary
                      (Kruskal only; ignored for `:wilson`)
    - tree_method:    `:kruskal` (default) or `:wilson` (uniform spanning tree)
    - n_parallel:     number of concurrent proposal attempts (`1` = serial).
                      With `n_parallel > 1`, tasks use independent RNGs; the
                      lowest task index among successes in a batch wins.
    - cut_method:     `:subtree_population` (default; O(N) cut search) or
                      `:edge_scan` (walk-and-sum each side per edge)
"""
function get_valid_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition,
    pop_constraint::PopulationConstraint,
    rng::AbstractRNG,
    num_tries::Int = 3;
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
    tree_method::Symbol = :kruskal,
    n_parallel::Int = 1,
    cut_method::Symbol = :subtree_population,
)
    n_parallel < 1 && throw(ArgumentError("n_parallel must be ≥ 1"))
    if cut_method !== :subtree_population && cut_method !== :edge_scan
        throw(ArgumentError(
            "cut_method must be :subtree_population or :edge_scan, got $(cut_method)",
        ))
    end

    if n_parallel == 1
        while true
            proposal = _try_valid_proposal(
                graph,
                partition,
                pop_constraint,
                rng,
                num_tries;
                tree_method = tree_method,
                region_surcharges = region_surcharges,
                cut_method = cut_method,
            )
            proposal !== nothing && return proposal
        end
    end

    while true
        seeds = [rand(rng, UInt64) for _ = 1:n_parallel]
        tasks = map(1:n_parallel) do i
            Threads.@spawn begin
                task_rng = Random.MersenneTwister(seeds[i])
                scratch = MSTScratch()
                _try_valid_proposal(
                    graph,
                    partition,
                    pop_constraint,
                    task_rng,
                    num_tries;
                    tree_method = tree_method,
                    region_surcharges = region_surcharges,
                    scratch = scratch,
                    cut_method = cut_method,
                )
            end
        end
        results = fetch.(tasks)
        best = nothing
        best_idx = typemax(Int)
        for (i, result) in enumerate(results)
            if result isa RecomProposal && i < best_idx
                best = result
                best_idx = i
            end
        end
        best !== nothing && return best
    end
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
    copy_parent::Bool = false,
)
    if copy_parent
        partition.parent = _copy_partition_fields(partition; parent = nothing)
    end

    partition.dist_populations[proposal.D₁] = proposal.D₁_pop
    partition.dist_populations[proposal.D₂] = proposal.D₂_pop

    for node in proposal.D₁_nodes
        partition.assignments[node] = proposal.D₁
    end
    for node in proposal.D₂_nodes
        partition.assignments[node] = proposal.D₂
    end

    # Replace (do not mutate shared BitSets from a CoW clone)
    partition.dist_nodes[proposal.D₁] = proposal.D₁_nodes
    partition.dist_nodes[proposal.D₂] = proposal.D₂_nodes

    update_partition_adjacency(partition, graph)
end

"""
    recom_chain_iter(graph::AbstractGraph,
                partition::AbstractPartition,
                pop_constraint::PopulationConstraint,
                num_steps::Int,
                scores::Array{S, 1};
                num_tries::Int=3,
                acceptance_fn::F=always_accept,
                rng::AbstractRNG=Random.default_rng(),
                no_self_loops::Bool=false,
                region_surcharges::Dict{String,Float64}=Dict{String,Float64}(),
                tree_method::Symbol=:kruskal,
                n_parallel::Int=1,
                cut_method::Symbol=:subtree_population) where {F<:Function,S<:AbstractScore}

Runs a Markov Chain for `num_steps` steps using ReCom. Returns an iterator
of `(Partition, score_vals)`. Note that `Partition` is mutable and will change
in-place with each iteration -- use `clone_for_update` if you wish to interact
with the `Partition` object outside of the for loop.

*Arguments:*
- graph:            `AbstractGraph`
- partition:        `AbstractPartition` with the plan information
- pop_constraint:   `PopulationConstraint`
- num_steps:        Number of steps to run the chain for
- scores:           Array of `AbstractScore`s to capture at each step
- num_tries:        num times to try getting a balanced cut from a subgraph
                    before giving up
- acceptance_fn:    A function generating a probability in [0, 1]
                    representing the likelihood of accepting the
                    proposal. Should accept a `Partition` as input.
- rng:              Random number generator. The user can pass in their
                    own; otherwise, we use the default RNG from Random. Must
                    implement the [AbstractRNG type](https://docs.julialang.org/en/v1/stdlib/Random/#Random.AbstractRNG)
                    (e.g. `Random.default_rng()` or `MersenneTwister(1234)`).
- no\\_self\\_loops: If this is true, then a failure to accept a new state
                    is not considered a self-loop; rather, the chain
                    simply generates new proposals until the acceptance
                    function is satisfied. BEWARE - this can create
                    infinite loops if the acceptance function is never
                    satisfied!
- region_surcharges: Optional region-boundary MST surcharges (Kruskal only)
- tree_method:      `:kruskal` or `:wilson`
- n_parallel:       concurrent proposal attempts per step (default `1`)
- cut_method:       `:subtree_population` (default) or `:edge_scan`
- progress_bar      If this is true, a progress bar will be printed to stdout.
"""
function recom_chain_iter(
    graph::AbstractGraph,
    partition::AbstractPartition,
    pop_constraint::PopulationConstraint,
    num_steps::Int,
    scores::Array{S,1};
    kwargs...,
) where {S<:AbstractScore}
    return _recom_chain_iter(graph, partition, pop_constraint, num_steps, scores; kwargs...)
end

function _recom_chain_iter end # ResumableFunctions workaround
@resumable function _recom_chain_iter(
    graph::BaseGraph,
    partition::Partition,
    pop_constraint::PopulationConstraint,
    num_steps::Int,
    scores::Array{S,1};
    num_tries::Int = 3,
    acceptance_fn::F = always_accept,
    rng::AbstractRNG = Random.default_rng(),
    no_self_loops::Bool = false,
    progress_bar = true,
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
    tree_method::Symbol = :kruskal,
    n_parallel::Int = 1,
    cut_method::Symbol = :subtree_population,
) where {F<:Function,S<:AbstractScore}
    if progress_bar
        iter = ProgressBar(1:num_steps)
    else
        iter = 1:num_steps
    end

    for steps_taken in iter
        step_completed = false
        while !step_completed
            proposal = get_valid_proposal(
                graph,
                partition,
                pop_constraint,
                rng,
                num_tries;
                region_surcharges = region_surcharges,
                tree_method = tree_method,
                n_parallel = n_parallel,
                cut_method = cut_method,
            )
            custom_acceptance = acceptance_fn !== always_accept
            update_partition!(partition, graph, proposal, custom_acceptance)
            if custom_acceptance && !satisfies_acceptance_fn(partition, acceptance_fn, rng)
                # go back to the previous partition
                partition = partition.parent
                # if user specifies this behavior, we do not increment the steps
                # taken if the acceptance function fails.
                if !no_self_loops
                    score_vals = score_partition_from_proposal(graph, partition, proposal, scores)
                    @yield partition, score_vals
                    step_completed = true
                end
            else
                score_vals = score_partition_from_proposal(graph, partition, proposal, scores)
                @yield partition, score_vals
                step_completed = true
            end
        end
    end
end

"""
    recom_chain(graph::AbstractGraph,
                partition::AbstractPartition,
                pop_constraint::PopulationConstraint,
                num_steps::Int,
                scores::Array{S, 1};
                num_tries::Int=3,
                acceptance_fn::F=always_accept,
                rng::AbstractRNG=Random.default_rng(),
                no_self_loops::Bool=false,
                region_surcharges::Dict{String,Float64}=Dict{String,Float64}(),
                tree_method::Symbol=:kruskal,
                n_parallel::Int=1,
                cut_method::Symbol=:subtree_population)::ChainScoreData where {F<:Function, S<:AbstractScore}

Runs a Markov Chain for `num_steps` steps using ReCom. Returns a `ChainScoreData`
object which can be queried to retrieve the values of every score at each
step of the chain.

*Arguments:*
- graph:            `AbstractGraph`
- partition:        `AbstractPartition` with the plan information
- pop_constraint:   `PopulationConstraint`
- num_steps:        Number of steps to run the chain for
- scores:           Array of `AbstractScore`s to capture at each step
- num_tries:        num times to try getting a balanced cut from a subgraph
                    before giving up
- acceptance_fn:    A function generating a probability in [0, 1]
                    representing the likelihood of accepting the
                    proposal. Should accept a `Partition` as input.
- rng:              Random number generator. The user can pass in their
                    own; otherwise, we use the default RNG from Random. Must
                    implement the [AbstractRNG type](https://docs.julialang.org/en/v1/stdlib/Random/#Random.AbstractRNG)
                    (e.g. `Random.default_rng()` or `MersenneTwister(1234)`).
- no\\_self\\_loops: If this is true, then a failure to accept a new state
                    is not considered a self-loop; rather, the chain
                    simply generates new proposals until the acceptance
                    function is satisfied. BEWARE - this can create
                    infinite loops if the acceptance function is never
                    satisfied!
- region_surcharges: Optional region-boundary MST surcharges (Kruskal only)
- tree_method:      `:kruskal` or `:wilson`
- n_parallel:       concurrent proposal attempts per step (default `1`)
- cut_method:       `:subtree_population` (default) or `:edge_scan`
- progress_bar      If this is true, a progress bar will be printed to stdout.
"""
function recom_chain(
    graph::AbstractGraph,
    partition::AbstractPartition,
    pop_constraint::PopulationConstraint,
    num_steps::Int,
    scores::Array{S,1};
    num_tries::Int = 3,
    acceptance_fn::F = always_accept,
    rng::AbstractRNG = Random.default_rng(),
    no_self_loops::Bool = false,
    progress_bar = true,
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
    tree_method::Symbol = :kruskal,
    n_parallel::Int = 1,
    cut_method::Symbol = :subtree_population,
)::ChainScoreData where {F<:Function,S<:AbstractScore}
    return _recom_chain(
        graph,
        partition,
        pop_constraint,
        num_steps,
        scores;
        num_tries = num_tries,
        acceptance_fn = acceptance_fn,
        rng = rng,
        no_self_loops = no_self_loops,
        progress_bar = progress_bar,
        region_surcharges = region_surcharges,
        tree_method = tree_method,
        n_parallel = n_parallel,
        cut_method = cut_method,
    )
end

function _recom_chain(
    graph::BaseGraph,
    partition::Partition,
    pop_constraint::PopulationConstraint,
    num_steps::Int,
    scores::Array{S,1};
    num_tries::Int = 3,
    acceptance_fn::F = always_accept,
    rng::AbstractRNG = Random.default_rng(),
    no_self_loops::Bool = false,
    progress_bar = true,
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
    tree_method::Symbol = :kruskal,
    n_parallel::Int = 1,
    cut_method::Symbol = :subtree_population,
)::ChainScoreData where {F<:Function,S<:AbstractScore}
    first_scores = score_initial_partition(graph, partition, scores)
    chain_scores = ChainScoreData(deepcopy(scores), [first_scores])

    for (_, score_vals) in _recom_chain_iter(
        graph,
        partition,
        pop_constraint,
        num_steps,
        scores;
        num_tries,
        acceptance_fn,
        rng,
        no_self_loops,
        progress_bar,
        region_surcharges,
        tree_method,
        n_parallel,
        cut_method,
    )
        push!(chain_scores.step_values, score_vals)
    end

    return chain_scores
end
