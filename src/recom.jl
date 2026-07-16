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
    get_valid_proposal(graph::AbstractGraph,
                       partition::AbstractPartition,
                       pop_constraint::PopulationConstraint,
                       rng::AbstractRNG,
                       num_tries::Int=3;
                       region_surcharges::Dict{String,Float64}=Dict{String,Float64}())

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
"""
function get_valid_proposal(
    graph::AbstractGraph,
    partition::AbstractPartition,
    pop_constraint::PopulationConstraint,
    rng::AbstractRNG,
    num_tries::Int = 3;
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
)
    use_weighted =
        !isempty(region_surcharges) || any(!iszero, edge_penalties(graph))

    while true
        D₁, D₂, sg_edges, sg_nodes = sample_subgraph(graph, partition, rng)
        sg_node_list = collect(sg_nodes)

        for _ = 1:num_tries
            mst_edges = if use_weighted
                weighted_kruskal_mst(
                    graph,
                    sg_edges,
                    sg_node_list,
                    rng;
                    region_surcharges = region_surcharges,
                )
            else
                random_kruskal_mst(graph, sg_edges, sg_node_list, rng)
            end

            # see if we can get a population-balanced cut in this mst
            proposal = get_balanced_proposal(
                graph,
                mst_edges,
                sg_nodes,
                partition,
                pop_constraint,
                D₁,
                D₂,
            )

            if proposal isa RecomProposal
                return proposal
            end
        end
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
                region_surcharges::Dict{String,Float64}=Dict{String,Float64}()) where {F<:Function,S<:AbstractScore}

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
- region_surcharges: Optional region-boundary MST surcharges
- progress_bar      If this is true, a progress bar will be printed to stdout.
"""
function recom_chain_iter end # this is a workaround (https://github.com/BenLauwens/ResumableFunctions.jl/issues/45)
@resumable function recom_chain_iter(
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
                region_surcharges::Dict{String,Float64}=Dict{String,Float64}())::ChainScoreData where {F<:Function, S<:AbstractScore}

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
- region_surcharges: Optional region-boundary MST surcharges
- progress_bar      If this is true, a progress bar will be printed to stdout.
"""
function recom_chain(
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
)::ChainScoreData where {F<:Function,S<:AbstractScore}
    first_scores = score_initial_partition(graph, partition, scores)
    chain_scores = ChainScoreData(deepcopy(scores), [first_scores])

    for (_, score_vals) in recom_chain_iter(
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
    )
        push!(chain_scores.step_values, score_vals)
    end

    return chain_scores
end
