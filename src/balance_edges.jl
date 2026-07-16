"""
    kruskal_mst(graph::AbstractGraph,
                edges::Array{Int, 1},
                nodes::Array{Int, 1},
                weights::Array{Float64, 1})::BitSet

Generates and returns a minimum spanning tree from the subgraph induced
by `edges` and `nodes`, using Kruskal's MST algorithm. The `edges` are weighted
by `weights`.

## Note:
The `graph` represents the entire graph of the plan, where as `edges` and
`nodes` represent only the sub-graph on which we want to draw the MST.

*Arguments:*
- graph: Underlying Graph object
- edges: Array of edges of the sub-graph
- nodes: Set of nodes of the sub-graph
- weights: Array of weights of `length(edges)` where `weights[i]` is the
           weight of `edges[i]`

*Returns* a BitSet of edges that form a mst.
"""
function kruskal_mst(
    graph::AbstractGraph,
    edges::Array{Int,1},
    nodes::Array{Int,1},
    weights::Array{Float64,1},
)::BitSet
    num_nodes = length(nodes)
    srcs = edge_src(graph)
    dsts = edge_dst(graph)

    # sort the edges arr by their weights
    sorted_indices = sortperm(weights)
    sorted_edges = edges[sorted_indices]

    mst = BitSet()
    connected_vs = DisjointSets{Int}(nodes)

    for edge in sorted_edges
        if !in_same_set(connected_vs, srcs[edge], dsts[edge])
            union!(connected_vs, srcs[edge], dsts[edge])
            push!(mst, edge)
            (length(mst) >= num_nodes - 1) && break
        end
    end
    return mst
end

"""
    build_mst_weights!(weights::Vector{Float64},
                       graph::AbstractGraph,
                       edges::Array{Int,1},
                       rng::AbstractRNG;
                       region_surcharges::Dict{String,Float64}=Dict{String,Float64}())

Fill `weights` (length = `length(edges)`) with:
`rand(rng) + edge_penalties[e] + Σ surcharge[col]` when endpoints differ in region `col`.
"""
function build_mst_weights!(
    weights::Vector{Float64},
    graph::AbstractGraph,
    edges::Array{Int,1},
    rng::AbstractRNG;
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
)
    length(weights) == length(edges) || throw(
        ArgumentError("weights length must match edges length"),
    )
    penalties = edge_penalties(graph)
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    surcharge_cols = collect(keys(region_surcharges))
    surcharge_vals = [region_surcharges[c] for c in surcharge_cols]
    region_vecs = [region_ids(graph, c) for c in surcharge_cols]

    @inbounds for i = 1:length(edges)
        e = edges[i]
        w = rand(rng) + penalties[e]
        u = srcs[e]
        v = dsts[e]
        for j = 1:length(surcharge_cols)
            if region_vecs[j][u] != region_vecs[j][v]
                w += surcharge_vals[j]
            end
        end
        weights[i] = w
    end
    return weights
end

"""
    weighted_kruskal_mst(graph::AbstractGraph,
                         edges::Array{Int,1},
                         nodes::Array{Int,1},
                         rng::AbstractRNG=Random.default_rng();
                         region_surcharges::Dict{String,Float64}=Dict{String,Float64}())::BitSet

Kruskal MST with random base weights plus graph edge penalties and optional
region-boundary surcharges.
"""
function weighted_kruskal_mst(
    graph::AbstractGraph,
    edges::Array{Int,1},
    nodes::Array{Int,1},
    rng::AbstractRNG = Random.default_rng();
    region_surcharges::Dict{String,Float64} = Dict{String,Float64}(),
)::BitSet
    weights = Vector{Float64}(undef, length(edges))
    build_mst_weights!(
        weights,
        graph,
        edges,
        rng;
        region_surcharges = region_surcharges,
    )
    return kruskal_mst(graph, edges, nodes, weights)
end

"""
    weighted_kruskal_mst(graph, edges, nodes, weights::AbstractVector{<:Real})

Compatibility overload matching the pre-0.2.0 / downstream API that passes
precomputed edge weights (e.g. UTGC `MST_FUNC`).
"""
function weighted_kruskal_mst(
    graph::AbstractGraph,
    edges::Array{Int,1},
    nodes::Array{Int,1},
    weights::AbstractVector{<:Real},
)::BitSet
    return kruskal_mst(graph, edges, nodes, Float64.(weights))
end

"""
    random_kruskal_mst(graph::AbstractGraph,
                       edges::Array{Int, 1},
                       nodes::Array{Int, 1},
                       rng::AbstractRNG=Random.default_rng())

Generates and returns a random minimum spanning tree from the subgraph induced
by `edges` and `nodes`, using Kruskal's MST algorithm.

## Note:
The `graph` represents the entire graph of the plan, where as `edges` and
`nodes` represent only the sub-graph on which we want to draw the MST.

*Arguments:*
- graph: Underlying Graph object
- edges: Array of edges of the sub-graph
- nodes: Set of nodes of the sub-graph
- rng: A random number generator that implements the [AbstractRNG type](https://docs.julialang.org/en/v1/stdlib/Random/#Random.AbstractRNG) (e.g. `Random.default_rng()` or `MersenneTwister(1234)`)

*Returns* a BitSet of edges that form a mst.
"""
function random_kruskal_mst(
    graph::AbstractGraph,
    edges::Array{Int,1},
    nodes::Array{Int,1},
    rng::AbstractRNG = Random.default_rng(),
)::BitSet
    weights = rand(rng, length(edges))
    return kruskal_mst(graph, edges, nodes, weights)
end
