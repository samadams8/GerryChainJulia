"""
    MSTScratch

Reusable buffers for Kruskal / weighted MST construction.
"""
const _NULL_REGION_ID = UInt32(0)

mutable struct MSTScratch
    weights::Vector{Float64}
    idx::Vector{Int}
    parent::Vector{Int}
    rank::Vector{Int}
    mst_edges::Vector{Int}
    max_node::Int
end

function MSTScratch(max_edges::Int = 0, max_node::Int = 0)
    return MSTScratch(
        Vector{Float64}(undef, max_edges),
        Vector{Int}(undef, max_edges),
        zeros(Int, max_node),
        zeros(Int, max_node),
        Int[],
        max_node,
    )
end

function _ensure_mst_scratch!(scratch::MSTScratch, n_edges::Int, max_node::Int)
    if length(scratch.weights) < n_edges
        resize!(scratch.weights, n_edges)
        resize!(scratch.idx, n_edges)
    end
    if scratch.max_node < max_node || length(scratch.parent) < max_node
        scratch.parent = zeros(Int, max_node)
        scratch.rank = zeros(Int, max_node)
        scratch.max_node = max_node
    end
    empty!(scratch.mst_edges)
    return scratch
end

"""
    SubtreeCutScratch

Reusable buffers for `get_balanced_proposal_subtree_population`.
"""
mutable struct SubtreeCutScratch
    adj::Vector{Vector{Int}}
    parent::Vector{Int}
    subpop::Vector{Int}
    order::Vector{Int}
    stack::Vector{Int}
    seen::BitVector
    max_node::Int
end

function SubtreeCutScratch(max_node::Int = 0)
    return SubtreeCutScratch(
        [Int[] for _ = 1:max_node],
        zeros(Int, max_node),
        zeros(Int, max_node),
        Int[],
        Int[],
        falses(max_node),
        max_node,
    )
end

function _ensure_subtree_cut_scratch!(scratch::SubtreeCutScratch, max_node::Int)
    if scratch.max_node < max_node || length(scratch.parent) < max_node
        old = scratch.max_node
        resize!(scratch.adj, max_node)
        for i = (old + 1):max_node
            scratch.adj[i] = Int[]
        end
        scratch.parent = zeros(Int, max_node)
        scratch.subpop = zeros(Int, max_node)
        scratch.seen = falses(max_node)
        scratch.max_node = max_node
    end
    # Clear adjacency lists touched on the previous run.
    @inbounds for u in scratch.order
        empty!(scratch.adj[u])
    end
    empty!(scratch.order)
    empty!(scratch.stack)
    fill!(scratch.parent, 0)
    # The `seen` bitvector is reset only when collecting a component.
    return scratch
end

@inline function _uf_find!(parent::Vector{Int}, x::Int)
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    return x
end

@inline function _uf_union!(parent::Vector{Int}, rank::Vector{Int}, a::Int, b::Int)
    ra, rb = _uf_find!(parent, a), _uf_find!(parent, b)
    ra == rb && return false
    if rank[ra] < rank[rb]
        parent[ra] = rb
    elseif rank[ra] > rank[rb]
        parent[rb] = ra
    else
        parent[rb] = ra
        rank[ra] += 1
    end
    return true
end

"""
    kruskal_mst!(scratch, graph, edges, nodes, weights) -> BitSet

In-place Kruskal using `scratch` buffers. `weights[1:length(edges)]` are used.
"""
function kruskal_mst!(
    scratch::MSTScratch,
    graph::AbstractGraph,
    edges::Vector{Int},
    nodes::Vector{Int},
    weights::AbstractVector{Float64},
)::BitSet
    n_edges = length(edges)
    num_nodes = length(nodes)
    max_node = isempty(nodes) ? 0 : maximum(nodes)
    _ensure_mst_scratch!(scratch, n_edges, max_node)
    srcs = edge_src(graph)
    dsts = edge_dst(graph)

    parent, rank = scratch.parent, scratch.rank
    @inbounds for v in nodes
        parent[v] = v
        rank[v] = 0
    end

    sortperm!(view(scratch.idx, 1:n_edges), view(weights, 1:n_edges))
    empty!(scratch.mst_edges)
    needed = num_nodes - 1
    @inbounds for t = 1:n_edges
        edge = edges[scratch.idx[t]]
        u, v = srcs[edge], dsts[edge]
        if _uf_union!(parent, rank, u, v)
            push!(scratch.mst_edges, edge)
            length(scratch.mst_edges) >= needed && break
        end
    end
    return BitSet(scratch.mst_edges)
end

"""
    kruskal_mst(graph::AbstractGraph,
                edges::Vector{Int},
                nodes::Vector{Int},
                weights::AbstractVector{<:Real})::BitSet

Build a Minimum Spanning Tree using Kruskal's algorithm and user-provided custom weights.
Allocates temporary scratch buffers.
"""
function kruskal_mst(
    graph::AbstractGraph,
    edges::Vector{Int},
    nodes::Vector{Int},
    weights::AbstractVector{<:Real},
)::BitSet
    scratch = MSTScratch(length(edges), isempty(nodes) ? 0 : maximum(nodes))
    return kruskal_mst!(scratch, graph, edges, nodes, Float64.(weights))
end

"""
    build_mst_weights!(weights::Vector{Float64},
                       graph::BaseGraph,
                       edges::Vector{Int},
                       rng::AbstractRNG)

Populate `weights` with random weights plus any pre-configured base weights
(surcharges and penalties) from the graph.
"""
function build_mst_weights!(
    weights::Vector{Float64},
    graph::BaseGraph,
    edges::Vector{Int},
    rng::AbstractRNG,
)
    length(weights) >= length(edges) || throw(
        ArgumentError("weights length must be at least edges length"),
    )
    base = graph._mst_base_weights[]
    # Splitting the loop based on whether base is nothing ensures type stability
    # by avoiding Union-splitting inside the hot loop.
    if base === nothing
        penalties = edge_penalties(graph)
        @inbounds for i = 1:length(edges)
            e = edges[i]
            weights[i] = rand(rng) + penalties[e]
        end
    else
        @inbounds for i = 1:length(edges)
            e = edges[i]
            weights[i] = rand(rng) + base[e]
        end
    end
    return weights
end

"""
    build_mst_weights!(scratch::MSTScratch,
                       graph::BaseGraph,
                       edges::Vector{Int},
                       rng::AbstractRNG)

Populate `scratch.weights` using scratch buffers.
"""
function build_mst_weights!(
    scratch::MSTScratch,
    graph::BaseGraph,
    edges::Vector{Int},
    rng::AbstractRNG,
)
    max_n = 0
    srcs, dsts = edge_src(graph), edge_dst(graph)
    @inbounds for e in edges
        max_n = max(max_n, srcs[e], dsts[e])
    end
    _ensure_mst_scratch!(scratch, length(edges), max_n)
    return build_mst_weights!(
        scratch.weights,
        graph,
        edges,
        rng,
    )
end

"""
    _kruskal_mst(graph::BaseGraph,
                 edges::Vector{Int},
                 nodes::Vector{Int},
                 rng::AbstractRNG;
                 scratch=nothing)

Internal Kruskal spanning tree implementation. Reads cached edge penalties
and region surcharges from `graph`. Fallback is uniform random if unconfigured.
"""
function _kruskal_mst(
    graph::BaseGraph,
    edges::Vector{Int},
    nodes::Vector{Int},
    rng::AbstractRNG = Random.default_rng();
    scratch::Union{MSTScratch,Nothing} = nothing,
)::BitSet
    max_node = isempty(nodes) ? 0 : maximum(nodes)
    if scratch === nothing
        scratch = MSTScratch(length(edges), max_node)
    end
    build_mst_weights!(
        scratch,
        graph,
        edges,
        rng,
    )
    return kruskal_mst!(scratch, graph, edges, nodes, scratch.weights)
end

"""
    wilson_ust(graph, edges, nodes, rng) -> BitSet

Wilson's algorithm: uniform random spanning tree on the induced subgraph.
"""
function wilson_ust(
    graph::AbstractGraph,
    edges::Vector{Int},
    nodes::Vector{Int},
    rng::AbstractRNG = Random.default_rng(),
)::BitSet
    length(nodes) <= 1 && return BitSet()

    node_set = BitSet(nodes)
    max_node = maximum(nodes)
    adj = [Tuple{Int,Int}[] for _ = 1:max_node]
    srcs, dsts = edge_src(graph), edge_dst(graph)
    for e in edges
        u, v = srcs[e], dsts[e]
        if (u in node_set) && (v in node_set)
            push!(adj[u], (v, e))
            push!(adj[v], (u, e))
        end
    end

    in_tree = falses(max_node)
    next_node = zeros(Int, max_node)
    next_edge = zeros(Int, max_node)
    root = nodes[rand(rng, 1:length(nodes))]
    in_tree[root] = true
    tree_edges = Int[]

    for start in nodes
        in_tree[start] && continue
        u = start
        while !in_tree[u]
            nbrs = adj[u]
            isempty(nbrs) && throw(ArgumentError("Induced subgraph is disconnected"))
            v, e = nbrs[rand(rng, 1:length(nbrs))]
            next_node[u] = v
            next_edge[u] = e
            u = v
        end
        u = start
        while !in_tree[u]
            push!(tree_edges, next_edge[u])
            in_tree[u] = true
            u = next_node[u]
        end
    end
    return BitSet(tree_edges)
end
