abstract type AbstractGraph end

struct BaseGraph <: AbstractGraph
    num_nodes::Int
    num_edges::Int
    total_pop::Int
    populations::Array{Int,1}  # Of length num_nodes.
    adj_matrix::SparseMatrixCSC{Int,Int}
    edge_src::Array{Int,1}  # Of length num_edges.
    edge_dst::Array{Int,1}  # Of length num_edges.
    neighbors::Array{Array{Int64,1},1}
    simple_graph::SimpleGraph  # The base SimpleGraph, if we need it.
    attributes::Array{Dict{String,Any}}
    edge_penalties::Vector{Float64}  # Of length num_edges; index = edge ID.
    region_cols::Dict{String,Vector{UInt32}}  # Dense region ID columns.
    _attr_cache::Dict{String,Vector{Float64}}  # Lazy dense attribute columns.
    _mst_base_weights::Base.RefValue{Union{Vector{Float64},Nothing}}
end

"""
    BaseGraph(num_nodes, num_edges, total_pop, populations, adj_matrix,
              edge_src, edge_dst, neighbors, simple_graph, attributes)

Compatibility constructor matching the pre-0.2.0 10-field layout.
Initializes `edge_penalties` to zeros and `region_cols` to an empty dict.

`attributes` may be any vector of string-keyed dicts (e.g. `Dict{String,String}`);
values are coerced to `Dict{String,Any}` for storage.
"""
function BaseGraph(
    num_nodes::Int,
    num_edges::Int,
    total_pop::Int,
    populations::Array{Int,1},
    adj_matrix::SparseMatrixCSC{Int,Int},
    edge_src::Array{Int,1},
    edge_dst::Array{Int,1},
    neighbors::Array{<:AbstractVector{<:Integer},1},
    simple_graph::SimpleGraph,
    attributes::AbstractVector{<:AbstractDict{<:AbstractString,<:Any}},
)
    attrs = Dict{String,Any}[Dict{String,Any}(d) for d in attributes]
    nbrs = Array{Array{Int64,1},1}([Int64.(n) for n in neighbors])
    return BaseGraph(
        num_nodes,
        num_edges,
        total_pop,
        populations,
        adj_matrix,
        edge_src,
        edge_dst,
        nbrs,
        simple_graph,
        attrs,
        zeros(Float64, num_edges),
        Dict{String,Vector{UInt32}}(),
        Dict{String,Vector{Float64}}(),
        Base.RefValue{Union{Vector{Float64},Nothing}}(nothing),
    )
end

# AbstractGraph accessors (default concrete implementations).

num_nodes(g::BaseGraph) = g.num_nodes
num_edges(g::BaseGraph) = g.num_edges
total_pop(g::BaseGraph) = g.total_pop
populations(g::BaseGraph) = g.populations
edge_src(g::BaseGraph) = g.edge_src
edge_dst(g::BaseGraph) = g.edge_dst
neighbors(g::BaseGraph) = g.neighbors
edge_penalties(g::BaseGraph) = g.edge_penalties

"""
    has_region(g::AbstractGraph, col::AbstractString)::Bool

Return whether dense region column `col` is registered on the graph.
"""
has_region(g::BaseGraph, col::AbstractString) = haskey(g.region_cols, String(col))

"""
    region_ids(g::AbstractGraph, col::AbstractString)::AbstractVector{UInt32}

Return the dense region-id vector for column `col`.
"""
function region_ids(g::BaseGraph, col::AbstractString)
    key = String(col)
    haskey(g.region_cols, key) ||
        throw(ArgumentError("Region column \"$key\" is not registered on this graph."))
    return g.region_cols[key]
end

"""
    read_table(filepath::AbstractString)::Shapefile.Table

Read table from shapefile. If a .shp and a .dbf file of the same name
are not found, then we throw an error.
"""
function read_table(filepath::AbstractString)::Shapefile.Table
    prefix = splitext(filepath)[1]
    if !(isfile(prefix * ".shp") || isfile(prefix * ".SHP")) ||
       !(isfile(prefix * ".dbf") || isfile(prefix * ".DBF"))
        throw(
            ArgumentError(
                "Error when processing filepath as shapefile: to read a graph from a shapefile, we require a .shp and .dbf file of the same name in the same folder.",
            ),
        )
    end
    return Shapefile.Table(prefix)
end

"""
    all_node_properties(table::Shapefile.Table)::Array{Dict{String, Any}}

*Returns* an Array of Dictionaries. Each dictionary corresponds to one
node in the graph.
"""
function all_node_properties(table::Shapefile.Table)::Array{Dict{String,Any}}
    properties = propertynames(table) # returns array of symbols
    string_keys = String.(properties) # convert by broadcasting

    # Internal function because we want to use both properties and values.
    function get_node_properties(row::Shapefile.Row)
        values = map(p -> getproperty(row, p), properties)
        return Dict(string_keys .=> values)
    end

    return get_node_properties.(table)
end


"""
    get_attribute_by_key(node_attributes::Array,
                         column_name::String,
                         process_value::Function=identity)::Array

*Returns* an array whose values correspond to the value of an attribute
for each node.

*Arguments:*
- node_attributes :   An array of Dict{String, Any}, where each dictionary
                      represents a mapping from attribute names to values
                      for a particular node.
- column_name     :   The name of the attribute (i.e., the key of the
                      attribute in the dictionaries)
- process_value   :   An optional argument that processes the raw value
"""
function get_attribute_by_key(
    node_attributes::Array,
    column_name::String,
    process_value::Function = identity,
)::Array
    return [process_value(n[column_name]) for n in node_attributes]
end

"""
    population_to_int(raw_value::Number)::Int

Tiny helper function to coerce population counts (whether they are
ints or floats) to int.
"""
function population_to_int(raw_value::Number)::Int
    return raw_value isa Int ? raw_value : convert(Int, round(raw_value))
end

"""
    get_node_coordinates(row::Shapefile.Row)::Vector{Vector{Vector{Vector{Float64}}}}

Construct an array of LibGEOS.Polygons from the given coordinates. The
coordinates are structured in the following way. Each element in the
outermost array represents one polygon. (One node can be made up of
multiple polygons).
[
    [  # One array = one polygon.
        [  # Points corresponding to outer ring.
            [1.0, 2.0],  # Single x,y coordinate of a point.
            ...
        ],
        [  # Points corresponding to a hole in polygon.
            [1.5, 1.7],
            ...
        ],
        ...  # Any other subsequent arrays would
             # correspond to other holes.
    ],
    ...  # Subsequent arrays correspond to other
         # polygons.
]
"""
function get_node_coordinates(row::Shapefile.Row)::Vector{Vector{Vector{Vector{Float64}}}}
    return LibGEOS.GeoInterface.coordinates(getfield(row, :geometry))
end


"""
    graph_from_shp(filepath::AbstractString,
                   pop_col::AbstractString,
                   adjacency::String="rook")::BaseGraph

Constructs BaseGraph from .shp file.
"""
function graph_from_shp(
    filepath::AbstractString,
    pop_col::AbstractString,
    adjacency::String = "rook";
    region_columns::Vector{String} = String[],
)::BaseGraph
    table = read_table(filepath)

    attributes = all_node_properties(table)
    coords = get_node_coordinates.(table)
    # These will be used in the adjacency method.
    node_polys = polygon_array.(coords)
    node_mbrs = min_bounding_rect.(coords)

    graph = simple_graph_from_polygons(node_polys, node_mbrs, adjacency)

    # Edge `i` connects nodes edge_src[i] and edge_dst[i].
    edge_src, edge_dst = edges_from_graph(graph)
    # Each entry in adj_matrix is the edge ID that connects the two nodes.
    adj_matrix = adjacency_matrix_from_graph(graph)
    neighbors = neighbors_from_graph(graph)

    populations = get_attribute_by_key(attributes, pop_col, population_to_int)
    total_pop = sum(populations)
    n_edges = ne(graph)
    region_cols = build_region_cols(attributes, region_columns)

    return BaseGraph(
        nv(graph),
        n_edges,
        total_pop,
        populations,
        adj_matrix,
        edge_src,
        edge_dst,
        neighbors,
        graph,
        attributes,
        zeros(Float64, n_edges),
        region_cols,
        Dict{String,Vector{Float64}}(),
        Base.RefValue{Union{Vector{Float64},Nothing}}(nothing),
    )
end

"""
    edges_from_graph(graph::SimpleGraph)

Extract edges from graph. Returns two arrays; the first contains the
indices of the source nodes and the second contains the indices
of the destination nodes.
"""
function edges_from_graph(graph::SimpleGraph)
    num_edges = ne(graph)

    # Edge `i` connects nodes edge_src[i] and edge_dst[i].
    edge_src = zeros(Int, num_edges)
    edge_dst = zeros(Int, num_edges)

    for (index, edge) in enumerate(edges(graph))
        edge_src[index] = src(edge)
        edge_dst[index] = dst(edge)
    end
    return edge_src, edge_dst
end

"""
    adjacency_matrix_from_graph(graph::SimpleGraph)

Extract sparse adjacency matrix from graph.
"""
function adjacency_matrix_from_graph(graph::SimpleGraph)
    # Each entry in adj_matrix is the edge ID that connects the two nodes.
    num_nodes = nv(graph)
    adj_matrix = spzeros(Int, num_nodes, num_nodes)
    for (index, edge) in enumerate(edges(graph))
        adj_matrix[src(edge), dst(edge)] = index
        adj_matrix[dst(edge), src(edge)] = index
    end
    return adj_matrix
end

"""
    neighbors_from_graph(graph::SimpleGraph)

Extract each node's neighbors from graph.
"""
function neighbors_from_graph(graph::SimpleGraph)
    # Each entry in the returned array is a list of neighbors of the corresponding node.
    neighbors = [Int[] for n = 1:nv(graph)]
    for (index, edge) in enumerate(edges(graph))
        push!(neighbors[src(edge)], dst(edge))
        push!(neighbors[dst(edge)], src(edge))
    end
    return neighbors
end

"""
    graph_from_json(filepath::AbstractString,
                    pop_col::AbstractString)::BaseGraph

*Arguments:*
- filepath:       file path to the .json file that contains the graph.
                  This file is expected to be generated by the `Graph.to_json()`
                  function of the Python implementation of Gerrychain. [1]
                  We assume that the JSON file has the structure of a dictionary
                  where (1) the key "nodes" yields an array of dictionaries
                  of node attributes, (2) the key "adjacency" yields an
                  array of edges (represented as dictionaries), and (3)
                  the key "id" within the edge dictionary indicates the
                  destination node of the edge.
- pop_col:        the node attribute key whose accompanying value is the
                  population of that node

[1]: https://github.com/mggg/GerryChain/blob/c87da7e69967880abc99b781cd37468b8cb18815/gerrychain/graph/graph.py#L38
"""
function graph_from_json(
    filepath::AbstractString,
    pop_col::AbstractString;
    region_columns::Vector{String} = String[],
)::BaseGraph
    raw_graph = JSON.parsefile(filepath)
    nodes = raw_graph["nodes"]
    num_nodes = length(nodes)

    # Get populations.
    populations = get_attribute_by_key(nodes, pop_col, population_to_int)
    total_pop = sum(populations)

    # Generate the base SimpleGraph.
    simple_graph = SimpleGraph(num_nodes)
    for (index, edges) in enumerate(raw_graph["adjacency"])
        for edge in edges
            if edge["id"] + 1 > index
                add_edge!(simple_graph, index, edge["id"] + 1)
            end
        end
    end

    num_edges = ne(simple_graph)

    # Edge `i` connects nodes edge_src[i] and edge_dst[i].
    edge_src, edge_dst = edges_from_graph(simple_graph)
    # Each entry in adj_matrix is the edge ID that connects the two nodes.
    adj_matrix = adjacency_matrix_from_graph(simple_graph)
    neighbors = neighbors_from_graph(simple_graph)

    # Get attributes.
    attributes = get_attributes(nodes)
    region_cols = build_region_cols(attributes, region_columns)

    return BaseGraph(
        num_nodes,
        num_edges,
        total_pop,
        populations,
        adj_matrix,
        edge_src,
        edge_dst,
        neighbors,
        simple_graph,
        attributes,
        zeros(Float64, num_edges),
        region_cols,
        Dict{String,Vector{Float64}}(),
        Base.RefValue{Union{Vector{Float64},Nothing}}(nothing),
    )
end

"""
    BaseGraph(filepath::AbstractString,
              pop_col::AbstractString;
              adjacency::String="rook",
              region_columns::Vector{String}=String[])::BaseGraph

Builds the BaseGraph object. This is the underlying network of our
districts, and its properties are immutable i.e they will not change
from step to step in our Markov Chains.

*Arguments:*
- filepath:       A path to a .json or .shp file which contains the
                  information needed to construct the graph.
- pop_col:        the node attribute key whose accompanying value is the
                  population of that node
- adjacency:      (Only used if the user specifies a filepath to a .shp
                  file.) Should be either "queen" or "rook"; "rook" by default.
- region_columns: Attribute column names to materialize as dense
                  `UInt32` region-id vectors for region-aware ReCom.
"""
function BaseGraph(
    filepath::AbstractString,
    pop_col::AbstractString;
    adjacency::String = "rook",
    region_columns::Vector{String} = String[],
)::BaseGraph
    extension = uppercase(splitext(filepath)[2])
    if uppercase(extension) == ".JSON"
        return graph_from_json(filepath, pop_col; region_columns = region_columns)
    elseif uppercase(extension) == ".SHP"
        return graph_from_shp(
            filepath,
            pop_col,
            adjacency;
            region_columns = region_columns,
        )
    else
        throw(
            DomainError(
                filepath,
                "Filepath must lead to valid JSON file or valid .shp/.dbf file.",
            ),
        )
    end
end

"""
    get_attributes(nodes::Array{Any, 1})

*Returns* an array of dicts `attributes` of length `length(nodes)` where
the attributes of the `nodes[i]` is at `attributes[i]` as a dictionary.
"""
function get_attributes(nodes::Array{Any,1})
    attributes = Array{Dict{String,Any}}(undef, length(nodes))
    for (index, node) in enumerate(nodes)
        attributes[index] = node
    end
    return attributes
end

"""
    induced_subgraph_edges(graph::BaseGraph,
                           vlist::Array{Int, 1})::Array{Int, 1}

*Returns* a list of edges of the subgraph induced by `vlist`, which is an array
of vertices.
"""
function induced_subgraph_edges(
    graph::BaseGraph,
    vlist::Array{Int,1},
)::Array{Int,1}
    n = num_nodes(graph)
    in_v = falses(n)
    @inbounds for v in vlist
        (v < 1 || v > n) && throw(ArgumentError("Vertex $v out of range 1:$n"))
        in_v[v] && throw(ArgumentError("Vertices in subgraph list must be unique"))
        in_v[v] = true
    end
    edges = Int[]
    sizehint!(edges, length(vlist))
    @inbounds for u in vlist
        for v in graph.neighbors[u]
            # undirected dual graph: emit each edge once
            if v > u && in_v[v]
                push!(edges, graph.adj_matrix[u, v])
            end
        end
    end
    return edges
end

"""
    get_subgraph_population(graph::AbstractGraph,
                            nodes::BitSet)::Int

*Arguments:*
- graph: Underlying graph object
- nodes: A Set of Ints

*Returns* the population of the subgraph induced by `nodes`.
"""
function get_subgraph_population(graph::AbstractGraph, nodes::BitSet)::Int
    pops = populations(graph)
    total = 0
    for node in nodes
        total += pops[node]
    end
    return total
end

"""
    _is_null_region_value(v) -> Bool

True for values that mean "no region": `missing`, `nothing`, or `""`.
Encoded as `UInt32(0)` by `encode_region_values`.
"""
function _is_null_region_value(v)::Bool
    return v === nothing || ismissing(v) || v == ""
end

"""
    encode_region_values(values)::Vector{UInt32}

Map arbitrary region labels to dense `UInt32` codes.
`missing`, `nothing`, and `""` map to sentinel `0` (no region).
All other labels get distinct 1-based dense codes.
"""
function encode_region_values(values)::Vector{UInt32}
    code_of = Dict{Any,UInt32}()
    coded = Vector{UInt32}(undef, length(values))
    next_code = UInt32(1)
    for (i, v) in enumerate(values)
        if _is_null_region_value(v)
            coded[i] = UInt32(0)
            continue
        end
        if !haskey(code_of, v)
            code_of[v] = next_code
            next_code += UInt32(1)
        end
        coded[i] = code_of[v]
    end
    return coded
end

"""
    build_region_cols(attributes, region_columns::Vector{String})

Build a dictionary of dense region-id columns from node attribute dicts.
"""
function build_region_cols(
    attributes,
    region_columns::Vector{String},
)::Dict{String,Vector{UInt32}}
    region_cols = Dict{String,Vector{UInt32}}()
    for col in region_columns
        values = get_attribute_by_key(attributes, col)
        region_cols[col] = encode_region_values(values)
    end
    return region_cols
end

"""
    add_region_column!(g::BaseGraph, name::AbstractString, values)

Register (or replace) a dense region column. `values` length must equal
`num_nodes(g)`. Non-integer values are encoded to dense codes via
`encode_region_values` (`missing` / `nothing` / `""` → `0`). Integer
inputs are cast to `UInt32` as-is (`0` remains the null-region sentinel).
"""
function add_region_column!(g::BaseGraph, name::AbstractString, values)
    length(values) == g.num_nodes || throw(
        ArgumentError(
            "Region column length ($(length(values))) must equal num_nodes ($(g.num_nodes)).",
        ),
    )
    if eltype(values) <: Integer
        g.region_cols[String(name)] = UInt32.(values)
    else
        g.region_cols[String(name)] = encode_region_values(values)
    end
    g._mst_base_weights[] = nothing
    return g
end

"""
    set_edge_penalty!(g::BaseGraph, u::Int, v::Int, w::Float64)

Set the penalty for the edge between nodes `u` and `v`.
"""
function set_edge_penalty!(g::BaseGraph, u::Int, v::Int, w::Float64)
    edge_id = g.adj_matrix[u, v]
    edge_id == 0 && throw(ArgumentError("No edge between nodes $u and $v."))
    g.edge_penalties[edge_id] = w
    g._mst_base_weights[] = nothing
    return g
end

"""
    set_edge_penalties_from_pairs!(g::BaseGraph, penalties)

Set edge penalties from a dictionary keyed by `(u, v)` node pairs
(order-insensitive) or from a vector indexed by edge id.
"""
function set_edge_penalties_from_pairs!(
    g::BaseGraph,
    penalties::AbstractDict{<:Tuple{Int,Int},<:Real},
)
    for ((u, v), w) in penalties
        set_edge_penalty!(g, u, v, Float64(w))
    end
    return g
end

function set_edge_penalties_from_pairs!(g::BaseGraph, penalties::AbstractVector{<:Real})
    length(penalties) == g.num_edges || throw(
        ArgumentError(
            "Penalty vector length ($(length(penalties))) must equal num_edges ($(g.num_edges)).",
        ),
    )
    for i = 1:g.num_edges
        g.edge_penalties[i] = Float64(penalties[i])
    end
    g._mst_base_weights[] = nothing
    return g
end

# Unified MST weights configuration.

"""
    configure_mst_weights!(graph::BaseGraph;
                           edge_penalties=nothing,
                           region_surcharges=Dict{String,Float64}())

Configure the deterministic component of MST edge weights and cache them
on the graph. The random jitter `rand(rng)` is added at tree-build time.

- `edge_penalties`: a vector of length `num_edges` or a `Dict{(u,v) => w}` or nothing.
  If omitted/nothing, the graph's current `edge_penalties` are used.
- `region_surcharges`: `Dict{String,Float64}` mapping region column name to
  surcharge added when an edge crosses a boundary in that column.
"""
function configure_mst_weights!(
    graph::BaseGraph;
    edge_penalties=nothing,
    region_surcharges::Dict{String,Float64}=Dict{String,Float64}(),
)
    if edge_penalties !== nothing
        set_edge_penalties_from_pairs!(graph, edge_penalties)
    end
    graph._mst_base_weights[] = _compute_mst_base_weights(graph, region_surcharges)
    return graph
end

"""
    _compute_mst_base_weights(graph::BaseGraph,
                              region_surcharges::Dict{String,Float64})::Vector{Float64}

Internal helper to precompute the deterministic part of the MST weights (penalties + surcharges)
across all edges of the graph.
"""
function _compute_mst_base_weights(
    graph::BaseGraph,
    region_surcharges::Dict{String,Float64},
)::Vector{Float64}
    penalties = edge_penalties(graph)
    base = copy(penalties)
    surcharge_cols = collect(keys(region_surcharges))
    if isempty(surcharge_cols)
        return base
    end
    srcs, dsts = edge_src(graph), edge_dst(graph)
    surcharge_vals = Float64[region_surcharges[c] for c in surcharge_cols]
    region_vecs = Vector{UInt32}[region_ids(graph, c) for c in surcharge_cols]
    @inbounds for e in 1:num_edges(graph)
        u, v = srcs[e], dsts[e]
        for j in 1:length(surcharge_cols)
            id_u, id_v = region_vecs[j][u], region_vecs[j][v]
            if id_u != _NULL_REGION_ID && id_v != _NULL_REGION_ID && id_u != id_v
                base[e] += surcharge_vals[j]
            end
        end
    end
    return base
end

# Lazy Float64 attribute column cache.

function _materialize_attribute_column(g::BaseGraph, key::String)::Vector{Float64}
    n = g.num_nodes
    col = Vector{Float64}(undef, n)
    @inbounds for i = 1:n
        raw = g.attributes[i][key]
        if raw isa Number
            col[i] = Float64(raw)
        else
            throw(
                ArgumentError(
                    "Attribute \"$key\" at node $i has type $(typeof(raw)); expected Number. " *
                    "Try coerce_aggregated_attributes! before scoring.",
                ),
            )
        end
    end
    return col
end

"""
    _attribute_vector(g::BaseGraph, key::AbstractString)::Vector{Float64}

Internal cached dense attribute column (do not mutate). Materializes on first use.
"""
function _attribute_vector(g::BaseGraph, key::AbstractString)::Vector{Float64}
    k = String(key)
    cache = g._attr_cache
    if haskey(cache, k)
        return cache[k]
    end
    col = _materialize_attribute_column(g, k)
    cache[k] = col
    return col
end

"""
    attribute_vector(g::BaseGraph, key::AbstractString)::Vector{Float64}

Return a **copy** of the dense attribute column for `key` (safe to own; do not
rely on mutating it to update the graph — use `set_attribute!` instead).
"""
function attribute_vector(g::BaseGraph, key::AbstractString)::Vector{Float64}
    return copy(_attribute_vector(g, key))
end

"""
    prefetch_attribute!(g::BaseGraph, key::AbstractString)

Force materialization of the attribute cache for `key`.
"""
function prefetch_attribute!(g::BaseGraph, key::AbstractString)
    _attribute_vector(g, key)
    return g
end

"""
    set_attribute!(g::BaseGraph, node::Int, key::AbstractString, value)

Update a single node attribute and invalidate the cached column for `key`.
"""
function set_attribute!(g::BaseGraph, node::Int, key::AbstractString, value)
    k = String(key)
    g.attributes[node][k] = value
    delete!(g._attr_cache, k)
    return g
end

"""
    set_attributes!(g::BaseGraph, key::AbstractString, values::AbstractVector)

Replace attribute `key` on all nodes and refresh the cache entry.
"""
function set_attributes!(g::BaseGraph, key::AbstractString, values::AbstractVector)
    length(values) == g.num_nodes || throw(
        ArgumentError(
            "values length ($(length(values))) must equal num_nodes ($(g.num_nodes))",
        ),
    )
    k = String(key)
    @inbounds for i = 1:g.num_nodes
        g.attributes[i][k] = values[i]
    end
    delete!(g._attr_cache, k)
    if values isa AbstractVector{<:Number}
        g._attr_cache[k] = Float64.(values)
    end
    return g
end
