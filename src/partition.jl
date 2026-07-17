abstract type AbstractPartition end

mutable struct Partition <: AbstractPartition
    num_dists::Int
    num_cut_edges::Int
    assignments::Array{Int,1}                  # of length(num_nodes)
    dist_populations::Array{Int,1}             # of length(num_districts)
    cut_edges::Array{Int,1}                    # of length(num_edges)
    dist_adj::SparseMatrixCSC{Int,Int}
    dist_nodes::Array{BitSet}
    parent::Union{AbstractPartition,Nothing}   # optional parent partition
end

# --- AbstractPartition accessors (default concrete implementations) ---

num_dists(p::Partition) = p.num_dists
num_cut_edges(p::Partition) = p.num_cut_edges
assignments(p::Partition) = p.assignments
dist_populations(p::Partition) = p.dist_populations
cut_edges(p::Partition) = p.cut_edges
dist_adj(p::Partition) = p.dist_adj
dist_nodes(p::Partition) = p.dist_nodes
Base.parent(p::Partition) = p.parent

"""
    _copy_partition_fields(p::Partition; parent=nothing)::Partition

Copy partition arrays; **share** `dist_nodes` BitSet references (callers must
`copy` a BitSet before mutating it). Used by `clone_for_update` and parent snapshots.
"""
function _copy_partition_fields(
    p::Partition;
    parent::Union{AbstractPartition,Nothing} = nothing,
)::Partition
    return Partition(
        p.num_dists,
        p.num_cut_edges,
        copy(p.assignments),
        copy(p.dist_populations),
        copy(p.cut_edges),
        copy(p.dist_adj),
        BitSet[p.dist_nodes[i] for i = 1:p.num_dists],  # share BitSets
        parent,
    )
end

"""
    PartitionBuffers

Optional reusable arrays for `clone_for_update`. After a clone, the buffers
are refreshed so the returned `Partition` owns the previous vectors.
"""
mutable struct PartitionBuffers
    assignments::Vector{Int}
    dist_populations::Vector{Int}
    cut_edges::Vector{Int}
end

function PartitionBuffers(p::Partition)
    return PartitionBuffers(
        similar(p.assignments),
        similar(p.dist_populations),
        similar(p.cut_edges),
    )
end

"""
    clone_for_update(p::AbstractPartition)::AbstractPartition

Return a clone suitable for mutation. Arrays are copied; district `BitSet`s are
shared until an update copies the districts it mutates. Graph is not involved.
"""
function clone_for_update(p::Partition)::Partition
    return _copy_partition_fields(p; parent = p)
end

function clone_for_update(p::Partition, buffers::PartitionBuffers)::Partition
    copyto!(buffers.assignments, p.assignments)
    copyto!(buffers.dist_populations, p.dist_populations)
    copyto!(buffers.cut_edges, p.cut_edges)
    new_p = Partition(
        p.num_dists,
        p.num_cut_edges,
        buffers.assignments,
        buffers.dist_populations,
        buffers.cut_edges,
        copy(p.dist_adj),
        BitSet[p.dist_nodes[i] for i = 1:p.num_dists],
        p,
    )
    buffers.assignments = similar(p.assignments)
    buffers.dist_populations = similar(p.dist_populations)
    buffers.cut_edges = similar(p.cut_edges)
    return new_p
end

"""
    Partition(graph::BaseGraph,
              assignment_col::AbstractString)::Partition

Partition represents a partition of the nodes of the graph.
It contains plan-specific information that will change each time we
change our plan.

*Arguments:*
- graph:          BaseGraph object that has the underlying network
                  structure of the plan.
- assignment_col: the key denoting the district assignment at the
                  node level
"""
function Partition(graph::BaseGraph, assignment_col::AbstractString)::Partition
    populations = graph.populations
    assignments = get_assignments(graph.attributes, assignment_col)
    num_districts = length(Set(assignments))

    # get cut_edges, district_adjacencies
    dist_adj, cut_edges = get_district_adj_and_cut_edges(graph, assignments, num_districts)
    num_cut_edges = sum(cut_edges)

    # get district populations
    dist_populations =
        get_district_populations(assignments, populations, graph.num_nodes, num_districts)

    # get district_nodes
    dist_nodes = get_district_nodes(assignments, graph.num_nodes, num_districts)

    # return Partition with no parent by default
    return Partition(
        num_districts,
        num_cut_edges,
        assignments,
        dist_populations,
        cut_edges,
        dist_adj,
        dist_nodes,
        nothing,
    )
end

"""
    get_assignments(node_attributes::Array,
                    assignment_col::AbstractString)::Array{Int, 1}

Extracts the assignments for each node in a graph. First attempts
to parse the values of the assignments column as an Int; if it fails,
it assumes every unique string assignment corresponds to to a unique
district.
*Returns* an Int Array of length `length(node_attributes)`.
"""
function get_assignments(
    node_attributes::Array,
    assignment_col::AbstractString,
)::Array{Int,1}
    assignment_to_num = Dict{String,Int}() # map unique strings to integers
    raw_assignments = get_attribute_by_key(node_attributes, assignment_col)
    processed_assignments = zeros(length(raw_assignments))
    for (i, raw_value) in enumerate(raw_assignments)
        if raw_value isa Int
            processed_assignments[i] = raw_value
        elseif raw_value isa String
            try
                processed_assignments[i] = parse(Int, raw_value)
            catch exception # if the String could not be read as an int
                if !haskey(assignment_to_num, raw_value)
                    assignment_to_num[raw_value] = length(assignment_to_num) + 1
                end
                processed_assignments[i] = assignment_to_num[raw_value]
            end
        else
            message = """
                District assignments should be Ints or Strings; type
                $(typeof(raw_value)) was found instead.
            """
            throw(DomainError(raw_value, message))
        end
    end
    return processed_assignments
end

"""
    get_district_nodes(assignments::Array{Int, 1},
                       num_nodes::Int,
                       num_districts::Int)::Array{Set{Int}, 1}

*Returns* an Array of Sets `district_nodes` where the nodes of the i'th
district will be at `district_nodes[i]` as a Set.
"""
function get_district_nodes(
    assignments::Array{Int,1},
    num_nodes::Int,
    num_districts::Int,
)::Array{Set{Int},1}
    district_nodes = [BitSet([]) for _ = 1:num_districts]
    for i = 1:num_nodes
        push!(district_nodes[assignments[i]], i)
    end
    return district_nodes
end

"""
    get_district_populations(assignments::Array{Int, 1},
                             populations::Array{Int, 1},
                             num_nodes::Int,
                             num_districts::Int)::Array{Int, 1}

*Returns* an Array of populations `dist_pops` where the population of
the i'th district is at `dist_pops[i]`.
"""
function get_district_populations(
    assignments::Array{Int,1},
    populations::Array{Int,1},
    num_nodes::Int,
    num_districts::Int,
)::Array{Int,1}
    district_populations = zeros(Int, num_districts)
    for i = 1:num_nodes
        district_populations[assignments[i]] += populations[i]
    end
    return district_populations
end

"""
    get_district_adj_and_cut_edges(graph::BaseGraph,
                                   assignments::Array{Int, 1},
                                   num_districts::Int)

*Returns:*
- district\\_adj: a `num_districts` x `num_districts` matrix where the
                elements are the number of cut edges between the districts
- cut\\_edges:    an Array of size(num\\_edges) where i'th element is 1
                if edge `i` is a cut\\_edge, 0 otherwise
"""
function get_district_adj_and_cut_edges(
    graph::BaseGraph,
    assignments::Array{Int,1},
    num_districts::Int,
)
    district_adj = zeros(Int, num_districts, num_districts)
    num_edges = ne(graph.simple_graph)
    cut_edges = zeros(Int, graph.num_edges)

    for edge in edges(graph.simple_graph)
        src_assignment = assignments[src(edge)]
        dst_assignment = assignments[dst(edge)]

        if src_assignment != dst_assignment
            # look up the index for this edge
            index = graph.adj_matrix[src(edge), dst(edge)]
            cut_edges[index] = 1

            district_adj[src_assignment, dst_assignment] += 1
            district_adj[dst_assignment, src_assignment] += 1
        end
    end
    return district_adj, cut_edges
end

"""
    sample_adjacent_districts_randomly(partition::Partition,
                                       rng::AbstractRNG)

Randomly sample two adjacent districts and return them.
"""
function sample_adjacent_districts_randomly(
    partition::AbstractPartition,
    rng::AbstractRNG,
)
    n = num_dists(partition)
    adj = dist_adj(partition)
    while true
        D₁ = rand(rng, 1:n)
        D₂ = rand(rng, 1:n)
        if adj[D₁, D₂] != 0
            return D₁, D₂
        end
    end
end

"""
    update_partition_adjacency(partition::Partition,
                               graph::BaseGraph)

Updates the district adjacency matrix and cut edges
to reflect the partition's assignments for each node.
"""
function update_partition_adjacency(partition::Partition, graph::AbstractGraph)
    partition.dist_adj = spzeros(Int, partition.num_dists, partition.num_dists)
    partition.num_cut_edges = 0
    srcs = edge_src(graph)
    dsts = edge_dst(graph)
    n_edges = num_edges(graph)

    for i = 1:n_edges
        src_assignment = partition.assignments[srcs[i]]
        dst_assignment = partition.assignments[dsts[i]]

        if src_assignment != dst_assignment
            partition.cut_edges[i] = 1
            partition.num_cut_edges += 1

            partition.dist_adj[src_assignment, dst_assignment] += 1
            partition.dist_adj[dst_assignment, src_assignment] += 1
        else
            partition.cut_edges[i] = 0
        end
    end
end
