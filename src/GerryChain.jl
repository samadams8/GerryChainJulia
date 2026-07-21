module GerryChain

using JSON
using HDF5
using SparseArrays
using LightGraphs
using Random
using DataStructures
using Statistics
using DelimitedFiles
using ProgressBars
using Shapefile: Shapefile
using LibGEOS: LibGEOS
using LibSpatialIndex: LibSpatialIndex
using Logging

export
    # Core abstract types and base structures.
    AbstractGraph,
    BaseGraph,
    AbstractPartition,
    Partition,
    clone_for_update,

    # Graph & Partition accessors.
    num_nodes,
    num_edges,
    total_pop,
    populations,
    edge_src,
    edge_dst,
    edge_penalties,
    has_region,
    region_ids,
    num_dists,
    num_cut_edges,
    assignments,
    dist_populations,
    cut_edges,
    dist_adj,
    dist_nodes,
    set_edge_penalty!,
    add_region_column!,
    attribute_vector,
    set_attribute!,
    set_attributes!,

    # Proposals.
    AbstractProposal,
    AbstractProposalConfiguration,
    ReComConfiguration,
    PopulationFlipConfiguration,
    propose,

    # Constraints.
    AbstractConstraint,
    satisfies_constraint,
    PopulationConstraint,
    ContiguityConstraint,
    within_population_bounds,
    within_percent_of_ideal_population,

    # Election & metrics.
    AbstractElection,
    Election,
    vote_count,
    vote_share,
    seats_won,
    mean_median,
    wasted_votes,
    efficiency_gap,

    # Chain iterators.
    AbstractChain,
    MarkovChain,
    CouponCollectorChain,
    coupon_collector_expectation

include("./graph.jl")
include("./partition.jl")
include("./balance_edges.jl")
include("./geo.jl")
include("./proposals.jl")
include("./constraints.jl")
include("./recom.jl")
include("./flip.jl")
include("./election.jl")
include("./Chain.jl")

end  # module GerryChain
