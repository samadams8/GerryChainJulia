module GerryChain
using JSON
using HDF5
using SparseArrays
using LightGraphs
using Random
using DataStructures
using Statistics
using DelimitedFiles
using PyPlot
using ProgressBars
using Shapefile: Shapefile
using LibGEOS: LibGEOS
using LibSpatialIndex: LibSpatialIndex
using Logging

export BaseGraph,
    AbstractPartition,
    Partition,
    get_attributes,
    get_district_nodes,
    get_district_populations,
    get_district_adj_and_cut_edges,
    get_subgraph_population,
    induced_subgraph_edges,
    update_partition_adjacency,
    clone_for_update,
    PartitionBuffers,

    # graph / partition accessors
    # (AbstractGraph / neighbors / kruskal_mst are intentionally not exported —
    # they clash with LightGraphs; use GerryChain.AbstractGraph, field access or
    # GerryChain.neighbors / GerryChain.kruskal_mst instead.)
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
    set_edge_penalties_from_pairs!,
    add_region_column!,
    attribute_vector,
    prefetch_attribute!,
    set_attribute!,
    set_attributes!,

    # balance edges
    configure_mst_weights!,
    build_mst_weights!,
    wilson_ust,
    MSTScratch,
    SubtreeCutScratch,
    kruskal_mst!,

    # proposals
    AbstractProposal,
    RecomProposal,
    FlipProposal,
    DummyProposal,

    # constraints
    AbstractConstraint,
    PopulationConstraint,
    ContiguityConstraint,
    satisfy_constraint,

    # recom
    update_partition!,
    recom_chain,
    recom_chain_iter,
    RecomChainIter,
    sample_subgraph,
    get_balanced_proposal,
    get_balanced_proposal_subtree_population,
    get_valid_proposal,

    # flip
    flip_chain,
    flip_chain_iter,
    FlipChainIter,

    # scores
    DistrictAggregate,
    DistrictScore,
    PlanScore,
    CompositeScore,
    AbstractScore,
    ChainScoreData,
    score_initial_partition,
    score_partition_from_proposal,
    eval_score_on_district,
    get_scores_at_step,
    eval_score_on_partition,
    save_scores_to_csv,
    save_scores_to_json,
    save_scores_to_hdf5,
    get_score_values,
    coerce_aggregated_attributes!,

    # acceptance functions
    always_accept,
    satisfies_acceptance_fn,

    # election
    AbstractElection,
    Election,
    ElectionTracker,
    vote_count,
    vote_share,
    seats_won,
    mean_median,
    wasted_votes,
    efficiency_gap,

    # plot
    score_boxplot,
    score_histogram

include("./graph.jl")
include("./partition.jl")
include("./balance_edges.jl")
include("./geo.jl")
include("./proposals.jl")
include("./constraints.jl")
include("./scores.jl")
include("./recom.jl")
include("./flip.jl")
include("./accept.jl")
include("./election.jl")
include("./plot.jl")

end # module
