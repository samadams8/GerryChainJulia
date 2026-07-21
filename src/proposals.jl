abstract type AbstractProposal end

"""
    AbstractProposalConfiguration

Abstract base type for all proposal configurations. Subclasses must implement
`propose(graph::AbstractGraph, partition::Partition, config::YourConfig)`.
"""
abstract type AbstractProposalConfiguration end

"""
    propose(graph::AbstractGraph, partition::Partition, config::AbstractProposalConfiguration) -> Partition

Apply the proposal configured by `config` on `partition` in the context of `graph`.
Returns a new `Partition` instance (typically generated via `clone_for_update` and updated).
"""
function propose end

struct RecomPayload <: AbstractProposal
    D₁::Int
    D₂::Int
    D₁_pop::Int
    D₂_pop::Int
    D₁_nodes::BitSet
    D₂_nodes::BitSet
end

struct FlipPayload <: AbstractProposal
    node::Int  # Node that is being flipped.
    D₁::Int  # Original district.
    D₂::Int  # New district.
    D₁_pop::Int
    D₂_pop::Int
    D₁_nodes::BitSet
    D₂_nodes::BitSet
end
