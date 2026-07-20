abstract type AbstractProposal end

struct RecomProposal <: AbstractProposal
    D₁::Int
    D₂::Int
    D₁_pop::Int
    D₂_pop::Int
    D₁_nodes::BitSet
    D₂_nodes::BitSet
end

struct FlipProposal <: AbstractProposal
    node::Int  # Node that is being flipped.
    D₁::Int  # Original district.
    D₂::Int  # New district.
    D₁_pop::Int
    D₂_pop::Int
    D₁_nodes::BitSet
    D₂_nodes::BitSet
end

struct DummyProposal <: AbstractProposal
    reason::String
end
