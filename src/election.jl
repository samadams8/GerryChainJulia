abstract type AbstractElection end

"""
    Election <: AbstractElection

Represents an election scenario with specified party column keys.
"""
struct Election <: AbstractElection
    name::String
    parties::Vector{String}
end

Election(parties::Vector{String}) = Election("Election", parties)

"""
    vote_count(graph::AbstractGraph, partition::AbstractPartition, party::String) -> Vector{Float64}

Computes the total votes for `party` in each district of `partition`.
"""
function vote_count(graph::AbstractGraph, partition::AbstractPartition, party::String)::Vector{Float64}
    n_dists = num_dists(partition)
    d_nodes = dist_nodes(partition)
    col = _attribute_vector(graph, party)
    counts = zeros(Float64, n_dists)
    for d in 1:n_dists
        s = 0.0
        @inbounds for node in d_nodes[d]
            s += col[node]
        end
        counts[d] = s
    end
    return counts
end

"""
    vote_share(graph::AbstractGraph, partition::AbstractPartition, party::String, total_pop_col::String) -> Vector{Float64}

Computes the vote share of `party` in each district of `partition`.
"""
function vote_share(
    graph::AbstractGraph, partition::AbstractPartition, party::String, total_pop_col::String
)::Vector{Float64}
    counts = vote_count(graph, partition, party)
    totals = vote_count(graph, partition, total_pop_col)
    return counts ./ totals
end

"""
    seats_won(party_votes::Vector{Float64}, other_votes::Vector{Float64}) -> Int

Computes the number of districts won by `party_votes` against `other_votes`.
Ties count as 0 seats won for both parties.
"""
function seats_won(party_votes::Vector{Float64}, other_votes::Vector{Float64})::Int
    seats = 0
    for i in 1:length(party_votes)
        if party_votes[i] > other_votes[i]
            seats += 1
        end
    end
    return seats
end

"""
    mean_median(vote_shares::Vector{Float64}) -> Float64

Computes the mean-median partisan bias metric for a party's district vote shares.
"""
function mean_median(vote_shares::Vector{Float64})::Float64
    return median(vote_shares) - mean(vote_shares)
end

"""
    wasted_votes(party1_votes::Float64, party2_votes::Float64) -> Tuple{Float64, Float64}

Computes the number of votes wasted by party 1 and party 2 in a single district.
"""
function wasted_votes(party1_votes::Float64, party2_votes::Float64)::Tuple{Float64,Float64}
    total = party1_votes + party2_votes
    if party1_votes > party2_votes
        w1 = party1_votes - total / 2
        w2 = party2_votes
    elseif party2_votes > party1_votes
        w1 = party1_votes
        w2 = party2_votes - total / 2
    else
        w1, w2 = party1_votes, party2_votes
    end
    return w1, w2
end

"""
    efficiency_gap(party1_district_votes::Vector{Float64}, party2_district_votes::Vector{Float64}) -> Float64

Computes the efficiency gap metric for party 1.
"""
function efficiency_gap(party1_district_votes::Vector{Float64}, party2_district_votes::Vector{Float64})::Float64
    w1_total = 0.0
    w2_total = 0.0
    total_votes = 0.0

    for i in 1:length(party1_district_votes)
        p1 = party1_district_votes[i]
        p2 = party2_district_votes[i]
        w1, w2 = wasted_votes(p1, p2)
        w1_total += w1
        w2_total += w2
        total_votes += (p1 + p2)
    end

    return (w1_total - w2_total) / total_votes
end
