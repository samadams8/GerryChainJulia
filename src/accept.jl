"""
    always_accept(args...) -> Float64

Acceptance function that accepts every proposal with probability 1.
Compatible with the `MarkovChain` accept signature `(graph, current_state, candidate)`.
"""
always_accept(args...; kwargs...) = 1.0

"""
    satisfies_acceptance_fn(partition::AbstractPartition,
                            acceptance_fn::Function,
                            rng::AbstractRNG=Random.default_rng()) -> Bool

Determines whether a partition should be accepted according to a user-provided
acceptance function. The acceptance function must return a probability in [0, 1];
this function samples from that probability to decide acceptance.

The acceptance function should accept a single `AbstractPartition` argument and
compare it to `partition.parent` if needed.
"""
function satisfies_acceptance_fn(
    partition::AbstractPartition,
    acceptance_fn::Function,
    rng::AbstractRNG = Random.default_rng(),
)::Bool
    @assert parent(partition) !== nothing "partition must have a valid parent"
    prob = acceptance_fn(partition)
    if !(prob isa Number) || prob < 0 || prob > 1
        throw(ArgumentError("Acceptance function must return value in [0, 1] range."))
    end
    return rand(rng) <= prob
end

