abstract type AbstractChain end

"""
    MarkovChain{G,S,P,C} <: AbstractChain

A mutable iterator over redistricting partitions driven by proposal strategy and constraints.
Follows the standard Julia iterator protocol.

# Fields
- `graph`                   – the underlying graph (e.g. `GerryChain.BaseGraph`)
- `proposal_config`         – `AbstractProposalConfiguration` object specifying the proposal strategy
- `constraints`             – `Tuple` of constraints (e.g. `PopulationConstraint`)
- `state`                   – current `GerryChain.Partition` of type `S`
- `step`                    – 0-based counter of steps taken
- `total_steps`             – maximum number of steps
- `max_constraint_attempts` – maximum constraint retry attempts before error (0 = unlimited)
- `show_progress`           – whether to display a ProgressMeter bar
"""
mutable struct MarkovChain{G,S,P,C} <: AbstractChain
    graph::G
    proposal_config::P
    constraints::C
    state::S
    step::Int
    total_steps::Int
    max_constraint_attempts::Int
    show_progress::Bool
end

_eval_constraint(c::AbstractConstraint, g, p) = satisfies_constraint(c, g, p)
_eval_constraint(c, g, p) = c(g, p)

function MarkovChain(
    graph::G,
    proposal_config::P,
    constraints,
    initial_state::S,
    total_steps::Int;
    show_progress::Bool = false,
    max_constraint_attempts::Int = 0,
) where {G,S,P}
    c_tuple = constraints isa Tuple ? constraints : (constraints isa AbstractVector ? Tuple(constraints) : (constraints,))
    return MarkovChain(
        graph,
        proposal_config,
        c_tuple,
        initial_state,
        0,
        total_steps,
        max_constraint_attempts,
        show_progress,
    )
end

Base.length(chain::MarkovChain) = chain.total_steps
Base.eltype(::Type{<:MarkovChain{G,S}}) where {G,S} = S

function Base.iterate(chain::MarkovChain{G,S}, progress = nothing) where {G,S}
    if chain.step == 0 && chain.show_progress
        progress = ProgressBar(1:chain.total_steps)
        set_description(progress, "Chain: ")
    end

    chain.step >= chain.total_steps && return nothing

    local candidate::S
    attempts = 0
    while true
        candidate = propose(chain.graph, chain.state, chain.proposal_config)::S
        all(_eval_constraint(c, chain.graph, candidate) for c in chain.constraints) && break
        attempts += 1
        if chain.max_constraint_attempts > 0 && attempts >= chain.max_constraint_attempts
            error(
                "MarkovChain: exceeded max_constraint_attempts=$(chain.max_constraint_attempts) without finding a valid proposal",
            )
        end
    end

    chain.state = candidate
    if chain.state.parent !== nothing
        chain.state.parent.parent = nothing
    end

    chain.step += 1
    chain.show_progress && update(progress)
    return chain.state, progress
end

# CouponCollectorChain implementation.

"""
    coupon_collector_expectation(n::Int) -> Float64

Expected number of draws to collect all `n` unique coupons (harmonic series formula).
"""
coupon_collector_expectation(n::Int) = n * sum(1.0 / i for i in 1:n)

"""
    CouponCollectorChain{G,S,P,C} <: AbstractChain

Wraps a `MarkovChain{G,S,P,C}` and burns `micro_steps_per_yield` micro-steps internally
for each macro-step it yields, decorrelating the sampled partitions.
"""
mutable struct CouponCollectorChain{G,S,P,C} <: AbstractChain
    chain::MarkovChain{G,S,P,C}
    micro_steps_per_yield::Int
    num_macro_steps::Int
    macro_step::Int
    show_progress::Bool
end

function CouponCollectorChain(
    graph,
    proposal_config,
    constraints,
    initial_state::S,
    num_macro_steps::Int,
    micro_steps_per_yield::Int;
    show_progress::Bool = false,
    max_constraint_attempts::Int = 0,
) where {S}
    inner = MarkovChain(
        graph,
        proposal_config,
        constraints,
        initial_state,
        num_macro_steps * micro_steps_per_yield;
        show_progress = false,
        max_constraint_attempts = max_constraint_attempts,
    )
    return CouponCollectorChain(
        inner, micro_steps_per_yield, num_macro_steps, 0, show_progress
    )
end

Base.length(c::CouponCollectorChain) = c.num_macro_steps
Base.eltype(::Type{<:CouponCollectorChain{G,S}}) where {G,S} = S

function Base.iterate(c::CouponCollectorChain{G,S}, progress = nothing) where {G,S}
    if c.macro_step == 0 && c.show_progress
        progress = ProgressBar(1:c.num_macro_steps)
        set_description(progress, "CouponCollector: ")
    end

    c.macro_step >= c.num_macro_steps && return nothing

    local state::S
    for _ in 1:c.micro_steps_per_yield
        result = iterate(c.chain, progress)
        result === nothing && return nothing
        state = result[1]::S
    end

    c.macro_step += 1
    c.show_progress && update(progress)
    return state, progress
end
