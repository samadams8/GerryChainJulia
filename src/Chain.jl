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

# ShortBurstChain implementation.

"""
    ShortBurstChain{G,S,P,C,F} <: AbstractChain

A chain that runs short bursts of Markov Chain steps, returning every individual micro-step
partition. At the end of each burst of length `burst_length`, the state of the inner chain is
reset/teleported back to the best-known partition encountered so far.

# Fields
- `chain`                  – underlying `MarkovChain`
- `score_fn`               – function `(graph, partition) -> Real` defining the objective score
- `burst_length`           – length of each burst (number of micro-steps)
- `num_bursts`             – total number of bursts to run
- `maximize`               – if `true`, maximize the score; if `false`, minimize it
- `current_burst`          – current burst index (0-based)
- `current_step_in_burst`  – step within the current burst (0-based)
- `best_state`             – the best-scoring partition encountered so far
- `best_score`             – the score of the best-scoring partition encountered so far
- `show_progress`          – whether to show a progress bar
"""
mutable struct ShortBurstChain{G,S,P,C,F} <: AbstractChain
    chain::MarkovChain{G,S,P,C}
    score_fn::F
    burst_length::Int
    num_bursts::Int
    maximize::Bool
    current_burst::Int
    current_step_in_burst::Int
    best_state::S
    best_score::Float64
    show_progress::Bool
end

function ShortBurstChain(
    graph,
    proposal_config,
    constraints,
    initial_state::S,
    num_bursts::Int,
    burst_length::Int,
    score_fn::F;
    maximize::Bool = true,
    show_progress::Bool = false,
    max_constraint_attempts::Int = 0,
) where {S, F}
    inner = MarkovChain(
        graph,
        proposal_config,
        constraints,
        initial_state,
        num_bursts * burst_length;
        show_progress = false,
        max_constraint_attempts = max_constraint_attempts,
    )
    best_score = Float64(score_fn(graph, initial_state))
    return ShortBurstChain(
        inner,
        score_fn,
        burst_length,
        num_bursts,
        maximize,
        0,
        0,
        initial_state,
        best_score,
        show_progress,
    )
end

Base.length(c::ShortBurstChain) = c.num_bursts * c.burst_length
Base.eltype(::Type{<:ShortBurstChain{G,S}}) where {G,S} = S

function Base.iterate(c::ShortBurstChain{G,S}, progress = nothing) where {G,S}
    if c.current_burst == 0 && c.current_step_in_burst == 0 && c.show_progress
        progress = ProgressBar(1:length(c))
        set_description(progress, "ShortBurst: ")
    end

    if c.current_burst >= c.num_bursts
        return nothing
    end

    result = iterate(c.chain, progress)
    result === nothing && return nothing
    state = result[1]::S

    score = Float64(c.score_fn(c.chain.graph, state))
    if c.maximize
        if score >= c.best_score
            c.best_score = score
            c.best_state = state
        end
    else
        if score <= c.best_score
            c.best_score = score
            c.best_state = state
        end
    end

    c.current_step_in_burst += 1
    if c.current_step_in_burst >= c.burst_length
        c.current_burst += 1
        c.current_step_in_burst = 0
        c.chain.state = c.best_state
    end

    c.show_progress && update(progress)
    return state, progress
end
