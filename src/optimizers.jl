"""
    greedy_accept(score_fn::Function; maximize::Bool=true) -> Function

Returns an acceptance function `(graph, current_state, candidate) -> Float64` that
greedily accepts candidates that improve or equal the `score_fn` value.
"""
function greedy_accept(score_fn::Function; maximize::Bool=true)
    return (graph, current_state, candidate) -> begin
        curr_val = score_fn(graph, current_state)
        cand_val = score_fn(graph, candidate)
        improved = maximize ? (cand_val >= curr_val) : (cand_val <= curr_val)
        return improved ? 1.0 : 0.0
    end
end

"""
    simulated_annealing_accept(score_fn::Function, temp::Float64; maximize::Bool=true) -> Function

Returns an acceptance function implementing the Metropolis criterion for simulated
annealing.
"""
function simulated_annealing_accept(score_fn::Function, temp::Float64; maximize::Bool=true)
    return (graph, current_state, candidate) -> begin
        curr_val = score_fn(graph, current_state)
        cand_val = score_fn(graph, candidate)
        diff = maximize ? (cand_val - curr_val) : (curr_val - cand_val)
        if diff >= 0
            return 1.0
        else
            return exp(diff / temp)
        end
    end
end
