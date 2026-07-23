@testset "ShortBurstChain" begin
    g = SimpleGraph(4)
    add_edge!(g, 1, 2)
    add_edge!(g, 2, 3)
    add_edge!(g, 3, 4)

    pops = [10, 10, 10, 10]
    adj_mat = spzeros(Int, 4, 4)
    adj_mat[1, 2] = adj_mat[2, 1] = 1
    adj_mat[2, 3] = adj_mat[3, 2] = 2
    adj_mat[3, 4] = adj_mat[4, 3] = 3

    edge_src = [1, 2, 3]
    edge_dst = [2, 3, 4]
    nbrs = [[2], [1, 3], [2, 4], [3]]
    attrs = [Dict{String,Any}("POP" => 10) for _ in 1:4]

    graph = BaseGraph(
        4, 3, 40, pops, adj_mat, edge_src, edge_dst, nbrs, g, attrs
    )

    dist_nodes = [BitSet([1, 2]), BitSet([3, 4])]
    dist_pops = [20, 20]
    assignments = [1, 1, 2, 2]
    dist_adj = spzeros(Int, 2, 2)
    dist_adj[1, 2] = dist_adj[2, 1] = 1
    cut_edges = [0, 1, 0]

    initial_partition = Partition(
        2, 1, assignments, dist_pops, cut_edges, dist_adj, dist_nodes, nothing
    )

    mutable struct TestProposalConfig <: GerryChain.AbstractProposalConfiguration
        states::Vector{Partition}
        index::Int
    end
    
    function GerryChain.propose(g, p, pc::TestProposalConfig)
        pc.index += 1
        return pc.states[pc.index]
    end

    state1 = Partition(2, 1, [1, 1, 2, 2], dist_pops, cut_edges, dist_adj, dist_nodes, initial_partition)
    state2 = Partition(2, 1, [1, 1, 2, 2], dist_pops, cut_edges, dist_adj, dist_nodes, state1)
    state3 = Partition(2, 1, [1, 2, 1, 2], dist_pops, cut_edges, dist_adj, dist_nodes, state2)
    state4 = Partition(2, 1, [1, 1, 2, 2], dist_pops, cut_edges, dist_adj, dist_nodes, state3)
    state5 = Partition(2, 1, [1, 1, 1, 2], dist_pops, cut_edges, dist_adj, dist_nodes, state4)
    state6 = Partition(2, 1, [1, 1, 2, 2], dist_pops, cut_edges, dist_adj, dist_nodes, state5)

    states_vec = [state1, state2, state3, state4, state5, state6]
    
    scores = Dict{Partition, Float64}(
        initial_partition => 2.0,
        state1 => 5.0,
        state2 => 3.0,
        state3 => 6.0,
        state4 => 4.0,
        state5 => 6.0,
        state6 => 2.0
    )
    score_fn = (g, p) -> scores[p]

    pc = TestProposalConfig(states_vec, 0)
    true_constraint = (g, p) -> true

    sbc = ShortBurstChain(
        graph,
        pc,
        (true_constraint,),
        initial_partition,
        2, # num_bursts
        3, # burst_length
        score_fn;
        maximize = true
    )

    @test length(sbc) == 6
    @test eltype(typeof(sbc)) == Partition

    yielded = collect(sbc)
    @test length(yielded) == 6
    @test yielded[1] === state1
    @test yielded[2] === state2
    @test yielded[3] === state3
    @test yielded[4] === state4
    @test yielded[5] === state5
    @test yielded[6] === state6

    # Verify best score and best state
    @test sbc.best_score == 6.0
    @test sbc.best_state === state5

    # Test minimizing:
    pc_min = TestProposalConfig(states_vec, 0)
    sbc_min = ShortBurstChain(
        graph,
        pc_min,
        (true_constraint,),
        initial_partition,
        2,
        3,
        score_fn;
        maximize = false
    )
    
    yielded_min = collect(sbc_min)
    @test sbc_min.best_score == 2.0
    @test sbc_min.best_state === state6
end
