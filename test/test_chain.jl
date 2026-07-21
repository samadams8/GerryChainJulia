@testset "MarkovChain and CouponCollectorChain" begin
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

    struct DummyProposalConfig <: AbstractProposalConfiguration end
    GerryChain.propose(g, p, ::DummyProposalConfig) = clone_for_update(p)
    dummy_config = DummyProposalConfig()

    true_constraint = (g, p) -> true

    @testset "MarkovChain basic iteration" begin
        mc = MarkovChain(
            graph,
            dummy_config,
            (true_constraint,),
            initial_partition,
            5
        )

        @test length(mc) == 5
        @test eltype(typeof(mc)) == Partition

        steps = 0
        for state in mc
            steps += 1
            @test state isa Partition
            if state.parent !== nothing
                @test state.parent.parent === nothing
            end
        end
        @test steps == 5
    end

    @testset "MarkovChain max_constraint_attempts" begin
        false_constraint = (g, p) -> false
        mc_fail = MarkovChain(
            graph,
            dummy_config,
            (false_constraint,),
            initial_partition,
            5;
            max_constraint_attempts = 10
        )

        @test_throws ErrorException collect(mc_fail)
    end

    @testset "CouponCollectorChain" begin
        @test @inferred(coupon_collector_expectation(1)) == 1.0
        @test @inferred(coupon_collector_expectation(2)) == 3.0

        ccc = CouponCollectorChain(
            graph,
            dummy_config,
            (true_constraint,),
            initial_partition,
            3,
            2
        )

        @test length(ccc) == 3
        @test eltype(typeof(ccc)) == Partition

        macro_steps = 0
        for state in ccc
            macro_steps += 1
            @test state isa Partition
        end
        @test macro_steps == 3
    end
end
