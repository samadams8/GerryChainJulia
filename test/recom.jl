@testset "Recom tests" begin
    graph = BaseGraph(square_grid_filepath, "population")

    function accept_on_third_try()
        counter = 0
        return p -> (counter += 1) < 3 ? 0.0 : 1.0
    end

    @testset "traverse_mst()" begin
        nodes = [1, 2, 3, 4, 5, 6, 7, 8]
        edges = [
            graph.adj_matrix[1, 5],
            graph.adj_matrix[5, 6],
            graph.adj_matrix[2, 6],
            graph.adj_matrix[2, 3],
            graph.adj_matrix[3, 7],
            graph.adj_matrix[3, 4],
            graph.adj_matrix[4, 8],
        ]
        mst = GerryChain.build_mst(graph, BitSet(nodes), BitSet(edges))
        stack = Stack{Int}()
        component_container = BitSet([])
        cut_edge = graph.adj_matrix[2, 3]
        component = GerryChain.traverse_mst(mst, 2, 3, stack, component_container)
        @test component == BitSet([1, 2, 5, 6])
        component = GerryChain.traverse_mst(mst, 1, 5, stack, component_container)
        @test component == BitSet([1])
    end

    @testset "recom_chain()" begin
        partition = Partition(graph, "assignment")
        # this is a dummy constraint
        pop_constraint = PopulationConstraint(graph, partition, 10.0)
        scores = [
            DistrictAggregate("electionD"),
            DistrictAggregate("electionR"),
            DistrictAggregate("purple"),
            DistrictAggregate("pink"),
        ]
        num_steps = 2 # test 2 steps for now

        function run_chain()
            try
                recom_chain(graph, partition, pop_constraint, num_steps, scores)
            catch ex
                return ex
            end
        end
        recom_chain(graph, partition, pop_constraint, num_steps, scores)
        # hacky way to run flip chain and test that it doesn't yield an exception
        @test !isa(run_chain(), Exception)
    end

    @testset "no_self_loops" begin
        partition = Partition(graph, "assignment")
        # this is a dummy constraint
        pop_constraint = PopulationConstraint(graph, partition, 10.0)
        scores = [
            DistrictAggregate("electionD"),
            DistrictAggregate("electionR"),
            DistrictAggregate("purple"),
            DistrictAggregate("pink"),
        ]
        num_steps = 1 # test 1 step (2 states) for now
        f = accept_on_third_try()
        chain_data = recom_chain(
            graph, partition, pop_constraint, num_steps, scores, acceptance_fn=f
        )
        @test get_scores_at_step(chain_data, 0) == get_scores_at_step(chain_data, 1)
        # acceptance function should still return 0, because the acceptance function
        # should only have been called once
        @test f(nothing) == 0.0

        f = accept_on_third_try() # reset f
        chain_data = recom_chain(
            graph,
            partition,
            pop_constraint,
            num_steps,
            scores,
            acceptance_fn=f,
            no_self_loops=true,
        )
        # acceptance function should now return 1, because the acceptance function
        # should have been called until it started returning 1.0
        @test f(nothing) == 1.0
    end

    @testset "region-aware recom_chain kwargs" begin
        graph = BaseGraph(square_grid_filepath, "population"; region_columns=["assignment"])
        partition = Partition(graph, "assignment")
        pop_constraint = PopulationConstraint(graph, partition, 10.0)
        scores = [DistrictAggregate("purple")]
        chain_data = recom_chain(
            graph,
            partition,
            pop_constraint,
            1,
            scores;
            region_surcharges=Dict("assignment" => 0.5),
            progress_bar=false,
        )
        @test length(chain_data.step_values) == 2  # initial + 1 step
    end

    @testset "wilson tree_method" begin
        partition = Partition(graph, "assignment")
        pop_constraint = PopulationConstraint(graph, partition, 10.0)
        scores = [DistrictAggregate("purple")]
        chain_data = recom_chain(
            graph,
            partition,
            pop_constraint,
            1,
            scores;
            tree_method=:wilson,
            progress_bar=false,
            rng=MersenneTwister(7),
        )
        @test length(chain_data.step_values) == 2
    end

    @testset "n_parallel proposals" begin
        partition = Partition(graph, "assignment")
        pop_constraint = PopulationConstraint(graph, partition, 10.0)
        scores = [DistrictAggregate("purple")]
        chain_data = recom_chain(
            graph,
            partition,
            pop_constraint,
            2,
            scores;
            n_parallel=4,
            progress_bar=false,
            rng=MersenneTwister(11),
        )
        @test length(chain_data.step_values) == 3

        # n_parallel=1 is seed-stable
        p1 = Partition(graph, "assignment")
        p2 = Partition(graph, "assignment")
        s = [DistrictAggregate("purple")]
        c1 = recom_chain(
            graph,
            p1,
            pop_constraint,
            2,
            s;
            n_parallel=1,
            progress_bar=false,
            rng=MersenneTwister(42),
        )
        c2 = recom_chain(
            graph,
            p2,
            pop_constraint,
            2,
            s;
            n_parallel=1,
            progress_bar=false,
            rng=MersenneTwister(42),
        )
        @test get_score_values(c1, "purple") == get_score_values(c2, "purple")
    end

    @testset "abstract-typed chain entrypoints" begin
        g::GerryChain.AbstractGraph = BaseGraph(square_grid_filepath, "population")
        p::GerryChain.AbstractPartition = Partition(g, "assignment")
        pop_c = PopulationConstraint(g, p, 10.0)
        scores = [DistrictAggregate("purple")]
        data = recom_chain(
            g, p, pop_c, 1, scores; progress_bar=false, rng=MersenneTwister(3)
        )
        @test length(data.step_values) == 2
    end

    @testset "get_balanced_proposal_subtree_population" begin
        partition = Partition(graph, "assignment")
        pop_constraint = PopulationConstraint(graph, partition, 10.0)

        D₁, D₂, sg_edges, sg_nodes = sample_subgraph(graph, partition, MersenneTwister(42))
        mst_edges = random_kruskal_mst(
            graph, sg_edges, collect(sg_nodes), MersenneTwister(7)
        )

        edge_scan = get_balanced_proposal(
            graph, mst_edges, sg_nodes, partition, pop_constraint, D₁, D₂
        )
        subtree = get_balanced_proposal_subtree_population(
            graph, mst_edges, sg_nodes, partition, pop_constraint, D₁, D₂
        )

        @test edge_scan isa Union{RecomProposal,DummyProposal}
        @test subtree isa Union{RecomProposal,DummyProposal}

        if edge_scan isa RecomProposal
            @test subtree isa RecomProposal
            @test satisfy_constraint(pop_constraint, subtree.D₁_pop, subtree.D₂_pop)
            @test subtree.D₁_nodes ∪ subtree.D₂_nodes == sg_nodes
            @test isempty(subtree.D₁_nodes ∩ subtree.D₂_nodes)
            @test subtree.D₁_pop == edge_scan.D₁_pop
            @test subtree.D₂_pop == edge_scan.D₂_pop
            @test subtree.D₁_nodes == edge_scan.D₁_nodes
            @test subtree.D₂_nodes == edge_scan.D₂_nodes
        else
            @test subtree isa DummyProposal
        end

        chain_data = recom_chain(
            graph,
            Partition(graph, "assignment"),
            pop_constraint,
            1,
            [DistrictAggregate("purple")];
            cut_method=:subtree_population,
            progress_bar=false,
            rng=MersenneTwister(19),
        )
        @test length(chain_data.step_values) == 2

        chain_edge = recom_chain(
            graph,
            Partition(graph, "assignment"),
            pop_constraint,
            1,
            [DistrictAggregate("purple")];
            cut_method=:edge_scan,
            progress_bar=false,
            rng=MersenneTwister(19),
        )
        @test length(chain_edge.step_values) == 2

        @test_throws ArgumentError get_valid_proposal(
            graph,
            Partition(graph, "assignment"),
            pop_constraint,
            MersenneTwister(1),
            1;
            cut_method=:bogus,
        )

        # SubtreeCutScratch reuse: same result as fresh call; reusable across MSTs
        cut_scratch = SubtreeCutScratch()
        fresh = get_balanced_proposal_subtree_population(
            graph, mst_edges, sg_nodes, partition, pop_constraint, D₁, D₂
        )
        reused1 = get_balanced_proposal_subtree_population(
            graph,
            mst_edges,
            sg_nodes,
            partition,
            pop_constraint,
            D₁,
            D₂;
            scratch=cut_scratch,
        )
        reused2 = get_balanced_proposal_subtree_population(
            graph,
            mst_edges,
            sg_nodes,
            partition,
            pop_constraint,
            D₁,
            D₂;
            scratch=cut_scratch,
        )
        @test typeof(fresh) == typeof(reused1) == typeof(reused2)
        if fresh isa RecomProposal
            @test reused1 isa RecomProposal
            @test reused2 isa RecomProposal
            @test reused1.D₁_pop == fresh.D₁_pop == reused2.D₁_pop
            @test reused1.D₂_pop == fresh.D₂_pop == reused2.D₂_pop
            @test reused1.D₁_nodes == fresh.D₁_nodes == reused2.D₁_nodes
            @test reused1.D₂_nodes == fresh.D₂_nodes == reused2.D₂_nodes
        else
            @test reused1 isa DummyProposal
            @test reused2 isa DummyProposal
        end

        mst2 = random_kruskal_mst(graph, sg_edges, collect(sg_nodes), MersenneTwister(99))
        other = get_balanced_proposal_subtree_population(
            graph, mst2, sg_nodes, partition, pop_constraint, D₁, D₂; scratch=cut_scratch
        )
        @test other isa Union{RecomProposal,DummyProposal}
    end

    @testset "recom_chain_iter" begin
        partition = Partition(graph, "assignment")
        pop_constraint = PopulationConstraint(graph, partition, 10.0)
        scores = [DistrictAggregate("purple")]
        num_steps = 3

        iter = recom_chain_iter(
            graph,
            partition,
            pop_constraint,
            num_steps,
            scores;
            progress_bar=false,
            rng=MersenneTwister(42),
        )

        @test iter isa RecomChainIter
        @test length(iter) == num_steps
        @test eltype(iter) == Tuple{Partition,Dict{String,Any}}

        # Iterate step-by-step
        results = collect(iter)
        @test length(results) == num_steps
        for (p, s) in results
            @test p isa Partition
            @test s isa Dict{String,Any}
            @test haskey(s, "purple")
            @test haskey(s, "dists")
        end

        # Test with progress_bar = true
        iter_pb = recom_chain_iter(
            graph,
            partition,
            pop_constraint,
            num_steps,
            scores;
            progress_bar=true,
            rng=MersenneTwister(42),
        )
        @test iter_pb isa ProgressBar
        results_pb = collect(iter_pb)
        @test length(results_pb) == num_steps
    end
end
