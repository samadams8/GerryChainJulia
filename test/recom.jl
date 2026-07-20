@testset "Recom tests" begin
    graph = BaseGraph(square_grid_filepath, "population")

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
        component = GerryChain.traverse_mst(mst, 2, 3, stack, component_container)
        @test component == BitSet([1, 2, 5, 6])
        component = GerryChain.traverse_mst(mst, 1, 5, stack, component_container)
        @test component == BitSet([1])
    end

    @testset "ReCom proposal generation with tolerance" begin
        partition = Partition(graph, "assignment")

        ideal_pop = total_pop(graph) / num_dists(partition)
        config = ReComConfiguration(ideal_pop, 0.1; rng=MersenneTwister(42))
        prop = propose(graph, partition, config)
        @test prop isa Partition
        @test num_dists(prop) == 4

        # Test MarkovChain integration
        mc = MarkovChain(
            graph,
            config,
            [(g, p) -> within_percent_of_ideal_population(g, p, 0.1)],
            always_accept,
            partition,
            3
        )

        steps = 0
        for state in mc
            steps += 1
            @test state isa Partition
        end
        @test steps == 3
    end

    @testset "get_balanced_proposal_subtree_population" begin
        partition = Partition(graph, "assignment")
        ideal_pop = total_pop(graph) / num_dists(partition)
        min_pop = Int(ceil(0.1 * ideal_pop))
        max_pop = Int(floor(2.0 * ideal_pop))

        D₁, D₂, sg_edges, sg_nodes = GerryChain.sample_subgraph(graph, partition, MersenneTwister(42))
        mst_edges = GerryChain._kruskal_mst(
            graph, sg_edges, collect(sg_nodes), MersenneTwister(7)
        )

        edge_scan = GerryChain.get_balanced_proposal(
            graph, mst_edges, sg_nodes, partition, min_pop, max_pop, D₁, D₂
        )
        subtree = GerryChain.get_balanced_proposal_subtree_population(
            graph, mst_edges, sg_nodes, partition, min_pop, max_pop, D₁, D₂
        )

        @test edge_scan isa Union{RecomProposal,GerryChain.DummyProposal}
        @test subtree isa Union{RecomProposal,GerryChain.DummyProposal}

        if edge_scan isa RecomProposal
            @test subtree isa RecomProposal
            @test subtree.D₁_nodes ∪ subtree.D₂_nodes == sg_nodes
            @test isempty(subtree.D₁_nodes ∩ subtree.D₂_nodes)
            @test subtree.D₁_pop == edge_scan.D₁_pop
            @test subtree.D₂_pop == edge_scan.D₂_pop
        else
            @test subtree isa GerryChain.DummyProposal
        end
    end
end
