@testset "Flip tests" begin
    graph = BaseGraph(square_grid_filepath, "population")
    partition = Partition(graph, "assignment")

    @testset "propose_random_flip()" begin
        flip_prop = GerryChain.propose_random_flip(graph, partition)
        pops = partition.dist_populations[[flip_prop.D₁, flip_prop.D₂]]
        @test sum(pops) == flip_prop.D₁_pop + flip_prop.D₂_pop
        @test flip_prop.D₁ != flip_prop.D₂
        neighbors = graph.neighbors[flip_prop.node]
        neighbor_districts = [partition.assignments[n] for n in neighbors]
        @test flip_prop.D₂ in neighbor_districts
    end

    @testset "PopulationFlip proposal generation" begin
        ideal_pop = total_pop(graph) / num_dists(partition)
        config = PopulationFlipConfiguration(ideal_pop, "population"; rng=MersenneTwister(42))
        prop = propose(graph, partition, config)
        @test prop isa Partition
        @test num_dists(prop) == 4

        mc = MarkovChain(
            graph,
            config,
            [(g, p) -> true],
            always_accept,
            partition,
            5
        )

        steps = 0
        for state in mc
            steps += 1
            @test state isa Partition
        end
        @test steps == 5
    end
end
