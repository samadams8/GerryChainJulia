@testset "Population Constraint Functions" begin
    graph = BaseGraph(square_grid_filepath, "population")
    partition = Partition(graph, "assignment")

    @test within_population_bounds(partition, 30, 50)
    @test !within_population_bounds(partition, 42, 50)

    @test within_percent_of_ideal_population(graph, partition, 0.1)

    unbalanced = clone_for_update(partition)
    unbalanced.dist_populations[1] = 50
    @test !within_percent_of_ideal_population(graph, unbalanced, 0.1)

    validator = population_constraint(0.1)
    @test validator(graph, partition)
end

@testset "Contiguity Constraint Function" begin
    graph = BaseGraph(cols_grid_filepath, "population")
    partition = Partition(graph, "assignment")

    discont_proposal =
        FlipProposal(5, 1, 0, 30, 42, BitSet([1, 9, 13]), BitSet([0, 4, 5, 8, 12]))

    @test !is_contiguous_flip(graph, partition, discont_proposal)

    cont_proposal =
        FlipProposal(1, 1, 0, 30, 42, BitSet([5, 9, 13]), BitSet([0, 1, 4, 8, 12]))

    @test is_contiguous_flip(graph, partition, cont_proposal)
end
