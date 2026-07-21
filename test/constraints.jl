@testset "Population Constraint Functions" begin
    graph = BaseGraph(square_grid_filepath, "population")
    partition = Partition(graph, "assignment")

    @test within_population_bounds(partition, 30, 50)
    @test !within_population_bounds(partition, 42, 50)

    @test within_percent_of_ideal_population(graph, partition, 0.1)

    unbalanced = clone_for_update(partition)
    unbalanced.dist_populations[1] = 50
    @test !within_percent_of_ideal_population(graph, unbalanced, 0.1)

    pop_c = PopulationConstraint(0.1)
    @test satisfies_constraint(pop_c, graph, partition)
    @test !satisfies_constraint(pop_c, graph, unbalanced)
end

@testset "Contiguity Constraint Function" begin
    graph = BaseGraph(cols_grid_filepath, "population")
    partition = Partition(graph, "assignment")

    cont_c = ContiguityConstraint()
    @test satisfies_constraint(cont_c, graph, partition)
end
