@testset "Optimizers tests" begin
    # Dummy score function that returns cut edges
    score_fn = (g, p) -> num_cut_edges(p)

    grid = SimpleGraph(4)
    add_edge!(grid, 1, 2)
    add_edge!(grid, 2, 3)
    add_edge!(grid, 3, 4)

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
        4, 3, 40, pops, adj_mat, edge_src, edge_dst, nbrs, grid, attrs
    )

    p1 = Partition(
        2, 1, [1, 1, 2, 2], [20, 20], [0, 1, 0], spzeros(Int, 2, 2), [BitSet([1, 2]), BitSet([3, 4])], nothing
    )
    p2 = Partition(
        2, 2, [1, 2, 1, 2], [20, 20], [1, 1, 1], spzeros(Int, 2, 2), [BitSet([1, 3]), BitSet([2, 4])], nothing
    )

    @testset "greedy_accept" begin
        greedy_max = greedy_accept(score_fn; maximize=true)
        # p2 cut_edges = 2, p1 cut_edges = 1 -> p2 improves max score
        @test greedy_max(graph, p1, p2) == 1.0
        @test greedy_max(graph, p2, p1) == 0.0

        greedy_min = greedy_accept(score_fn; maximize=false)
        @test greedy_min(graph, p1, p2) == 0.0
        @test greedy_min(graph, p2, p1) == 1.0
    end

    @testset "simulated_annealing_accept" begin
        sa = simulated_annealing_accept(score_fn, 1.0; maximize=true)
        @test sa(graph, p1, p2) == 1.0
        @test sa(graph, p2, p1) < 1.0
        @test sa(graph, p2, p1) > 0.0
    end
end
