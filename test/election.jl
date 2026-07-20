@testset "Election tests" begin
    graph = BaseGraph(square_grid_filepath, "population")
    partition = Partition(graph, "assignment")

    @testset "vote_count and vote_share" begin
        d_counts = vote_count(graph, partition, "electionD")
        r_counts = vote_count(graph, partition, "electionR")

        @test d_counts[1] == 6.0
        @test r_counts[1] == 6.0

        d_shares = vote_share(graph, partition, "electionD", "population")
        @test d_shares[1] ≈ 6.0 / 41.0
    end

    @testset "seats_won" begin
        d_votes = [10.0, 15.0, 5.0, 8.0]
        r_votes = [8.0, 5.0, 15.0, 10.0]

        @test seats_won(d_votes, r_votes) == 2
        @test seats_won(r_votes, d_votes) == 2
    end

    @testset "wasted_votes and efficiency_gap" begin
        w1, w2 = wasted_votes(60.0, 40.0)
        @test w1 == 10.0
        @test w2 == 40.0

        d_votes = [60.0, 60.0]
        r_votes = [40.0, 40.0]
        eg = efficiency_gap(d_votes, r_votes)
        @test eg == -0.3
    end

    @testset "mean_median" begin
        shares = [0.4, 0.5, 0.6]
        @test mean_median(shares) ≈ 0.0
    end
end
