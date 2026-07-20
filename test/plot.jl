@testset "Plotting tests" begin
    # 2D matrix of per-district scores over 10 chain steps (10 x 4)
    district_scores = rand(10, 4)
    # 1D vector of plan-wide score over 10 steps
    plan_scores = rand(10)

    @testset "score_boxplot()" begin
        ax = score_boxplot(district_scores)
        @test ax !== nothing
    end

    @testset "score_histogram()" begin
        ax = score_histogram(plan_scores; comparison_scores=[("Plan A", 0.5)])
        @test ax !== nothing
    end
end
