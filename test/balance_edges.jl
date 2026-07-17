
""" Refer to test/graph.jl to see the test graph being loaded
"""

# test random_weighted_kruskal_mst
@testset "Random Kruskal MST" begin
    graph = BaseGraph(square_grid_filepath, "population")

    rng = MersenneTwister(1234)
    nodes = [1, 2, 3, 4, 5, 6, 7, 8]
    edges = [
        graph.adj_matrix[1, 2],
        graph.adj_matrix[2, 3],
        graph.adj_matrix[3, 4],
        graph.adj_matrix[5, 6],
        graph.adj_matrix[6, 7],
        graph.adj_matrix[7, 8],
        graph.adj_matrix[1, 5],
        graph.adj_matrix[2, 6],
        graph.adj_matrix[3, 7],
        graph.adj_matrix[4, 8],
    ]

    mst = random_kruskal_mst(graph, edges, nodes, rng)
    @test length(mst) == length(nodes) - 1
    @test begin # are there loops in the tree?
        # find by union-find algorithm
        connected_vs = DisjointSets{Int}(nodes)
        cycle_found = false
        for edge in mst
            if in_same_set(connected_vs, graph.edge_src[edge], graph.edge_dst[edge])
                cycle_found = true
                break
            else
                union!(connected_vs, graph.edge_src[edge], graph.edge_dst[edge])
            end
        end
        !cycle_found
    end



end

@testset "Correctness of MST" begin
    graph = BaseGraph(square_grid_filepath, "population")

    nodes = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]
    heavy_edges = [
        Edge(1, 2),
        Edge(2, 3),
        Edge(3, 4),
        Edge(4, 8),
        Edge(7, 8),
        Edge(6, 7),
        Edge(5, 6),
        Edge(5, 9),
        Edge(9, 10),
        Edge(10, 11),
        Edge(11, 12),
        Edge(12, 16),
        Edge(15, 16),
        Edge(14, 15),
        Edge(13, 14),
    ]

    is = Array{Int}([]) # store edge indices
    weights = Array{Float64}([])
    correct_mst = BitSet()

    for (i, e) in enumerate(LightGraphs.edges(graph.simple_graph))
        push!(is, i)
        if e in heavy_edges
            push!(weights, 0.0)
            push!(correct_mst, i)
        else
            push!(weights, 1.0)
        end
    end

    @test correct_mst == GerryChain.kruskal_mst(graph, is, nodes, weights)
end

@testset "Weighted / region-aware Kruskal MST" begin
    # Triangle: nodes 1-2-3 with edges e12, e23, e13
    # Huge penalty on e13 should exclude it from the MST.
    simple = SimpleGraph(3)
    add_edge!(simple, 1, 2)
    add_edge!(simple, 2, 3)
    add_edge!(simple, 1, 3)
    edge_src_arr, edge_dst_arr = GerryChain.edges_from_graph(simple)
    adj = GerryChain.adjacency_matrix_from_graph(simple)
    nbrs = GerryChain.neighbors_from_graph(simple)
    attrs = [Dict{String,Any}("pop" => 1) for _ = 1:3]
    tri = BaseGraph(
        3,
        3,
        3,
        [1, 1, 1],
        adj,
        edge_src_arr,
        edge_dst_arr,
        nbrs,
        simple,
        attrs,
        zeros(Float64, 3),
        Dict{String,Vector{UInt32}}(),
        Dict{String,Vector{Float64}}(),
    )
    e12 = adj[1, 2]
    e23 = adj[2, 3]
    e13 = adj[1, 3]
    edges = [e12, e23, e13]
    nodes = [1, 2, 3]

    set_edge_penalty!(tri, 1, 3, 1e9)
    rng = MersenneTwister(1)
    mst = weighted_kruskal_mst(tri, edges, nodes, rng)
    @test length(mst) == 2
    @test !(e13 in mst)
    @test e12 in mst
    @test e23 in mst

    # Reset penalties; huge region surcharge should prefer in-region edges
    fill!(tri.edge_penalties, 0.0)
    add_region_column!(tri, "county", UInt32[1, 1, 2])
    # edge 1-2 is in-region; 2-3 and 1-3 cross
    rng = MersenneTwister(2)
    mst2 = weighted_kruskal_mst(
        tri,
        edges,
        nodes,
        rng;
        region_surcharges = Dict("county" => 1e9),
    )
    @test length(mst2) == 2
    @test e12 in mst2  # in-region edge must be included

    # Null-region sentinel (0): no surcharge on edges with a null endpoint
    add_region_column!(tri, "muni", ["A", missing, "B"])  # nodes 1=A, 2=null, 3=B
    fill!(tri.edge_penalties, 0.0)
    weights_null = zeros(3)
    build_mst_weights!(
        weights_null,
        tri,
        edges,
        MersenneTwister(3);
        region_surcharges = Dict("muni" => 100.0),
    )
    weights_base = zeros(3)
    build_mst_weights!(
        weights_base,
        tri,
        edges,
        MersenneTwister(3);
        region_surcharges = Dict{String,Float64}(),
    )
    # edges touch null node 2 (e12, e23) → no surcharge; only e13 (A↔B) boosted
    @test weights_null[1] ≈ weights_base[1]  # e12
    @test weights_null[2] ≈ weights_base[2]  # e23
    @test weights_null[3] ≈ weights_base[3] + 100.0  # e13

    # Legacy weights-vector overload (UTGC MST_FUNC)
    weights = [0.1, 0.2, 1e9]
    mst3 = weighted_kruskal_mst(tri, edges, nodes, weights)
    @test mst3 == GerryChain.kruskal_mst(tri, edges, nodes, Float64.(weights))
    @test !(e13 in mst3)
end

@testset "Wilson UST" begin
    graph = BaseGraph(square_grid_filepath, "population")
    nodes = [1, 2, 3, 4, 5, 6, 7, 8]
    edges = [
        graph.adj_matrix[1, 2],
        graph.adj_matrix[2, 3],
        graph.adj_matrix[3, 4],
        graph.adj_matrix[5, 6],
        graph.adj_matrix[6, 7],
        graph.adj_matrix[7, 8],
        graph.adj_matrix[1, 5],
        graph.adj_matrix[2, 6],
        graph.adj_matrix[3, 7],
        graph.adj_matrix[4, 8],
    ]
    rng = MersenneTwister(99)
    tree = wilson_ust(graph, edges, nodes, rng)
    @test length(tree) == length(nodes) - 1
    connected_vs = DisjointSets{Int}(nodes)
    cycle_found = false
    for edge in tree
        if in_same_set(connected_vs, graph.edge_src[edge], graph.edge_dst[edge])
            cycle_found = true
            break
        else
            union!(connected_vs, graph.edge_src[edge], graph.edge_dst[edge])
        end
    end
    @test !cycle_found
    # connected: one component
    roots = Set(find_root!(connected_vs, n) for n in nodes)
    @test length(roots) == 1
end

@testset "MSTScratch reuse" begin
    graph = BaseGraph(square_grid_filepath, "population")
    nodes = collect(1:16)
    edges = collect(1:graph.num_edges)
    weights = rand(MersenneTwister(1), length(edges))
    scratch = MSTScratch(length(edges), maximum(nodes))
    mst1 = kruskal_mst!(scratch, graph, edges, nodes, weights)
    mst2 = GerryChain.kruskal_mst(graph, edges, nodes, weights)
    @test mst1 == mst2
    # reuse scratch
    mst3 = kruskal_mst!(scratch, graph, edges, nodes, weights)
    @test mst3 == mst2
end
