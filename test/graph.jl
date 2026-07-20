using LightGraphs
using Random
using DataStructures

"""
    the test graph being loaded is labeled as

    01 - 02 - 03 - 04
     |    |   |    |
    05 - 06 - 07 - 08
     |    |   |    |
    09 - 10 - 11 - 12
     |    |   |    |
    13 - 14 - 15 - 16

    the initial assignment is

    01 - 01 - 02 - 02
     |    |   |    |
    01 - 01 - 02 - 02
     |    |   |    |
    03 - 03 - 04 - 04
     |    |   |    |
    03 - 03 - 04 - 04

    the population distribution is

    20 - 10 - 10 - 10
     |    |   |    |
    01 - 10 - 20 - 01
     |    |   |    |
    20 - 10 - 10 - 10
     |    |   |    |
    01 - 10 - 01 - 20
"""

@testset "Graph tests" begin
    @testset "Bad extension" begin
        # file should have a .json or .shp extension
        @test_throws DomainError BaseGraph("nonexistent.txt", "population")
    end

    @testset "Reading node attributes from shapefile" begin
        # file should have a .json or .shp extension
        table = GerryChain.read_table(square_shp_filepath)
        node_attributes = GerryChain.all_node_properties(table)
        correct_attributes = [ # refer to maps/make_simple_shp.py
            Dict("assignment" => 1, "population" => 2),
            Dict("assignment" => 2, "population" => 4),
            Dict("assignment" => 3, "population" => 6),
            Dict("assignment" => 4, "population" => 8),
        ]
        @test node_attributes == correct_attributes
    end

    @testset "BaseGraph from shp() - rook adjacency" begin
        graph = BaseGraph(square_shp_filepath, "population")
        @test graph.num_nodes == 4
        @test graph.num_edges == 4
        @test graph.total_pop == 20
        @test graph.populations == [2, 4, 6, 8]
        @test graph.adj_matrix[1, 2] != 0
        @test graph.adj_matrix[1, 3] != 0
        # upper left corner and bottom right corner should be non-adjacent
        @test graph.adj_matrix[1, 4] == 0
    end

    @testset "BaseGraph from shp() - queen adjacency" begin
        graph = BaseGraph(square_shp_filepath, "population", adjacency = "queen")
        @test graph.num_nodes == 4
        @test graph.num_edges == 6 # queen adjacency means all 6 edges
        @test graph.total_pop == 20

        @test graph.populations == [2, 4, 6, 8]
        # with queen adjacency, all squares should be adjacent to each other
        @test graph.adj_matrix[1, 2] != 0
        @test graph.adj_matrix[1, 3] != 0
        @test graph.adj_matrix[1, 4] != 0
        @test graph.adj_matrix[2, 3] != 0
        @test graph.adj_matrix[2, 4] != 0
        @test graph.adj_matrix[3, 4] != 0
    end

    graph = BaseGraph(square_grid_filepath, "population")

    @test graph.num_nodes == 16
    @test graph.num_edges == 24
    @test graph.total_pop == 164

    @testset "Populations" begin
        for i in [1, 7, 9, 16]
            @test graph.populations[i] == 20
        end
        for i in [5, 8, 13, 15]
            @test graph.populations[i] == 1
        end
        for i in [2, 3, 4, 6, 10, 11, 12, 14]
            @test graph.populations[i] == 10
        end
    end

    # test adjacencies
    @test graph.adj_matrix[1, 2] != 0
    @test graph.adj_matrix[1, 5] != 0
    @test graph.adj_matrix[1, 6] == 0

    @testset "Graph Adjacency Symmetry" begin
        for i = 1:graph.num_nodes
            for j = 1:graph.num_nodes
                @test graph.adj_matrix[i, j] == graph.adj_matrix[j, i]
            end
        end
    end

    @testset "Check type of district assignments - get_assignments()" begin
        partition = Partition(graph, "assignment")
        # assignment that is a Float should throw an error
        foreach(d -> d["assignment"] *= 1.0, graph.attributes) # convert Int to Float
        @test_throws DomainError GerryChain.get_assignments(graph.attributes, "assignment")
    end

    # test the edge arrays
    @test graph.edge_src[graph.adj_matrix[10, 11]] in (10, 11)
    @test graph.edge_dst[graph.adj_matrix[11, 10]] in (10, 11)
    @test graph.edge_src[graph.adj_matrix[9, 13]] in (9, 13)
    @test graph.edge_dst[graph.adj_matrix[13, 9]] in (13, 9)
    @test length(graph.edge_src) == graph.num_edges
    @test length(graph.edge_dst) == graph.num_edges

    # test the node neighbors
    @test sort(graph.neighbors[1]) == [2, 5]
    @test sort(graph.neighbors[6]) == [2, 5, 7, 10]
    @test sort(graph.neighbors[14]) == [10, 13, 15]

    # test the simple graph
    @test LightGraphs.nv(graph.simple_graph) == graph.num_nodes
    @test LightGraphs.ne(graph.simple_graph) == graph.num_edges

    # test induced_subgraph
    @testset "induced_subgraph_edges()" begin
        induced_edges = induced_subgraph_edges(graph, [1, 2, 3, 4])
        @test sort(induced_edges) == sort([1, 3, 5])

        induced_vertices = Set{Int}()
        for edge in induced_edges
            push!(induced_vertices, graph.edge_src[edge], graph.edge_dst[edge])
        end
        @test induced_vertices == Set{Int}([1, 2, 3, 4])

        @test induced_subgraph_edges(graph, Int[]) == Int[]
        @test isempty(induced_subgraph_edges(graph, [5]))  # isolated in induced sense if no self-loop
    end
    @test_throws ArgumentError induced_subgraph_edges(graph, [1, 1, 4])

    # get_subgraph_population()
    @test get_subgraph_population(graph, BitSet([1, 2, 3, 4])) == 50
    @test get_subgraph_population(graph, BitSet([5])) == 1
    @test get_subgraph_population(graph, BitSet([1, 6, 11, 16])) == 60

    # test that attributes can be accessed
    @test graph.attributes[1]["purple"] == 15
    @test graph.attributes[1]["pink"] == 5

    @testset "edge_penalties and region columns" begin
        @test length(edge_penalties(graph)) == graph.num_edges
        @test all(iszero, edge_penalties(graph))

        u, v = 1, 2
        set_edge_penalty!(graph, u, v, 5.0)
        eid = graph.adj_matrix[u, v]
        @test edge_penalties(graph)[eid] == 5.0

        set_edge_penalties_from_pairs!(graph, Dict((2, 3) => 7.5, (3, 4) => 1.0))
        @test edge_penalties(graph)[graph.adj_matrix[2, 3]] == 7.5
        @test edge_penalties(graph)[graph.adj_matrix[3, 4]] == 1.0

        @test !has_region(graph, "county")
        add_region_column!(graph, "county", [1, 1, 2, 2, 1, 1, 2, 2, 3, 3, 4, 4, 3, 3, 4, 4])
        @test has_region(graph, "county")
        @test length(region_ids(graph, "county")) == graph.num_nodes
        @test region_ids(graph, "county")[1] == region_ids(graph, "county")[2]
        @test region_ids(graph, "county")[1] != region_ids(graph, "county")[3]

        # cross-boundary detection
        srcs, dsts = edge_src(graph), edge_dst(graph)
        cross = 0
        for e = 1:num_edges(graph)
            if region_ids(graph, "county")[srcs[e]] != region_ids(graph, "county")[dsts[e]]
                cross += 1
            end
        end
        @test cross > 0

        # null / empty labels encode as sentinel 0
        add_region_column!(graph, "muni", ["A", missing, "", "A", nothing, "B", "B", "A",
                                           "A", "B", "A", "B", "A", "B", "A", "B"])
        muni_ids = region_ids(graph, "muni")
        @test muni_ids[2] == UInt32(0)
        @test muni_ids[3] == UInt32(0)
        @test muni_ids[5] == UInt32(0)

        # configure_mst_weights! tests
        configure_mst_weights!(graph; region_surcharges=Dict("county" => 10.0))
        # Ensure base weights cache exists
        @test graph._mst_base_weights[] !== nothing
        # Test that cross-boundary edges in "county" get surcharge, in-boundary do not
        for e = 1:num_edges(graph)
            u, v = srcs[e], dsts[e]
            expected = edge_penalties(graph)[e]
            if region_ids(graph, "county")[u] != region_ids(graph, "county")[v]
                expected += 10.0
            end
            @test graph._mst_base_weights[][e] == expected
        end

        # Test invalidation of cache on set_edge_penalty!
        set_edge_penalty!(graph, 1, 2, 99.0)
        @test graph._mst_base_weights[] === nothing

        @test muni_ids[1] == muni_ids[4]
        @test muni_ids[1] != UInt32(0)
        @test muni_ids[1] != muni_ids[6]

        g2 = BaseGraph(
            square_grid_filepath,
            "population";
            region_columns = ["assignment"],
        )
        @test has_region(g2, "assignment")
        @test length(region_ids(g2, "assignment")) == g2.num_nodes
    end

    @testset "10-arg BaseGraph compatibility constructor" begin
        g10 = BaseGraph(
            graph.num_nodes,
            graph.num_edges,
            graph.total_pop,
            graph.populations,
            graph.adj_matrix,
            graph.edge_src,
            graph.edge_dst,
            graph.neighbors,
            graph.simple_graph,
            graph.attributes,
        )
        @test length(edge_penalties(g10)) == g10.num_edges
        @test all(iszero, edge_penalties(g10))
        @test isempty(g10.region_cols)

        # UTGC-style attributes: Vector{Dict{String,String}}
        attrs_str = [Dict("county" => "A"), Dict("county" => "B")]
        simple = SimpleGraph(2)
        add_edge!(simple, 1, 2)
        src, dst = GerryChain.edges_from_graph(simple)
        adj = GerryChain.adjacency_matrix_from_graph(simple)
        nbrs = GerryChain.neighbors_from_graph(simple)
        g_str = BaseGraph(
            2, 1, 2, [1, 1], adj, src, dst, nbrs, simple, attrs_str
        )
        @test g_str.attributes isa Array{Dict{String,Any}}
        @test g_str.attributes[1]["county"] == "A"
        @test length(edge_penalties(g_str)) == 1
    end

    @testset "lazy attribute cache" begin
        g = BaseGraph(square_grid_filepath, "population")
        col1 = GerryChain._attribute_vector(g, "purple")
        col2 = GerryChain._attribute_vector(g, "purple")
        @test col1 === col2
        public = attribute_vector(g, "purple")
        @test public == col1
        @test public !== col1

        partition = Partition(g, "assignment")
        score = DistrictAggregate("purple")
        expected = sum(g.attributes[n]["purple"] for n in partition.dist_nodes[1])
        @test eval_score_on_district(g, partition, score, 1) == expected

        old = g.attributes[1]["purple"]
        set_attribute!(g, 1, "purple", old + 1)
        @test !haskey(g._attr_cache, "purple")
        col3 = GerryChain._attribute_vector(g, "purple")
        @test col3[1] == old + 1
        @test col3 !== col1
    end
end
