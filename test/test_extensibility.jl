@testset "Extensibility stubs" begin
    # Minimal stub graph implementing the AbstractGraph accessor contract
    struct StubGraph <: GerryChain.AbstractGraph
        n_nodes::Int
        n_edges::Int
        pop_total::Int
        pops::Vector{Int}
        srcs::Vector{Int}
        dsts::Vector{Int}
        nbrs::Vector{Vector{Int}}
        penalties::Vector{Float64}
        regions::Dict{String,Vector{UInt32}}
    end

    GerryChain.num_nodes(g::StubGraph) = g.n_nodes
    GerryChain.num_edges(g::StubGraph) = g.n_edges
    GerryChain.total_pop(g::StubGraph) = g.pop_total
    GerryChain.populations(g::StubGraph) = g.pops
    GerryChain.edge_src(g::StubGraph) = g.srcs
    GerryChain.edge_dst(g::StubGraph) = g.dsts
    GerryChain.neighbors(g::StubGraph) = g.nbrs
    GerryChain.edge_penalties(g::StubGraph) = g.penalties
    GerryChain.has_region(g::StubGraph, col::AbstractString) = haskey(g.regions, String(col))
    GerryChain.region_ids(g::StubGraph, col::AbstractString) = g.regions[String(col)]

    mutable struct StubPartition <: GerryChain.AbstractPartition
        n_dists::Int
        n_cut::Int
        asg::Vector{Int}
        pops::Vector{Int}
        cuts::Vector{Int}
        adj::SparseMatrixCSC{Int,Int}
        nodes::Vector{BitSet}
        par::Union{GerryChain.AbstractPartition,Nothing}
    end

    GerryChain.num_dists(p::StubPartition) = p.n_dists
    GerryChain.num_cut_edges(p::StubPartition) = p.n_cut
    GerryChain.assignments(p::StubPartition) = p.asg
    GerryChain.dist_populations(p::StubPartition) = p.pops
    GerryChain.cut_edges(p::StubPartition) = p.cuts
    GerryChain.dist_adj(p::StubPartition) = p.adj
    GerryChain.dist_nodes(p::StubPartition) = p.nodes
    Base.parent(p::StubPartition) = p.par

    function GerryChain.clone_for_update(p::StubPartition)
        return StubPartition(
            p.n_dists,
            p.n_cut,
            copy(p.asg),
            copy(p.pops),
            copy(p.cuts),
            copy(p.adj),
            BitSet[copy(s) for s in p.nodes],
            p,
        )
    end

    stub_g = StubGraph(
        3,
        2,
        30,
        [10, 10, 10],
        [1, 2],
        [2, 3],
        [[2], [1, 3], [2]],
        zeros(2),
        Dict("county" => UInt32[1, 1, 2]),
    )
    stub_p = StubPartition(
        2,
        1,
        [1, 1, 2],
        [20, 10],
        [0, 1],
        sparse([0 1; 1 0]),
        [BitSet([1, 2]), BitSet([3])],
        nothing,
    )

    @test stub_g isa GerryChain.AbstractGraph
    @test stub_p isa GerryChain.AbstractPartition
    @test Partition <: GerryChain.AbstractPartition
    @test BaseGraph <: GerryChain.AbstractGraph

    @test num_nodes(stub_g) == 3
    @test num_edges(stub_g) == 2
    @test total_pop(stub_g) == 30
    @test populations(stub_g) == [10, 10, 10]
    @test has_region(stub_g, "county")
    @test region_ids(stub_g, "county")[3] == UInt32(2)
    @test num_dists(stub_p) == 2
    @test assignments(stub_p) == [1, 1, 2]

    # PopulationConstraint accepts abstracts
    pop_c = PopulationConstraint(stub_g, stub_p, 0.1)
    @test pop_c isa PopulationConstraint
    @test satisfy_constraint(pop_c, 15, 15)

    # clone_for_update isolation on stub
    cloned = clone_for_update(stub_p)
    @test parent(cloned) === stub_p
    cloned.asg[1] = 2
    @test stub_p.asg[1] == 1
    @test cloned.nodes[1] !== stub_p.nodes[1]
end

@testset "clone_for_update on Partition" begin
    graph = BaseGraph(square_grid_filepath, "population")
    partition = Partition(graph, "assignment")
    original_asg = copy(partition.assignments)
    original_dist1 = copy(partition.dist_nodes[1])

    cloned = clone_for_update(partition)
    @test parent(cloned) === partition
    @test cloned.assignments !== partition.assignments
    @test cloned.dist_nodes[1] !== partition.dist_nodes[1]
    @test cloned.dist_nodes[1] == partition.dist_nodes[1]

    cloned.assignments[1] = 99
    @test partition.assignments[1] == original_asg[1]
    push!(cloned.dist_nodes[1], 99)
    @test !(99 in partition.dist_nodes[1])
    @test partition.dist_nodes[1] == original_dist1

    # parent snapshot path avoids deepcopy
    proposal = FlipProposal(
        1,
        partition.assignments[1],
        2,
        partition.dist_populations[partition.assignments[1]] - graph.populations[1],
        partition.dist_populations[2] + graph.populations[1],
        setdiff(partition.dist_nodes[partition.assignments[1]], 1),
        union(partition.dist_nodes[2], 1),
    )
    before = copy(partition.assignments)
    update_partition!(partition, graph, proposal, true)
    @test partition.parent isa Partition
    @test partition.parent.assignments == before
    @test partition.parent.parent === nothing
end

@testset "LightGraphs coexistence (no dual-export clash)" begin
    # These must resolve without Ambiguous/UndefVarError after using both packages.
    @test LightGraphs.AbstractGraph isa Type
    @test GerryChain.AbstractGraph isa Type
    @test LightGraphs.AbstractGraph !== GerryChain.AbstractGraph

    g = SimpleGraph(2)
    add_edge!(g, 1, 2)
    @test neighbors(g, 1) == [2]  # LightGraphs.neighbors

    bg = BaseGraph(square_grid_filepath, "population")
    @test GerryChain.neighbors(bg) === bg.neighbors
    @test length(GerryChain.kruskal_mst(bg, [1], [bg.edge_src[1], bg.edge_dst[1]], [0.0])) <= 1
end
