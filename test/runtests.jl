using GerryChain
using LibGEOS: LibGEOS
using LibSpatialIndex: LibSpatialIndex
using Test
using LightGraphs
using JSON
using HDF5
using Logging
using SparseArrays
using Random
using DataStructures
using ProgressBars

const testdir = dirname(@__FILE__)
square_grid_filepath = joinpath(testdir, "maps", "test_grid_4x4.json")
cols_grid_filepath = joinpath(testdir, "maps", "cols_grid_4x4.json")
square_shp_filepath = joinpath(testdir, "maps", "simple_squares.shp")

tests = [
    "graph",
    "partition",
    "balance_edges",
    "geo",
    "constraints",
    "flip",
    "recom",
    "scores",
    "accept",
    "election",
    "plot",
    "test_extensibility",
]

@testset "GerryChainJulia" begin
    for t in tests
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end
