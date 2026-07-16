# Extending GerryChain

GerryChain is designed as an extensible redistricting kernel. Downstream packages
can implement custom graphs and partitions without forking internals.

## Abstract types

| Type | Default concrete type |
|------|-----------------------|
| `AbstractGraph` | `BaseGraph` |
| `AbstractPartition` | `Partition` |
| `AbstractProposal` | `RecomProposal`, `FlipProposal`, `DummyProposal` |
| `AbstractConstraint` | `PopulationConstraint`, `ContiguityConstraint` |
| `AbstractScore` | `DistrictAggregate`, `DistrictScore`, `PlanScore`, `CompositeScore` |

## Accessor contract

Implement these methods on custom subtypes:

**Graph:** `num_nodes`, `num_edges`, `total_pop`, `populations`, `edge_src`,
`edge_dst`, `neighbors`, `edge_penalties`, `has_region`, `region_ids`.

**Partition:** `num_dists`, `num_cut_edges`, `assignments`, `dist_populations`,
`cut_edges`, `dist_adj`, `dist_nodes`, `Base.parent` (or the `parent` field),
and `clone_for_update`.

## Copy-on-write updates

Prefer returning a new partition from proposal helpers:

```julia
new_p = clone_for_update(partition)   # copies arrays/BitSets; parent = partition
update_partition!(new_p, graph, proposal, false)
```

`clone_for_update` does **not** copy the graph and does not recurse into an
existing parent chain. The in-place `update_partition!(..., copy_parent=true)`
path stores a field-wise parent snapshot (no recursive `deepcopy`).

## Region-aware / weighted ReCom

Register dense region columns at load time or afterward:

```julia
graph = BaseGraph(path, "population"; region_columns = ["COUNTYID", "MUNIID"])
set_edge_penalty!(graph, u, v, 10.0)
# or: set_edge_penalties_from_pairs!(graph, Dict((u, v) => w, ...))

recom_chain(
    graph, partition, pop_constraint, num_steps, scores;
    region_surcharges = Dict("COUNTYID" => 1.0, "MUNIID" => 0.5),
)
```

MST edge weights are:

`rand(rng) + edge_penalties[e] + Σ surcharge[col]` when endpoints differ in region `col`.

## LightGraphs coexistence

GerryChain depends on LightGraphs. Several names exist in both packages
(`AbstractGraph`, `neighbors`, `kruskal_mst`). Those symbols are **not**
re-exported from GerryChain so `using GerryChain` + `using LightGraphs` stays
unambiguous:

- Use `GerryChain.AbstractGraph` for the dual-graph abstract type
- Use `graph.neighbors` (field) or `GerryChain.neighbors(graph)` for adjacency lists
- Use `GerryChain.kruskal_mst` for the dual-graph MST helper

`BaseGraph`, `weighted_kruskal_mst`, and `random_kruskal_mst` remain exported.

## Downstream compatibility (0.1.3-style call sites)

- The pre-0.2.0 **10-arg** `BaseGraph(...)` constructor still works; it fills
  `edge_penalties` with zeros and `region_cols` with an empty dict.
- `weighted_kruskal_mst(graph, edges, nodes, weights::AbstractVector)` remains
  available for callers that pass precomputed weights (distinct from the
  RNG-based region-aware method).
