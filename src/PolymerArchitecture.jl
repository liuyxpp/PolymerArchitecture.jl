module PolymerArchitecture

using Polymer

using LightGraphs
using LightGraphs.Experimental

include("subtree.jl")
export
    Subtree,
    front_vertex,
    front_vertices,
    front_edges,
    find_leaf_neighbor,
    find_neighbors!,
    find_neighbors,
    induced_subtree,
    merge_subtrees

include("equivalent.jl")
export
    equivalent_blockend,
    equivalent_block,
    is_symmetric_subtree,
    is_isomorphic_subtree,
    is_equivalent_subtree,
    count_isomorphic_subtree,
    all_isomorphic_subtree

include("process.jl")
export
    group_isomorphic_subtrees,
    process_isomorphic_subtree_groups,
    process_isolated_isomorphic_subtrees,
    process_equivalent_subtrees,
    process_leaf

include("semi_equivalent.jl")
export
    group_isomorphic_branches,
    group_semi_equivalent_subtrees

end # module