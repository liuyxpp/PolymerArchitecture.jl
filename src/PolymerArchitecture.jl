module PolymerArchitecture

using Polymer

using Graphs
using Graphs.Experimental

include("types.jl")
export
    AbstractPolymerArchitecture,
    NonCyclicArchitecture,
    CyclicArchitecture,
    LinearArchitecture,
    BranchedArchitecture,
    StarArchitecture,
    CombArchitecture,
    GeneralBranchedArchitecture,
    RingArchitecture
export
    iscyclicchain,
    isnoncyclicchain,
    islinearchain

include("edge.jl")
include("graph.jl")
export
    BlockCopolymerGraph,
    SPECIECOLORS
export
    build_graph,
    chaintype

include("subtree.jl")
export
    Subtree,
    front_vertex,
    front_vertices,
    pre_front_vertices,
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
    is_leaf_node,
    is_leaf_edge,
    all_leafs,
    group_isomorphic_subtrees,
    process_isomorphic_subtree_groups,
    process_isolated_isomorphic_subtrees,
    process_equivalent_subtrees,
    process_leaf,
    can_merge_eq,
    merge_elements_eq,
    all_distinct_eq,
    find_all_equivalent_subtrees

include("semi_equivalent.jl")
export
    group_isomorphic_branches,
    group_semi_equivalent_subtrees,
    process_semi_equivalent_subtree_group,
    process_semi_equivalent_subtrees

import SnoopPrecompile

SnoopPrecompile.@precompile_all_calls begin
    AB = diblock_chain()
    ABG = BlockCopolymerGraph(AB)
end

end # module