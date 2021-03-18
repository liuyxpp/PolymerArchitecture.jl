module PolymerArchitecture

using Polymer

using LightGraphs
using LightGraphs.Experimental

include("subtree.jl")
export
    Subtree,
    front_vertex,
    front_edges,
    find_neighbors!,
    find_neighbors,
    induced_subtree

include("equivalent.jl")
export
    equivalent_blockend,
    equivalent_block,
    is_symmetric_tree

end # module