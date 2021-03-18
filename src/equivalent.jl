"""
    equivalent_blockend(BlockCopolymerGraph, v1, v2; vmap1, vmap2)
    equivalent_blockend(BlockCopolymerGraph, v1, v2, vmap1, vmap2)

Return whether two blockends are equivalent. We say two blockends are equivalent if the dgrees of their corresponding vertices in the graph are identical.
"""
function equivalent_blockend(bcg::BlockCopolymerGraph, v1, v2; vmap1=collect(1:nv(bcg.graph)), vmap2=collect(1:nv(bcg.graph)))
    # convert vertices to BlockCopolymerGraph vertices
    v1n = vmap1[v1]
    v2n = vmap2[v2]
    # check the degree
    (degree(bcg.graph, v1n) == degree(bcg.graph, v2n)) && return true
    # Otherwise
    return false
end

equivalent_blockend(bcg::BlockCopolymerGraph, v1, v2, vmap1, vmap2) = equivalent_blockend(bcg, v1, v2; vmap1=vmap1, vmap2=vmap2)

"""
    equivalent_block(BlockCopolymerGraph, e1, e2; vmap1, vmap2)
    equivalent_block(BlockCopolymerGraph, e1, e2, vmap1, vmap2)

Return whether two blocks are equivalent. We say two blocks are equivalent if the species and lengths of the blocks are identical.

Edges `e1` and `e2` can be in the induced subtrees of the input `bcg.graph`.
The maps from subtrees to original graph are provided by `vmap1` and `vmap2` for `e1` and `e2`, respectively.
By default, `vmap1` and `vmap2` are identical maps, which means `e1` and `e2` are both indexed in the input `bcg.graph`.
"""
function equivalent_block(bcg::BlockCopolymerGraph, e1::Edge, e2::Edge; vmap1=collect(1:nv(bcg.graph)), vmap2=collect(1:nv(bcg.graph)))
    # convert edges to BlockCopolymerGraph edges
    v11, v12 = src(e1), dst(e1)
    e1n = (vmap1[v11], vmap1[v12]) |> Polymer._sort_tuple2
    # convert to BlockCopolymerGraph
    v21, v22 = src(e2), dst(e2)
    e2n = (vmap2[v21], vmap2[v22]) |> Polymer._sort_tuple2
    # compare if this two edges are blocks with identical properties
    block1 = bcg.edge2block[e1n]
    block2 = bcg.edge2block[e2n]
    # 1. Edges correspond to same block
    (block1 == block2) && return true
    # 2. Edges correspond to two blocks but with same specie and length
    (specie(block1) == specie(block2) && block1.f == block2.f) && return true
    # Otherwise
    return false
end

equivalent_block(bcg::BlockCopolymerGraph, e1::Edge, e2::Edge, vmap1, vmap2) = equivalent_block(bcg, e1, e2; vmap1=vmap1, vmap2=vmap2)

"""
    is_symmetric_tree(BlockCopolymerGraph, subtree, v1, v2)

Return whether the subtree is symmetric.

We say a tree is symmetric about two vertices v1 and v2, if its depth first trasversal path starting from v1 and starting from v2 are isomorphic, in the sense that all blocks are equivalent, i.e. the specie and length of the block are identical, and all blockend are equivalent, i.e. the degree of the corresponding vertex is identical.

## Arguments
* `bcg`: a graph representation of a BlockCopolymer object.
* `subtree`: an induced subtree of the `bcg.graph`.
* `v1`, `v2`: are two vertices in the `subtree` which are indexed in the orginal `bcg.graph`.
"""
function is_symmetric_subtree(bcg::BlockCopolymerGraph, subtree::Subtree, v1, v2)
    # If tree only has two vertices (one common edge), then it is symmetric.
    (nv(subtree.graph) == 2) && return true

    v1i = subtree.rvmap[v1]
    v2i = subtree.rvmap[v2]
    vmap = subtree.vmap

    dfs1 = dfs_tree(subtree.graph, v1i)
    dfs2 = dfs_tree(subtree.graph, v2i)
    return has_isomorph(dfs1, dfs2;
            vertex_relation=(v1,v2)->equivalent_blockend(bcg,v1,v2,vmap,vmap),
            edge_relation=(e1,e2)->equivalent_block(bcg,e1,e2,vmap,vmap))
end

"""
    is_isomorphic_tree(BlockCopolymerGraph, subtree1, subtree2)

Return whether a pair of subtrees, `subtree1` and `subtree2`, are isomorphic.

We say a pair of subtrees are isomorphic, named isomorphic subtrees, if their depth-first trasversal paths (directed graphs) starting from their front vertices are isomorphic, in the sense that all blocks are equivalent and all blockends are equivalent.
"""
function is_isomorphic_subtree(bcg::BlockCopolymerGraph, subtree1::Subtree, subtree2::Subtree)
    dfs1 = dfs_tree(subtree1.graph, subtree1.vi)
    dfs2 = dfs_tree(subtree2.graph, subtree2.vi)
    vmap1, vmap2 = subtree1.vmap, subtree2.vmap
    return has_isomorph(dfs1, dfs2;
            vertex_relation=(v1,v2)->equivalent_blockend(bcg,v1,v2,vmap1,vmap2),
            edge_relation=(e1,e2)->equivalent_block(bcg,e1,e2,vmap1,vmap2))
end

"""
    is_equivalent_subtree(BlockCopolymerGraph, subtree1, subtree2; check=false)

Return whether `subtree1` and `subtree2` are equivalent.

We say two subtrees are equivalent, named as equivalent subtrees, if they are isomorphic AND the vertices not in these two subtrees forms an induced subtree that is symmetric about the two front vertices of `subtree1` and `subtree2`.

Note that in this function, subtree1 and subtree2 are implicitly assumed to be isomorphic. Because the only way to obtain subtree1 and subtree2 is by connecting one or more equivalent subtrees in our applications.

In other applications, one may have to check whether `subtree1` and `subtree2` are isomorphic.
"""
function is_equivalent_subtree(bcg::BlockCopolymerGraph, subtree1::Subtree, subtree2::Subtree; check=false)
    if check
        is_isomorphic_subtree(bcg, subtree1, subtree2) || return false
    end

    subtree = induced_subtree(bcg, subtree1, subtree2)
    return is_symmetric_subtree(bcg, subtree, subtree1.v, subtree2.v)
end

"""
    count_isomorphic_subtree(BlockCopolymerGraph, subtree)

Return the number of subtrees which are isomorphic to the input `subtree`.

Note for a same set of vertices, there may be multiple isomorphic trees, thus the number will always be larger or equal to `length(all_isomorphic_subtree(bcg, subtree))`.
"""
function count_isomorphic_subtree(bcg::BlockCopolymerGraph, subtree::Subtree)
    vmap1 = collect(1:nv(bcg.graph))
    vmap2 = subtree.vmap
    return count_induced_subgraphisomorph(bcg.graph, subtree.graph;
            vertex_relation=(v1,v2)->equivalent_blockend(bcg,v1,v2,vmap1,vmap2),
            edge_relation=(e1,e2)->equivalent_block(bcg,e1,e2,vmap1,vmap2))
end

"""
    all_isomorphic_subtree(BlockCopolymerGraph, subtree)

Find all subtrees which are isomorphic to the input `subtree` in `bcg.graph`.

Note that we keep one subtree for the same set of vertices.
"""
function all_isomorphic_subtree(bcg::BlockCopolymerGraph, subtree::Subtree)
    T = eltype(bcg.graph)
    vmap1 = collect(1:nv(bcg.graph))
    vmap2 = subtree.vmap
    # ilist is a Channel. After collect, it is a list of list of 2-element tuples.
    # [[(3, 1),(8, 2),(5, 3)], [(6, 1),(4, 2),(9, 3)]
    # In each tuple, the first one corresponds to the vertex indexed in the first argument of `all_induced_subgraphisomorph`, the second one to the second argument.
    # Here, [3,8,5] is isomorphic to [6,4,9] in the first argument graph. And a mapping relation exists, such as 3 <=> 6, 8 <=> 4, and 5 <=> 9.
    ilist = all_induced_subgraphisomorph(bcg.graph, subtree.graph;
        vertex_relation=(v1,v2)->equivalent_blockend(bcg,v1,v2,vmap1,vmap2),
        edge_relation=(e1,e2)->equivalent_block(bcg,e1,e2,vmap1,vmap2))

    # We collect all vertices in the first argument, i.e. bcg.graph.
    # The front vertex is the one mapping to `subtree.v`.
    vertices_list = []
    front_vertex_list = T[]
    for tree in collect(ilist)
        vs = T[]
        for (v1, v2) in tree
            push!(vs, v1)
            (vmap2[v2] == subtree.v) && push!(front_vertex_list, v1)
        end
        push!(vertices_list, vs)
    end

    # Construct `Subtree` objects. Only one subtree is constructed for a same set of vertices.
    subtrees = [subtree]
    for i in 1:length(vertices_list)
        newtree = true
        vs = vertices_list[i]
        fv = front_vertex_list[i]
        for tree in subtrees
            (Set(vs) == Set(tree.vmap)) && (newtree=false; break)
        end
        newtree && push!(subtrees, Subtree(bcg, vs, fv))
    end

    return subtrees
end