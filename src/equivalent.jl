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
function is_symmetric_tree(bcg::BlockCopolymerGraph, subtree::Subtree, v1, v2)
	# If tree only has two vertices (one common edge), then it is symmetric.
	(nv(tree) == 2) && return true

	dfs1 = dfs_tree(tree, v1)
	dfs2 = dfs_tree(tree, v2)
	return has_isomorph(dfs1, dfs2;
				vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,vmap,vmap,bcg),
				edge_relation=(e1,e2)->equivalent_block(e1,e2,vmap,vmap,bcg))
end