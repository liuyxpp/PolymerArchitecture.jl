"""
    group_isomorphic_subtrees(BlockCopolymerGraph, subtree)

Group isomorphic subtrees according to the front vertex. Those subtrees share same front vertex are grouped into one group.
"""
function group_isomorphic_subtrees(subtrees::AbstractVector{T}) where T<:AbstractSubtree
    unique_vertices = [subtree.v for subtree in subtrees] |> unique
    subtree_groups = []
    for v in unique_vertices
        group = T[]
        for subtree in subtrees
            (subtree.v == v) && push!(group, subtree)
        end
        push!(subtree_groups, group)
    end

    return subtree_groups
end

"""
    process_isomorphic_subtree_groups(groups)


"""
# function process_isomorphic_subtree_groups(groups, front_vertices, ntrees)
# 	equivalent_subtrees = []
# 	equivalent_front_vertices = []
# 	semi_equivalent_subtrees = []
# 	semi_equivalent_front_vertices = []
# 	isolated_isomorphic_subtrees = []
# 	isolated_isomorphic_vertices = []
	
# 	sum(ntrees) < 2 && error("At least two isomorphic subtrees in the groups.")
	
# 	unique_ntrees = ntrees |> unique
# 	for n in unique_ntrees
# 		Qgroup = []
# 		Qv = []
# 		for i in 1:length(groups)
# 			if ntrees[i] == n
# 				push!(Qgroup, groups[i])
# 				push!(Qv, front_vertices[i])
# 			end
# 		end
# 		if length(Qgroup) == 1
# 			if n == 1
# 				push!(semi_equivalent_subtrees, Qgroup)
# 				push!(semi_equivalent_front_vertices, Qv)
# 			else
# 				push!(equivalent_subtrees, Qgroup)
# 				push!(equivalent_front_vertices, Qv)
# 			end
# 		else
# 			Qgroup_new = []
# 			for j in 1:length(Qgroup)
# 				new_subtree = []
# 				for subtree in Qgroup[j]
# 					append!(new_subtree, subtree)
# 				end
# 				push!(Qgroup_new, unique(new_subtree))
# 			end
# 			push!(isolated_isomorphic_subtrees, Qgroup_new)
# 			push!(isolated_isomorphic_vertices, Qv)
# 		end
# 	end
	
# 	return equivalent_subtrees, equivalent_front_vertices, semi_equivalent_subtrees, semi_equivalent_front_vertices, isolated_isomorphic_subtrees, isolated_isomorphic_vertices
# end