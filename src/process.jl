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

Process isomorphic subtree groups which are obtained from `group_isomorphic_subtrees`.
"""
function process_isomorphic_subtree_groups(bcg::BlockCopolymerGraph, subtree_groups)
    equivalent_subtrees = []
    semi_equivalent_subtrees = []
    isolated_isomorphic_subtrees = []

    ntrees = [length(group) for group in subtree_groups]

    sum(ntrees) < 2 && error("At least two isomorphic subtrees in the groups.")

    for n in unique(ntrees)
        # Put all groups with same number of isomorphic subtrees into a Qgroup
        Qgroups = []
        for i in 1:length(subtree_groups)
            if ntrees[i] == n
                push!(Qgroups, subtree_groups[i])
            end
        end
        # Process Qgroups below
        if length(Qgroups) == 1
            if n == 1
                # one group with one tree, it can be only semi equivalent
                append!(semi_equivalent_subtrees, Qgroups)
            else
                # one group with more than one subtrees, these subtrees are equivalent to each other because they are connected by one common vertex.
                append!(equivalent_subtrees, Qgroups)
            end
        else
            # Combine all subtrees in a Qgroup into one new subtree.
            # We therefore can form a new group from Qgroups.
            Qgroup_new = []
            for Qgroup in Qgroups
                new_subtree = merge_subtrees(bcg, Qgroup)
                push!(Qgroup_new, new_subtree)
            end
            # These new subtrees are also isomophic, but now they do not share common front vertex.
            push!(isolated_isomorphic_subtrees, Qgroup_new)
        end
    end

    return equivalent_subtrees, semi_equivalent_subtrees, isolated_isomorphic_subtrees
end