"""
    group_isomorphic_subtrees(BlockCopolymerGraph, subtree)

Group isomorphic subtrees according to the front vertex. Those subtrees share same front vertex are grouped into one group.
"""
function group_isomorphic_subtrees(subtrees)
    T = typeof(first(subtrees))
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
    process_isomorphic_subtree_groups(BlockCopolymerGraph, groups)

Process isomorphic subtree groups which are obtained from `group_isomorphic_subtrees`.
"""
function process_isomorphic_subtree_groups(bcg::BlockCopolymerGraph, subtree_groups)
    equivalent_subtrees = []
    semi_equivalent_subtrees = []
    isolated_isomorphic_subtrees = []

    ntrees = [length(group) for group in subtree_groups]

    # sum(ntrees) < 2 && error("At least two isomorphic subtrees in the groups.")

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

"""
    process_isolated_isomorphic_subtrees(BlockCopolymerGraph, subtrees)

Process isolated isomorphic subtrees. Isolated isomorphic subtrees do not share any common front vertices. Thus they are isolated, i.e. at least one edge apart.
"""
function process_isolated_isomorphic_subtrees(bcg::BlockCopolymerGraph, subtrees)
    equivalent_subtrees = []
    semi_equivalent_subtrees = []

    m = length(subtrees)

    if m == 1
        # Here subtrees = [subtree]
        push!(semi_equivalent_subtress, subtrees)
    end

    if m == 2
        subtree1, subtree2 = subtrees
        issym = is_equivalent_subtree(bcg, subtree1, subtree2)
        if issym
            push!(equivalent_subtrees, [subtree1, subtree2])
        else
            push!(semi_equivalent_subtrees, [subtree1, subtree2])
        end
    end

    processed = [false for _ in 1:m]
    if m > 2
        for i in 1:m
            equiv_subtrees_subset = []
            semiequiv_subtrees_subset = []
            processed[i] && continue
            subtree1 = subtrees[i]
            for j in (i+1):m
                processed[j] && continue
                subtree2 = subtrees[j]
                issym = is_equivalent_subtree(bcg, subtree1, subtree2)
                if issym
                    push!(equiv_subtrees_subset, subtree2)
                    processed[j] = true
                end
            end
            if length(equiv_subtrees_subset) > 0
                push!(equiv_subtrees_subset, subtree1)
                push!(equivalent_subtrees, equiv_subtrees_subset)
                processed[i] = true
            end
        end
        for i in 1:m
            if !processed[i]
                push!(semiequiv_subtrees_subset, subtrees[i])
                processed[i] = true
            end
        end
        push!(semi_equivalent_subtrees, semiequiv_subtrees_subset)
    end

    return equivalent_subtrees, semi_equivalent_subtrees
end

"""
    process_equivalent_subtrees(BlockCopolymerGraph, subtrees)

For each equivalent subtree, there should be only one branch (emanating from the front vertex of that subtree) contains all other equivalent subtrees.
Therefore, we can always propagate one edge.

The special case is when there are only two equivalent subtrees, and their front vertices form an edge. In such case, we can still propagate one edge, and the front vertices of resulted two equivalent subtrees exchange their positions.
"""
function process_equivalent_subtrees(bcg::BlockCopolymerGraph, subtrees)
    equivalent_subtrees = []
    isomorphic_subtrees = []

    m = length(subtrees)

    (m == 1) && error("At least two equivalent subtrees must be provided. Only got one.")

    # Finalize equivalent trees when there are only two equivalent trees and they are connected by an edge. No more isomorphic subtrees are returned.
    if m == 2
        subtree1, subtree2 = subtrees
        v1, v2 = subtree1.v, subtree2.v
        if has_edge(bcg.graph, v1, v2)
            # Note we have to include all other branches emanating from v1 and v2 into our two equivalent trees.
            # In this case, the equivalent trees are final, meaning that they can not be expanded further.
            tree1, tree2 = induced_subtree(bcg, Edge(v1, v2))
            vertices1 = copy(tree1.vmap)
            vertices2 = copy(tree2.vmap)
            push!(vertices1, v2)
            push!(vertices2, v1)
            new_subtree1 = Subtree(bcg, vertices1, v2)
            new_subtree2 = Subtree(bcg, vertices2, v1)
            append!(equivalent_subtrees, [new_subtree1, new_subtree2])
            return equivalent_subtrees, isomorphic_subtrees
        end
    end

    front_vertices = [subtree.v for subtree in subtrees]
    for i in 1:m
        subtree = subtrees[i]
        v = subtree.v
        # induced all branches including subtree itself from front vertex v
        branches = induced_subtree(bcg, v)
        to_merge_branches = []  # store all branches to be merged
        exbranch = nothing  # the branch contains all other front vertices excpet v. This branch should not be merged.
        for branch in branches
            # a branch is not a exbranch by default
            exclude = false
            for vf in front_vertices
                # the front vertex v is in all branches, should be skipped.
                (vf == v) && continue
                # If any front vertex other than v is in the branch, we are in the exbran. Set the flag to be true and break this loop.
                (vf ∈ branch.vmap) && (exclude = true; break)
            end
            exclude ? (exbranch=branch) : push!(to_merge_branches, branch)
        end
        # Find the new front vertex in the exbranch which forms an edge with v.
        # Note that the internal front vertex (exbranch.vi) is used. And we have to convert the internal index to the original BlockCopolymerGraph.
        vnext = exbranch.vmap[first(neighbors(exbranch.graph, exbranch.vi))]
        # Merge all branches excluding the exbranch.
        merged_subtree = merge_subtrees(bcg, to_merge_branches)
        # Add the new front vertex into the merged subtree.
        new_vertices = copy(merged_subtree.vmap)
        push!(new_vertices, vnext)
        # Generate the new subtree and add it to the list.
        new_subtree = Subtree(bcg, new_vertices, vnext)
        push!(isomorphic_subtrees, new_subtree)
    end

    return equivalent_subtrees, isomorphic_subtrees
end

"""
    process_leaf(BlockCopolymerGraph, leaf_vertex)

Process the whole BlockCopolymerGraph to identify equivalent and semi-equivalent subtrees starting from a leaf vertex.
"""
function process_leaf(bcg::BlockCopolymerGraph, v)
    equivalent_subtrees = []
    semi_equivalent_subtrees = []

    vnext = find_leaf_neighbor(bcg, v)
    subtree = Subtree(bcg.graph, [v, vnext], vnext)

    n = count_isomorphic_subtree(bcg, subtree)
    (n == 1) && return [], [], [], []

    subtrees = all_isomorphic_subtree(bcg, subtree)
    isomorphic_subtrees_list = [subtrees]

    while length(isomorphic_subtrees_list) > 0
        equiv_subtrees_list = []
        equiv_vertices_list = []
        for isomorphic_subtrees in isomorphic_subtrees_list
            sg = group_isomorphic_subtrees(isomorphic_subtrees)
            es, ss, is = process_isomorphic_subtree_groups(bcg, sg)
            append!(equivalent_subtrees, es)
            append!(semi_equivalent_subtrees, ss)
            for isolated_subtrees in is
                es, ss = process_isolated_isomorphic_subtrees(bcg,
                                                isolated_subtrees)
                append!(semi_equivalent_subtrees, ss)
                append!(equiv_subtrees_list, es)
            end
        end

        isomorphic_subtrees_list = []
        for equiv_subtrees in equiv_subtrees_list
            es, is = process_equivalent_subtrees(bcg, equiv_subtrees)
            (length(es) > 0) && push!(equivalent_subtrees, es)
            (length(is) > 0) && push!(isomorphic_subtrees_list, is)
        end
    end

    return equivalent_subtrees, semi_equivalent_subtrees
end