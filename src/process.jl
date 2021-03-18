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
        isomorphic_vertices = copy(subtree.vmap)
        new_front_vertex = nothing
        for vx in neighbors(bcg.graph, v)
            (vx ∈ subtree.vmap) && continue
            subtreex, _ = induced_subtree(bcg, Edge(vx, v))
            tomerge = true
            for vi in front_vertices
                # we are in the branch containing all other equivalent subtrees
                (vi ∈ subtreex.vmap) && (tomerge = false; break)
            end
            tomerge ? append!(isomorphic_vertices, subtree.vmap) : new_front_vertex = vx
        end
        push!(isomorphic_vertices, new_front_vertex)
        new_subtree = Subtree(bcg, isomorphic_vertices, new_front_vertex)
        push!(isomorphic_subtrees, new_subtree)
    end

    return equivalent_subtrees, isomorphic_subtrees
end