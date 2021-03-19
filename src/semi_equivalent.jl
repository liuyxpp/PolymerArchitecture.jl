"""
    group_isomorphic_branches(BlockCopolymerGraph, branches)

Group all isomorphic branches detached from a common vertex. Note that this is different from `group_isomorphic_subtrees`.
"""
function group_isomorphic_branches(bcg::BlockCopolymerGraph, branches)
    branch_groups = []

    n = length(branches)
    (n == 1) && return [branches]

    grouped = [false for _ in 1:n]
    for i in 1:n
        grouped[i] && continue
        subtree1 = branches[i]
        bgroup = [subtree1]
        for j in (i+1):n
            grouped[j] && continue
            subtree2 = branches[j]
            if is_isomorphic_subtree(bcg, subtree1, subtree2)
                push!(bgroup, subtree2)
                grouped[j] = true
            end
        end
        grouped[i] = true
        push!(branch_groups, bgroup)
    end

    return branch_groups
end

"""
    group_semi_equivalent_subtrees(BlockCopolymerGraph, subtrees)

Group semi-equivalent subtrees according to the number of branches emanating from the front vertices.
"""
function group_semi_equivalent_subtrees(bcg::BlockCopolymerGraph, subtrees)
    tree_groups = []

    nbranch = [degree(bcg.graph, i) for i in front_vertices(subtrees)]
    nbranch_unique = nbranch |> unique |> sort

    for n in nbranch_unique
        tgroup = []
        for i in 1:length(subtrees)
            (nbranch[i] == n) && push!(tgroup, subtrees[i])
        end
        push!(tree_groups, tgroup)
    end

    return tree_groups
end

"""
    process_semi_equivalent_subtree_group(BlockCopolymerGraph, subtrees)
    process_semi_equivalent_subtree_group(::Val{1}, BlockCopolymerGraph, subtrees)
    process_semi_equivalent_subtree_group(::Val{2}, BlockCopolymerGraph, subtrees)
    process_semi_equivalent_subtree_group(::Val{N}, BlockCopolymerGraph, subtrees) where N

Process a group of semi-equivalent subtrees with same number of branches on their front vertices. A list of semi-equivalent subtrees and a list of isomorphic subtrees are returned. The list may be empty. The semi-equivalent subtrees returned are final, meaning that they will not be processed further.
"""
process_semi_equivalent_subtree_group

"""
When there is only one subtree in the group, no process is needed and simply return it as a final semi-equivalent subtree. Note that, to make this subtree useful, we should combine it with other equivalent or semi-equivalent subtrees in the same leaf process.
"""
function process_semi_equivalent_subtree_group(::Val{1}, bcg::BlockCopolymerGraph, subtrees)
    return [subtrees], [front_vertices], [], []
end

"""
When there are two subtrees in the group, we have to perform following processing:

1. 
"""
function process_semi_equivalent_subtree_group(::Val{2}, bcg::BlockCopolymerGraph, subtrees)
    semi_equivalent_subtrees = []
    isomorphic_subtrees = []

    subtree1, subtree2 = deepcopy(subtrees)
    v1, v2 = subtree1.v, subtree2.v
    # all subtrees have same number of branches
    nbranch = degree(bcg.graph, v1)

    if nbranch == 2
        # find the edge ouside current subtree but contains the front vertex.
        e1, e2 = nothing, nothing
        for vx in neighbors(bcg.graph, v1)
            (vx ∉ subtree1.vmap) && (e1 = Edge(vx, v1))
        end
        for vx in neighbors(bcg.graph, v2)
            (vx ∉ subtree2.vmap) && (e2 = Edge(vx, v2))
        end
        if equivalent_block(bcg, e1, e2)
            vertices1 = [src(e1)]
            vertices2 = [src(e2)]
            push!(vertices1, subtree1.vmap)
            push!(vertices2, subtree2.vmap)
            new_subtree1 = Subtree(bcg, vertices1, src(e1))
            new_subtree2 = Subtree(bcg, vertices2, src(e2))
            push!(isomorphic_subtrees, [new_subtree1, new_subtree2])
        else
            push!(semi_equivalent_subtrees, subtrees)
        end
        return semi_equivalent_subtrees, isomorphic_subtrees
    end

    # when nbranch > 2
    # Find all  that is in the subtree but not the front vertex.
    v1pre = subtree1.vmap[subtree1.vi == 1 ? 2 : subtree1.vi - 1]
    v2pre = subtree2.vmap[subtree2.vi == 1 ? 2 : subtree2.vi - 1]
    isomorphic_branches1 = []
    isomorphic_branches2 = []
    branch12 = nothing  # branch detached from v1 containing v2
    branch21 = nothing  # branch detached from v2 containing v1
    branches1 = induced_subtree(bcg, v1)
    newbranches1 = []
    for branch in branches1
        # skip the subtree1 branch
        (v1pre ∈ branch.vmap) && continue
        # skip the branch that contains subtree2.
        (v2 ∈ branch.vmap) && (branch12 = branch; continue)
        push!(newbranches1, branch)
    end
    branches2 = induced_subtree(bcg, v2)
    newbranches2 = []
    for branch in branches2
        # skip the subtree2 branch
        (v2pre ∈ branch.vmap) && continue
        # skip the branch that contains subtree1.
        (v1 ∈ branch.vmap) && (branch21 = branch; continue)
        push!(newbranches2, branch)
    end
    bgroups1 = group_isomorphic_branches(bcg, newbranches1)
    bgroups2 = group_isomorphic_branches(bcg, newbranches2)
    for i in 1:length(bgroups1)
        bgroup1 = bgroups1[i]
        b1 = first(bgroup1)
        for j in 1:length(bgroups2)
            bgroup2 = bgroups2[j]
            b2 = first(bgroup2)
            if is_isomorphic_tree(first(bgroup1), v1, first(bgroup2), v2, bcg)
                # Only including minimum number of isomorphic branches
                p = min(length(bgroup1), length(bgroup2))
                append!(isomorphic_branches1, bgroup1[1:p])
                append!(isomorphic_branches2, bgroup2[1:p])
                # Since we have grouped all isomorphic branches together,
                # there can be only one group in v1 that is isomorphic to v2
                break
            end
        end
    end
    # merge all isomorphic branches into subtree
    tree1 = copy(subtree1)
    for ib1 in isomorphic_branches1
        append!(tree1, ib1)
    end
    # remove multiple v1
    tree1 = unique(tree1)
    tree2 = copy(subtree2)
    for ib2 in isomorphic_branches2
        append!(tree2, ib2)
    end
    # remove multiple v2
    tree2 = unique(tree2)
    if length(isomorphic_branches1) == length(isomorphic_branches2) && length(isomorphic_branches1) == length(newbranches1)
        # all branches are included, we will then step an edge
        for vx in neighbors(bcg.graph, v1)
            (vx ∈ branch12) && (e1 = Edge(vx, v1); break)
        end
        for vx in neighbors(bcg.graph, v2)
            (vx ∈ branch21) && (e2 = Edge(vx, v2); break)
        end
        if equivalent_block(e1, e2, bcg)
            push!(tree1, src(e1))
            push!(tree2, src(e2))
            push!(isomorphic_subtrees, [tree1, tree2])
            push!(isomorphic_front_vertices, [src(e1), src(e2)])
        else
            push!(semi_equivalent_subtrees, [tree1, tree2])
            push!(semi_equivalent_vertices, [v1, v2])
        end
    else
        push!(semi_equivalent_subtrees, [tree1, tree2])
        push!(semi_equivalent_vertices, [v1, v2])
    end
    return semi_equivalent_subtrees, isomorphic_subtrees
end