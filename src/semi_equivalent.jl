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
    _find_branches(BlockCopolymerGraph, subtree1, subtree2)

Return a list of branches that are detached from subtree1.v, excluding all branches which contain vertices in either subtree1 or subtree2, and also return the branch contains the whole subtree2.
"""
function _find_branches(bcg::BlockCopolymerGraph, subtree1, subtree2)
    v1, v2 = subtree1.v, subtree2.v
    # Find all vertices which form edges with the front vertex.
    vspre1 = pre_front_vertices(subtree1)
    branches = induced_subtree(bcg, v1)
    newbranches = []
    branch12 = nothing
    for branch in branches
        # skip the subtree1 branches
        # Note that subtree1 can be broken into several branches.
        branch_in_subtree1 = false
        for v in vspre1
            (v ∈ branch.vmap) && (branch_in_subtree1 = true; break)
        end
        branch_in_subtree1 && continue
        # skip the branch that contains subtree2.
        (v2 ∈ branch.vmap) && (branch12 = branch; continue)
        push!(newbranches, branch)
    end

    return newbranches, branch12
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

1. when number of branches = 2
2. when number of branches > 2
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
    isomorphic_branches1 = [subtree1]
    isomorphic_branches2 = [subtree2]
    # branch12: branch detached from v1 containing v2
    newbranches1, branch12 = _find_branches(bcg, subtree1, subtree2)
    # branch21: branch detached from v2 containing v1
    newbranches2, branch21 = _find_branches(bcg, subtree2, subtree1)
    bgroups1 = group_isomorphic_branches(bcg, newbranches1)
    bgroups2 = group_isomorphic_branches(bcg, newbranches2)
    for i in 1:length(bgroups1)
        bgroup1 = bgroups1[i]
        for j in 1:length(bgroups2)
            bgroup2 = bgroups2[j]
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
    tree1 = merge_subtrees(bcg, isomorphic_branches1)
    tree2 = merge_subtrees(bcg, isomorphic_branches2)
    # isomorphic_branches1 contains subtree1 but newbranches1 not, thus +1 here to check.
    if length(isomorphic_branches1) == length(newbranches1) + 1
        # all branches are included, we will then step an edge
        for vx in neighbors(bcg.graph, v1)
            (vx ∈ branch12.vmap) && (e1 = Edge(vx, v1); break)
        end
        for vx in neighbors(bcg.graph, v2)
            (vx ∈ branch21.vmap) && (e2 = Edge(vx, v2); break)
        end
        if equivalent_block(bcg, e1, e2)
            # Add the new front vertex into the merged subtree.
            new_vertices1 = copy(tree1.vmap)
            new_vertices2 = copy(tree2.vmap)
            push!(new_vertices1, src(e1))
            push!(new_vertices2, src(e2))
            # Generate the new subtree and add it to the list.
            newtree1 = Subtree(bcg, new_vertices1, src(e1))
            newtree2 = Subtree(bcg, new_vertices2, src(e2))
            push!(isomorphic_subtrees, [newtree1, newtree2])
        else
            push!(semi_equivalent_subtrees, [tree1, tree2])
        end
    else
        push!(semi_equivalent_subtrees, [tree1, tree2])
    end
    return semi_equivalent_subtrees, isomorphic_subtrees
end

"""
This is for N-element group where N > 2. There are N subtrees in the group.
"""
function process_semi_equivalent_subtree_group(::Val{N}, bcg::BlockCopolymerGraph, subtrees) where N
    semi_equivalent_subtrees = []
    isomorphic_subtrees = []

    for i in 1:length(subtrees)
        subtree1 = subtrees[i]
        for j in (i+1):length(subtrees)
            subtree2 = subtrees[j]
            ss, is = process_semi_equivalent_subtree_group(Val(2), bcg, [subtree1, subtree2])
            (length(ss) > 0) && append!(semi_equivalent_subtrees, ss)
            (length(is) > 0) && append!(isomorphic_subtrees, is)
        end
    end

    return semi_equivalent_subtrees, isomorphic_subtrees
end

function process_semi_equivalent_subtree_group(bcg::BlockCopolymerGraph, subtrees)
    m = length(subtrees)
    # Simply dispatch accroding to the number of subtrees in the group
    return process_semi_equivalent_subtree_group(Val(m), bcg, subtrees)
end

"""
    process_semi_equivalent_subtrees(BlockCopolymerGraph, subtrees)

Process semi-equivalent subtrees which are obtained from `process_leaf`. After processing, the size of each semi-equivalent subtree should be maximized.
"""
function process_semi_equivalent_subtrees(bcg::BlockCopolymerGraph, subtrees)
    semi_equivalent_subtrees = []

    isomorphic_subtrees = [subtrees]
    while length(isomorphic_subtrees) > 0
        tree_groups = []
        for i in 1:length(isomorphic_subtrees)
            trees = isomorphic_subtrees[i]
            tgroups = group_semi_equivalent_subtrees(bcg, trees)
            append!(tree_groups, tgroups)
        end

        isomorphic_subtrees = []
        for i in 1:length(tree_groups)
            ss, is = process_semi_equivalent_subtree_group(bcg, tree_groups[i])
            (length(ss) > 0) && append!(semi_equivalent_subtrees, ss)
            (length(is) > 0) && append!(isomorphic_subtrees, is)
        end
    end

    return semi_equivalent_subtrees
end
