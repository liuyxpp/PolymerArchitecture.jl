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