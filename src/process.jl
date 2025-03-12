is_leaf_node(g::AbstractGraph, v) = length(neighbors(g, v)) == 1

function is_leaf_edge(g::AbstractGraph, e::AbstractEdge)
    return is_leaf_node(g, src(e)) || is_leaf_node(g, dst(e))
end

function all_leafs(g::AbstractGraph)
    return [v for v in vertices(g) if is_leaf_node(g, v)]
end

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
        push!(semi_equivalent_subtrees, subtrees)
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
        semiequiv_subtrees_subset = []
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
    (n == 1) && return [], []

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

function can_merge(g::BlockCopolymerGraph, ts1::Set{<:Subtree}, ts2::Set{<:Subtree})
    length(ts1) == length(ts2) || return false
    t1 = first(ts1)
    t1 ∈ ts2 && return true
    t2 = first(ts2)
    return is_equivalent_subtree(g, t1, t2, check=true)
end

function merge_elements(g::BlockCopolymerGraph, ts1::Set{<:Subtree}, ts2::Set{<:Subtree})
    can_merge(g, ts1, ts2) || return (ts1, ts2)

    out = deepcopy(ts1)
    for t2 in ts2
        t2 ∈ ts1 || push!(out, t2)
    end

    return (out,)
end

function all_distinct(g::BlockCopolymerGraph, vts::Vector{Set{<:Subtree}})
    isempty(vts) && return true
    length(vts) == 1 && return true

    for i in 1:length(vts)
        for j in i+1:length(vts)
            can_merge(g, vts[i], vts[j]) && return false
        end
    end

    return true
end

function merge_elements(g::BlockCopolymerGraph, vts::Vector{Set{<:Subtree}})
    (length(vts) <= 1) && return vts

    temp = Set{<:Subtree}[]
    out = deepcopy(vts)
    merged = fill(false, length(vts))
    while !all_distinct(g, out)
        for i in 1:length(out)
            merged[i] && continue
            for j in i+1:length(out)
                merged[j] && continue
                res = merge_elements(g, out[i], out[j])
                if length(res) == 1
                    push!(temp, res[1])
                    merged[i] = true
                    merged[j] = true
                end
            end
        end
        for i in 1:length(out)
            !merged[i] && push!(temp, vts[i])
        end
        empty!(out)
        append!(out, temp)
        empty!(temp)
        merged = fill(false, length(out))
    end

    return out
end

function find_all_equivalent_subtrees(g::BlockCopolymerGraph)
    leafs = all_leafs(g)
    equivalent_trees = Set{<:Subtree}[]
    for leaf in leafs
        es, ss = process_leaf(g, leaf)
        isempty(es) && continue
        for etrees in es
            isempty(etrees) && continue
            push!(equivalent_trees, Set([t for t in etrees]))
        end
    end

    equivalent_trees = merge_elements(g, equivalent_trees)

    return equivalent_trees
end

find_all_equivalent_subtrees(p::BlockCopolymer) = find_all_equivalent_subtrees(BlockCopolymerGraph(p))

function list_edges_by_order!(edges::Vector{<:Edge}, st::Subtree, v, vp)
    for vx in neighbors(st.graph, v)
        (vx == vp) && continue
        e = Edge(st.vmap[vx], st.vmap[v])
        push!(edges, e)
        is_leaf_node(st.graph, vx) || list_edges_by_order!(edges, st, vx, v)
    end
end

"""
    list_edges_by_order(st::Subtree)

List all edges in a subtree via a depth-first traversal starting from the front vertex of the input subtree. The return value is a vector of `Edge` instance. The node values in the returned edges are in the convention of original BlockCopolymerGraph.  
"""
function list_edges_by_order(st::Subtree)
    edges = Edge[]
    list_edges_by_order!(edges, st, st.vi, nothing)
    return edges
end

"""
    vertex_mapping(g::BlockCopolymerGraph, es1::Subtree, es2::Subtree)

Return a mapping of nodes for two equivalent subtrees. Since equivalent subtrees are isomorphic. Given two subtrees G and H, this method identify a mapping f: u -> v, where u is a node in G and v is a node in v.
"""
function vertex_mapping(g::BlockCopolymerGraph, es1::Subtree, es2::Subtree)
    vertex_relation = (v1, v2) -> equivalent_blockend(g, v1, v2, es1.vmap, es2.vmap)
    edge_relation = (e1, e2) -> equivalent_block(g, e1, e2, es1.vmap, es2.vmap)
    ichannel = all_isomorph(es1.graph, es2.graph; vertex_relation, edge_relation)
    internal_mapping = collect(ichannel)[1]
    T = eltype(g)
    vmap = Dict{T,T}()
    for (v1, v2) in internal_mapping
        vmap[es1.vmap[v1]] = es2.vmap[v2]
    end
    return vmap
end

"""
    apply_mapping(vmap, edges::Vector{<:Edge})

Apply the mapping identified by the method `vertex_mapping` to convert a list of edges in a subtree to another list of edges in its equivalent subtree. The input argument `vmap` is computed from `vertex_mapping` using subtree1 and subtree2.
"""
function apply_mapping(vmap, edges::Vector{<:Edge})
    return [Edge(vmap[e.src], vmap[e.dst]) for e in edges]
end

"""
    unique_blocks(eqblock_list::Vector{Vector{Pair{V,V}}}) where V

Remove identical equivalent blocks in the list. This case occurs when two equivalent subtrees have a common edge. For example, if the common edge is Edge(7,8), it will produce [[7=>8, 8=>7], [8=>7, 7=>8]]. Therefore, we have to remove one of these to eliminate duplication.
"""
function unique_blocks(eqblock_list::Vector{Vector{Pair{V,V}}}) where V
    same_blocks = Pair{V, V}[]
    unique_eqblock_list = Vector{Pair{V,V}}[]
    for eqblocks in eqblock_list
        if length(eqblocks) != 2
            push!(unique_eqblock_list, eqblocks)
            continue
        end
        b1, b2 = eqblocks
        if (b1 == reverse(b2))
            if b1 ∉ same_blocks
                push!(same_blocks, b1)
                push!(same_blocks, b2)
                push!(unique_eqblock_list, eqblocks)
            end
        else
            push!(unique_eqblock_list, eqblocks)
        end
    end

    return unique_eqblock_list
end

"""
    group_equivalent_blocks(bc::BlockCopolymer)

Group equivalent blocks into a vector of block groups. For blocks in the same group, they share a single propagator. The typical usage of this method is to reduce the computational cost of SCFT calculations.
"""
#=
function group_equivalent_blocks(bcg::BlockCopolymerGraph)
    if nv(bcg) == 2
        v1, v2 = vertices(bcg)
        return [[v1=>v2, v2=>v1]]
    end

    V = eltype(bcg)
    es_list = find_all_equivalent_subtrees(bcg)
    eqblock_list = Vector{Pair{V,V}}[]
    for es in es_list
        ess = collect(es)
        idx0 = length(eqblock_list) + 1
        edges = list_edges_by_order(ess[1])
        for e in edges
            v1, v2 = src(e), dst(e)
            push!(eqblock_list, [v1=>v2])
            push!(eqblock_list, [v2=>v1])
        end
        length(ess) == 1 && continue
        for et in ess[2:end]
            idx = idx0
            vmap = vertex_mapping(bcg, ess[1], et)
            edgesx = apply_mapping(vmap, edges)
            for e in edgesx
                v1, v2 = src(e), dst(e)
                push!(eqblock_list[idx], v1=>v2)
                push!(eqblock_list[idx+1], v2=>v1)
                idx += 2 
            end
        end
    end

    eqblock_list = unique_blocks(eqblock_list)
    flat_eqblock_list = isempty(eqblock_list) ? [] : reduce(vcat, eqblock_list)
    for e in keys(bcg.edge2block)
        v1, v2 = e
        if (v1=>v2) ∉ flat_eqblock_list
            push!(eqblock_list, [v1=>v2])
            push!(eqblock_list, [v2=>v1])
        end
    end

    return eqblock_list
end
=#

function bfs_path(graph, start, goal; print_level = false)
    predecessor = Dict{Int, Int}() #前驱顶点
    visited = Set{Int}() #已访问顶点
    queue = [start] #顶点队列
    push!(visited, start)
    
    level = 0  # 当前搜索层次
    print_level && println("Level $level: $queue")  # 当前层次的顶点
    
    while !isempty(queue)
        level += 1
        current_size = length(queue)
        for i in 1:current_size
            current = popfirst!(queue) # 移除并获取第一个顶点 current
            if current == goal
                break
            end
            for neighbor in neighbors(graph, current)
                if !(neighbor in visited)
                    push!(queue, neighbor)
                    push!(visited, neighbor)
                    predecessor[neighbor] = current # 前驱顶点为 current
                end
            end
        end
        print_level && println("Level $level: $queue")  # 输出当前层次的节点
    end
    
    # 重建路径
    if !(haskey(predecessor, goal))
        return nothing
    end
    
    path = [goal]
    while path[end] != start
        push!(path, predecessor[path[end]])
    end
    reverse(path)
end

function dependency_list(graph, pair::Pair{Int, Int})
    src, dst = pair[1], pair[2]
    list = Pair{Int, Int}[]
    for neighbor in neighbors(graph, src)
        !(neighbor == dst) && push!(list, neighbor=>src)
    end
    return list
end

function equiv_pair(graph, pair1, pair2)
    block1 = graph.edge2block[Polymer._sort_tuple2((pair1[1], pair1[2]))]
    block2 = graph.edge2block[Polymer._sort_tuple2((pair2[1], pair2[2]))]
    if block1.segment == block2.segment && block1.f == block2.f
        return true
    end
    return false
end

function is_equivalent_dependency_sets(graph, list1, list2, equivalent_list)
    n1, n2 = length(list1), length(list2)
    correct_matrix = fill(-1, n1, n2)
    list1 = collect(list1)
    list2 = collect(list2)

    for index1 in 1:n1
        for index2 in 1:n2
            test_set = Set([list1[index1], list2[index2]])
            for equiv_set in equivalent_list
                if issubset(test_set, equiv_set) || length(test_set) == 1
                    if correct_matrix[index1, index2] == -1
                        correct_matrix[index1, index2] = 1
                        correct_matrix[index1, :] .= replace(correct_matrix[index1, :], -1 => 0)
                        correct_matrix[:, index2] .= replace(correct_matrix[:, index2], -1 => 0)
                        break
                    end
                end
            end
        end
    end

    if sum(correct_matrix) == n1
        return true
    end
    return false
end

function is_equivalent_blocks(graph, pair1, pair2, equivalent_list)
    list1 = dependency_list(graph, pair1)
    list2 = dependency_list(graph, pair2)
    if length(list1) == length(list2) && equiv_pair(graph, pair1, pair2)
        if isempty(list1)
            return true
        elseif is_equivalent_dependency_sets(graph, list1, list2, equivalent_list)
            return true
        end
    end
    return false
end

function find_computation_sequence(graph, goal_vertex)    
    forward_paths = []
    backward_paths = []
    all_blocks = []
    visited = Set()
    
    for edge in collect(edges(graph))
        v1, v2 = edge.src, edge.dst
        push!(all_blocks, v1=>v2)
        push!(all_blocks, v2=>v1)
    end
    
    nblocks = div(length(all_blocks), 2)
    
    for leaf_vertex in all_leafs(graph)
        forward_path = bfs_path(graph, leaf_vertex, goal_vertex;)
        backward_path = bfs_path(graph, goal_vertex, leaf_vertex;)
        !isnothing(forward_path) && push!(forward_paths, forward_path)
        !isnothing(backward_path) && push!(backward_paths, backward_path)
    end
    
    index_forward = fill(1,length(forward_paths))
    index_backward = fill(1,length(backward_paths))
    
    sequence = Dict{Int, Set{Pair{Int, Int}}}()
    layer = 1
    
    while length(visited) < nblocks && layer <= nblocks
        for i in 1:length(forward_paths)
            path = forward_paths[i]
            if index_forward[i] < length(path)
                pair = path[index_forward[i]] => path[index_forward[i]+1]
                dep_set = Set(dependency_list(graph, pair))
                if isempty(dep_set) || dep_set ∩ visited == dep_set
                    push!(visited, pair)
                    if !haskey(sequence, layer)
                        sequence[layer] = Set{Pair{Int, Int}}()
                    end
                    push!(sequence[layer], pair)
                    index_forward[i] += 1
                end
            end
        end
        layer += 1
    end
    
    while length(visited) < nblocks * 2 && layer <= nblocks * 2
        for i in 1:length(backward_paths)
            path = backward_paths[i]
            if index_backward[i] < length(path)
                pair = path[index_backward[i]] => path[index_backward[i]+1]
                dep_set = Set(dependency_list(graph, pair))
                if isempty(dep_set) || dep_set ∩ visited == dep_set
                    push!(visited, pair)
                    if !haskey(sequence, layer)
                        sequence[layer] = Set{Pair{Int, Int}}()
                    end
                    push!(sequence[layer], pair)
                    index_backward[i] += 1
                end
            end
        end
        layer += 1
    end
    
    output_sequence = Dict(key => collect(value) for (key, value) in sequence)
    return all_blocks, visited, output_sequence
end

function find_equiv_blocks(graph, goal_vertex)
    all_blocks, visited, sequence = find_computation_sequence(graph, goal_vertex)  
    equivalent_list = []

    for seq_index in sort(collect(keys(sequence)))
        blocks = collect(sequence[seq_index])
        for i in 1:(length(blocks)-1)
            for j in i+1:length(blocks)
                pair1 = blocks[i]
                pair2 = blocks[j]
                if is_equivalent_blocks(graph, pair1, pair2, equivalent_list)
                    push!(equivalent_list, Set([pair1, pair2]))
                    visited = filter!(!=(pair1), visited)
                    visited = filter!(!=(pair2), visited)
                end
            end
        end
    end
    
    for i in 1:(length(visited)-1)
        for j in i+1:length(visited)
            pair1 = collect(visited)[i]
            pair2 = collect(visited)[j]
            if is_equivalent_blocks(graph, pair1, pair2, equivalent_list)
                push!(equivalent_list, Set([pair1, pair2]))
            end
        end
    end

    merged_sets = []

    for current_set in equivalent_list
        merged = false
        for i in 1:length(merged_sets)
            if !isempty(current_set ∩ merged_sets[i])
                merged_sets[i] = merged_sets[i] ∪ current_set
                merged = true
                break
            end
        end
        if !merged
            push!(merged_sets, current_set)
        end
    end
    
    all_equiv_set = Set()
    for set in merged_sets
        all_equiv_set = all_equiv_set ∪ set
    end
    
    for pair in all_blocks
        if !((pair) in all_equiv_set)
            push!(merged_sets, Set([pair]))
        end
    end
    
    output_vector = [collect(set) for set in merged_sets]
    return output_vector 
    # n-element Vector{Vector{Pair{Int64, Int64}}}
end

function find_computation_sequence(graph)    
    groups_dict = Dict()
    subtree_arrays = [collect(set) for set in find_all_equivalent_subtrees(graph)]
    front_vertices = Set()
    for subtree_array in subtree_arrays
        front_vertex = subtree_array[1].v
        push!(front_vertices, front_vertex)
    end
    front_vertices = collect(front_vertices)
    
    for goal in front_vertices
        groups_dict[goal] = length(find_equiv_blocks(graph, goal))
    end

    sequence = find_computation_sequence(graph, argmin(groups_dict))[3]

    return sequence
end

function group_equiv_blocks(graph)
    groups_dict = Dict()
    front_vertices = collect(values(graph.joint2node))

    for goal in front_vertices
        groups_dict[goal] = length(find_equiv_blocks(graph, goal))
    end

    groups = find_equiv_blocks(graph, argmin(groups_dict))

    return groups
end
group_equiv_blocks(bc::BlockCopolymer) = group_equiv_blocks(BlockCopolymerGraph(bc))

group_equivalent_blocks(bcg::BlockCopolymerGraph) = group_equiv_blocks(bcg::BlockCopolymerGraph)
group_equivalent_blocks(bc::BlockCopolymer) = group_equiv_blocks(bc::BlockCopolymer)