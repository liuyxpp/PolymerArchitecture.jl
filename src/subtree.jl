abstract type AbstractSubtree end

"""
    Subtree

An induced subgraph of the original graph. It has a special vertex which is named *front vertex*. A `Subtree` instance can be used to represent the sub-architecture of a non-cyclic block copolymer chain.

## Fields
* `graph`: a LightGraphs object that describes the subtree.
* `vmap`: an `AbstractVector` object whose indices are vertices of `Subtree.graph` and values are their corresponding vertices in the original graph.
* `vmap`: an `AbstractVector` object which is a reverse of `vmap`.
* `vi`: the internal index (i.e. vertex in `Subtree.graph`) of the front vertex. Value 0 means the front vertex is not set.
* `v`: the index of the front vertex in the original graph. Value 0 means the front vertex is not set.
"""
struct Subtree{T, G<:AbstractGraph{T}, M<:AbstractVector{T}} <: AbstractSubtree
    graph::G
    vmap::M
    rvmap::Dict{T,T}
    vi::T
    v::T
end

function _subtree(graph::AbstractGraph{T}, vmap::AbstractVector{T}, v::T=zero(T)) where T
    rvmap = Dict{T,T}(vmap[i] => i for i in 1:length(vmap))
    vi = (v == 0) ? zero(T) : rvmap[v]
    return Subtree(graph, vmap, rvmap, vi, v)
end

function Subtree(graph::AbstractGraph{T}, vertices, v::T=zero(T)) where T
    g, vm = induced_subgraph(graph, vertices)
    return _subtree(g, vm, v)
end

Subtree(bcg::BlockCopolymerGraph, vertices, v=0) = Subtree(bcg.graph, vertices, v)

"""
    front_vertex(subtree)

Return the front vertex of a Subtree object. Note that the vertex is indexed in the orginal graph (the graph that the subtree is induced).
"""
front_vertex(subtree::Subtree) = subtree.v

"""
    front_vertices(subtrees)

Return a list of front vertices corresponding to the list of subtrees.
"""
front_vertices(subtrees) = [subtree.v for subtree in subtrees]

"""
    front_edge(subtree)

Return front edges of a Subtree object. A front edge is an edge in the subtree that contains the front vertex. Note that all edges are indexed in the original graph (the graph that the subtree is induced).
"""
function front_edges(subtree::Subtree)
    edges = Edge{typeof(subtree.v)}[]
    for vx in neighbors(subtree.graph, subtree.vi)
        push!(edges, Edge(subtree.v, subtree.vmap[vx]))
    end
    return edges
end

"""
    find_leaf_neighbor(BlockCopolymerGraph, v)

Return the vertex which forms an edge with the leaf vertex. Note that there should be only one neighbor of the leaf vertex.
"""
find_leaf_neighbor(bcg::BlockCopolymerGraph, v) = first(neighbors(bcg.graph, v))

"""
    find_neighbors(vertices, graph, v, vfrom)

Find those neighbors and all their neighborhood that are connected to vertex v, excluding the one vfrom and its all neighborhood.

## Arguments
* `vertices`: contains all vertices of the neighbors found.
* `graph`: find the neighbors in this graph.
* `v`: neighbors of this vertex to be found.
* `vfrom`: (v, vfrom) is an edge of `graph`.
"""
function find_neighbors!(vertices::AbstractVector{T}, graph::AbstractGraph{T}, v::T, vfrom::T) where T
    # Do nothing when (v, vfrom) is not an edge in graph.
    has_edge(graph, v, vfrom) || return vertices

    for vx in neighbors(graph, v)
        (vx == vfrom) && continue
        push!(vertices, vx)
        if degree(graph, vx) == 1
            continue
        else
            find_neighbors!(vertices, graph, vx, v)
        end
    end

    return vertices
end

function find_neighbors(graph::AbstractGraph{T}, v::T, vfrom::T) where T
    vertices = T[]
    find_neighbors!(vertices, graph, v, vfrom)
    return vertices
end

"""
    induced_subtree(BlockCopolymerGraph, edge)
    induced_subtree(LightGraphs, edge)

Make two induced subgraphs of the input graph by cutting the edge e.
Note that the input edge is not included in the induced subgraphs.
"""
function induced_subtree(graph::AbstractGraph, e::AbstractEdge)
    v1, v2 = src(e), dst(e)
    vertices1 = find_neighbors(graph, v1, v2)
    push!(vertices1, v1)
    vertices2 = find_neighbors(graph, v2, v1)
    push!(vertices2, v2)

    subtree1 = Subtree(graph, vertices1, v1)
    subtree2 = Subtree(graph, vertices2, v2)

    return subtree1, subtree2
end

induced_subtree(bcg::BlockCopolymerGraph, e::AbstractEdge) = induced_subtree(bcg.graph, e)

"""
    induced_subtree(BlockCopolymerGraph, vertex)
    induced_subtree(LightGraphs, vertex)

Make induced subgraphs of the input graph by disconnecting all branches from vertex `v`.
Note that vertex `v` is included in all induced subtrees.
"""
function induced_subtree(graph::AbstractGraph, v)
    subtrees = Subtree[]
    for vx in neighbors(graph, v)
        vertices = find_neighbors(graph, vx, v)
        append!(vertices, [v, vx])
        subtree = Subtree(graph, vertices, v)
        push!(subtrees, subtree)
    end
    return subtrees
end

induced_subtree(bcg::BlockCopolymerGraph, v) = induced_subtree(bcg.graph, v)

"""
    induced_subtree(BlockCopolymerGraph, edge1, edge2)
    induced_subtree(LightGraphs, edge1, edge2)

Cut the input graph at edges `e1` and `e2`.
Return the subtree wich contains both edges' vertices.
"""
function induced_subtree(graph::AbstractGraph, e1::AbstractEdge, e2::AbstractEdge)
    # Make sure e1 and e2 has are two distinct edges
    length(unique([src(e1), dst(e1), src(e2), dst(e2)])) < 3 && return nothing

    # Make two subgraphs by cutting edge e1
    subtree1_src, subtree1_dst = induced_subtree(graph, e1)

    # Find out edge e2 is in which subtree1 (src or dst).
    # index of vmap1 is vertex of subtree1
    # vlaue of vmap1 is vertex of graph
    issrc = src(e2) ∈ subtree1_src.vmap
    subtree1 = issrc ? subtree1_src : subtree1_dst
    v1 = issrc ? src(e1) : dst(e1)

    # Make two subgraphs of subtree1 by cutting edge e2
    # Note that we have to map e2 of bcg to subtree1 edge.
    # key of rvmap1 is vertex of graph
    # value of rvmap1 is vertex of subtree1
    rvmap1 = subtree1.rvmap
    subtree2_src, subtree2_dst = induced_subtree(subtree1.graph, Edge(rvmap1[src(e2)], rvmap1[dst(e2)]))
    # index of vmap2_* is vertex of subgraph of subtree1
    # value of vmap2_* is vertex of subtree1

    # Find out v1 from edge e1 is in which subtree2 (src or dst)
    issrc = rvmap1[subtree1.v] ∈ subtree2_src.vmap
    subtree2 = issrc ? subtree2_src : subtree2_dst
    v2 = issrc ? src(e2) : dst(e2)

    # Remaping from vertex of subtree1 to orginal graph
    vertices = [subtree1.vmap[v] for v in subtree2.vmap]
    subtree = Subtree(graph, vertices, v2)

    return subtree, v1, v2
end

induced_subtree(bcg::BlockCopolymerGraph, e1::AbstractEdge, e2::AbstractEdge) = induced_subtree(bcg.graph, e1, e2)

"""
    induced_subtree(BlockCopolymerGraph, subtree1, subtree2)
    induced_subtree(graph, subtree1, subtree2)

Return an induced subtree whose vertices are not in `subtree1` and `subtree2` except their front vertices. The front vertex of `subtree1` is chosen to be the front vertex of the returned subtree.

`graph` should be the original graph that `subtree1` and `subtree2` are both induced from.
"""
function induced_subtree(graph::AbstractGraph, subtree1::Subtree, subtree2::Subtree)
    v1, v2 = subtree1.v, subtree2.v
    vs = setdiff(vertices(graph), subtree1.vmap ∪ subtree2.vmap)
    (v1 != 0) && push!(vs, v1)
    (v2 != 0) && push!(vs, v2)
    return Subtree(graph, vs, v1)
end

induced_subtree(bcg::BlockCopolymerGraph, subtree1::Subtree, subtree2::Subtree) = induced_subtree(bcg.graph, subtree1, subtree2)

function merge_subtrees(bcg::BlockCopolymerGraph, subtrees)
    first_tree = first(subtrees)
    vs = typeof(first_tree.v)[]
    for tree in subtrees
        append!(vs, tree.vmap)
    end
    return Subtree(bcg, unique(vs), first_tree.v)
end