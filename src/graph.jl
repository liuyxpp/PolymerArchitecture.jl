import Graphs: AbstractGraph

"""
    BlockCopolymerGraph <: AbstractGraph{Int}

A graph representation of a `Polymer.BlockCopolymerGraph` object by treating each block end as a vetex and each block as an edge. The block length is assigned to the weights for corresponding edge. see [`build_graph`](@ref) for more details.
"""
struct BlockCopolymerGraph <: AbstractGraph{Int}
    graph::SimpleGraph
    block2edge::Dict{PolymerBlock, Tuple{Int, Int}}
    edge2block::Dict{Tuple{Int, Int}, PolymerBlock}
    joint2node::Dict{BranchPoint, Int}
    node2joint::Dict{Int, BranchPoint}
    free2node::Dict{FreeEnd, Int}
    node2free::Dict{Int, FreeEnd}
    distmx::Matrix{<:Real}

    function BlockCopolymerGraph(c::BlockCopolymer)
        g, d1, d2, d3 = build_graph(c)
        rd1 = reverse_dict(d1)
        rd2 = reverse_dict(d2)
        rd3 = reverse_dict(d3)
        distmx = zeros(typeof(first(c.blocks).f), nv(g), nv(g))
        for (e, block) in rd1
            distmx[e[1], e[2]] = block.f
            distmx[e[2], e[1]] = block.f
        end
        new(g, d1, rd1, d2, rd2, d3, rd3, distmx)
    end
end

Polymer.species(bcg::BlockCopolymerGraph) = [specie(b) for b in keys(bcg.block2edge)] |> unique |> sort

"""
    build_graph(c::BlockCopolymer)

Build the graph representation of a BlockCopolymer instance, i.e. a block polymer chain. A LightGraphs undirected Graph instance as well as two Dicts are returned.

The first Dict instance is a collection of `PolymerBlock => (Int, Int)` pair where each integer number is the vertex index in the Graph instance. Note that the order of these two integers in the tuple is important. The first one corresponds to the block E1 end and the second one to the block E2 end, respectively. By denfinition, it is a one-to-one mapping. Therefore, we can safely reverse this Dict to get another Dict of `(Int, Int) => PolymerBlock`.

The Second Dict instance is a collection of `BranchPoint => Int` pair. Similarly, the integer is also the node index in the Graph instance. By denfinition, it is also a one-to-one mapping. We can safely reverse this Dict to get another Dict of `Int => BranchPoint`.

A block polymer chain consists of one or more polymer blocks, each of which consists of a number of identical monomers. Polymer blocks are connected by convalent bonds. Each polymer block has two ends. When the end is not connected to other blocks, we call it a free end. Otherwise, we call it a branch point or a joint. We map the architecture of a block copolymer chain to an undirected graph. Free ends and branch points of blocks are nodes of the graph, and blocks themselves are edges of the graph. Therefore, a block is uniquely determined by two node ids.

We use Graphs.jl package to describe the graph. Thus the node id is just an integer. The first added node id is 1.
"""
function build_graph(c::BlockCopolymer)
    e1 = 0 # the node id for the first block end.
    e2 = 0 # the node id for the second block end.
    g = Graph()
    d1 = Dict{PolymerBlock, Tuple{Int, Int}}()
    d2 = Dict{BranchPoint, Int}()
    d3 = Dict{FreeEnd, Int}()
    for b in c.blocks
        if isfreeblockend(b.E1)
            # Each free end is a distinct node in the graph.
            add_vertex!(g)
            e1 = nv(g)
            d3[b.E1] = e1
        else
            # For branch point, we have to make sure if it is already added to the graph.
            # If added, we just use its node id.
            # If not, we have to create a new node corresponds to this branch point.
            if b.E1 ∈ keys(d2)
                e1 = d2[b.E1]
            else
                add_vertex!(g)
                e1 = nv(g)
                d2[b.E1] = e1
            end
        end
        # Do the same thing for the second block end.
        if isfreeblockend(b.E2)
            add_vertex!(g)
            e2 = nv(g)
            d3[b.E2] = e2
        else
            if b.E2 ∈ keys(d2)
                e2 = d2[b.E2]
            else
                add_vertex!(g)
                e2 = nv(g)
                d2[b.E2] = e2
            end
        end
        # A block always corresponds to a new edge in the graph.
        add_edge!(g, e1, e2)
        d1[b] = Polymer._sort_tuple2((e1, e2))
    end
    return g, d1, d2, d3
end

function weighttype(bcg::BlockCopolymerGraph)
    b = first(keys(bcg.block2edge))
    return typeof(b.f)
end

# Implement the AbstractGraph interface.
Base.eltype(bcg::BlockCopolymerGraph) = eltype(bcg.graph)
Graphs.has_edge(bcg::BlockCopolymerGraph, s, d) = has_edge(bcg.graph, s, d)
Graphs.has_vertex(bcg::BlockCopolymerGraph, v) = has_vertex(bcg.graph, v)
Graphs.inneighbors(bcg::BlockCopolymerGraph, v) = inneighbors(bcg.graph, v)
Graphs.ne(bcg::BlockCopolymerGraph) = ne(bcg.graph)
Graphs.nv(bcg::BlockCopolymerGraph) = nv(bcg.graph)
Graphs.outneighbors(bcg::BlockCopolymerGraph, v) = outneighbors(bcg.graph, v)
Graphs.vertices(bcg::BlockCopolymerGraph) = vertices(bcg.graph)
Graphs.is_directed(bcg::BlockCopolymerGraph) = false
Graphs.is_directed(bcg::Type{BlockCopolymerGraph}) = false

function Graphs.edgetype(bcg::BlockCopolymerGraph)
    T = weighttype(bcg)
    return PolymerBlockEdge{T}
end

function Graphs.edges(bcg::BlockCopolymerGraph)
    es = edgetype(bcg)[]
    for e in edges(bcg.graph)
        s, d = src(e), dst(e)
        block = bcg.edge2block[Polymer._sort_tuple2((s, d))]
        push!(es, PolymerBlockEdge(s, d, block))
    end

    return (e for e in es)
end

Graphs.weights(bcg::BlockCopolymerGraph) = bcg.distmx

Base.show(io::IO, bcg::BlockCopolymerGraph) = print(io, "BlockCopolymerGraph as a graph object: $(bcg.graph)")

"""
    chaintype(g::BlockCopolymer)
    chaintype(g::BlockCopolymerGraph)

Return the trait of `PolymerArchitecture`.

We know our chain is always a connected graph. Therefore, it is easy to check whether it has cycles in it. For acyclic connected graph, it is merely a tree. And a tree has exactly (n - 1) edges where n is the number of nodes (vertices).
"""
function chaintype(g::BlockCopolymerGraph)
    # Chain has cycle(s)/ring(s) if and only if its number of edges is different than n - 1 where n is the number of nodes in the graph.
    if ne(g.graph) != nv(g.graph) - 1
        return RingArchitecture()
    end
    ds = degree(g.graph)
    # Linear chain cannot have branch point or can have branch point(s) with degree less than or equal to 2.
    if maximum(ds) <= 2
        return LinearArchitecture()
    end
    # Comb chain should only have branch points of degree 3 in the backbone.
    if sum(ds .== 3) + sum(ds .== 1) == length(ds) && sum(ds .== 3) > 1
        return CombArchitecture()
    end
    # Star chain should have a unique branch point with maximum degree (n>2).
    # Other branch points (if any) should have maximum degree less than n.
    if sum(ds .== maximum(ds)) > 1
        return GeneralBranchedArchitecture()
    else
        return StarArchitecture()
    end
end

chaintype(bcp::BlockCopolymer) = chaintype(BlockCopolymerGraph(bcp))

iscyclicchain(bcp::BlockCopolymer) = iscyclicchain(chaintype(bcp))
isnoncyclicchain(bcp::BlockCopolymer) = !iscyclicchain(bcp)
islinearchain(bcp::BlockCopolymer) = islinearchain(chaintype(bcp))