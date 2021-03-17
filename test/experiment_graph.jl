### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 5b913a92-7f20-11eb-1f10-b3059d07c947
using Polymer

# ╔═╡ 71fa9742-7f20-11eb-369a-ef11d25edb97
using LightGraphs

# ╔═╡ 088b1254-7f21-11eb-0ab9-234397a05e9d
using TikzGraphs, TikzPictures

# ╔═╡ 79fac120-7f2e-11eb-2007-a3bde9fad3f7
using Pkg

# ╔═╡ f20867c8-8017-11eb-2b31-9b24a0235101
using LightGraphs.Experimental

# ╔═╡ a3376b3a-7f2e-11eb-3c1b-b1fdf40461a5
using Revise

# ╔═╡ 110d5194-824b-11eb-19c4-4991b80b89bc


# ╔═╡ 77040be2-7f20-11eb-2132-7d12bbe3751e
function chainAB3A3()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb0 = BranchPoint(:EB0)
	eb1 = BranchPoint(:EB1)
	eb2 = BranchPoint(:EB2)
	eb3 = BranchPoint(:EB3)
	A = PolymerBlock(:A, sA, 0.1, FreeEnd(:A), eb0)
	B1 = PolymerBlock(:B1, sB, 0.2, eb1, eb0)
	B2 = PolymerBlock(:B2, sB, 0.2, eb2, eb0)
	B3 = PolymerBlock(:B3, sB, 0.2, eb3, eb0)
	A1 = PolymerBlock(:A1, sA, 0.1, eb1, FreeEnd(:A1))
	A2 = PolymerBlock(:A2, sA, 0.1, eb2, FreeEnd(:A2))
	A3 = PolymerBlock(:A3, sA, 0.1, eb3, FreeEnd(:A3))
	return BlockCopolymer(:AB3A3, [A, B1, B2, B3, A1, A2, A3])
end

# ╔═╡ 88554b48-7f36-11eb-15b2-6d5cb42186d1
function chainAB3A6()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb0 = BranchPoint(:EB0)
	eb1 = BranchPoint(:EB1)
	eb2 = BranchPoint(:EB2)
	eb3 = BranchPoint(:EB3)
	A = PolymerBlock(:A, sA, 0.4, FreeEnd(:A), eb0)
	A1 = PolymerBlock(:A1, sA, 0.08, eb1, FreeEnd(:A1))
	A2 = PolymerBlock(:A2, sA, 0.08, eb2, FreeEnd(:A2))
	A3 = PolymerBlock(:A3, sA, 0.08, eb3, FreeEnd(:A3))
	A4 = PolymerBlock(:A4, sA, 0.02, eb1, FreeEnd(:A4))
	A5 = PolymerBlock(:A5, sA, 0.02, eb2, FreeEnd(:A5))
	A6 = PolymerBlock(:A6, sA, 0.02, eb3, FreeEnd(:A6))
	B1 = PolymerBlock(:B1, sB, 0.12, eb1, eb0)
	B2 = PolymerBlock(:B2, sB, 0.12, eb2, eb0)
	B3 = PolymerBlock(:B3, sB, 0.06, eb3, eb0)
	return BlockCopolymer(:AB3A3, [A, B1, B2, B3, A1, A2, A3, A4, A5, A6])
end

# ╔═╡ 0e1c3dee-7fc9-11eb-0b44-ad794b13a6ca
Symbol("EB"*string(2))

# ╔═╡ f0e9ac34-7fc8-11eb-27a3-b95c83218a18
function branchpoints(n; prefix="EB")
	bps = BranchPoint[]
	return [BranchPoint(Symbol(prefix*string(i))) for i in 1:n]
end

# ╔═╡ a166207e-7fc9-11eb-2ada-77f2740283b4
function freeends(n; prefix="A")
	fes = FreeEnd[]
	return [FreeEnd(Symbol(prefix*string(i))) for i in 1:n]
end

# ╔═╡ 36cfce44-824a-11eb-367e-9b18ea4df6bf
function chainABA()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb = branchpoints(2)
	fe = freeends(2)
	A1 = PolymerBlock(:A1, sA, 0.4, eb[1], fe[1])
	A2 = PolymerBlock(:A2, sA, 0.4, eb[2], fe[2])
	B = PolymerBlock(:B, sB, 0.2, eb[1], eb[2])
	return BlockCopolymer(:ABA, [A1, A2, B])
end

# ╔═╡ 997c616a-824a-11eb-0a7a-e7dcb8a66df0
chainABAg = BlockCopolymerGraph(chainABA())

# ╔═╡ 8c2f0162-824e-11eb-3a4d-7110b677cbea
function starAB4A8()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb = branchpoints(5)
	fe = freeends(8)
	A = PolymerBlock(:A, sA, 0.32, FreeEnd(:A), eb[1])
	A1 = PolymerBlock(:A1, sA, 0.04, eb[2], fe[1])
	A2 = PolymerBlock(:A2, sA, 0.08, eb[3], fe[2])
	A3 = PolymerBlock(:A3, sA, 0.04, eb[4], fe[3])
	A4 = PolymerBlock(:A4, sA, 0.08, eb[5], fe[4])
	A5 = PolymerBlock(:A5, sA, 0.06, eb[2], fe[5])
	A6 = PolymerBlock(:A6, sA, 0.02, eb[3], fe[6])
	A7 = PolymerBlock(:A7, sA, 0.06, eb[4], fe[7])
	A8 = PolymerBlock(:A8, sA, 0.02, eb[5], fe[8])
	B1 = PolymerBlock(:B1, sB, 0.1, eb[2], eb[1])
	B2 = PolymerBlock(:B2, sB, 0.04, eb[3], eb[1])
	B3 = PolymerBlock(:B3, sB, 0.1, eb[4], eb[1])
	B4 = PolymerBlock(:B3, sB, 0.04, eb[5], eb[1])
	return BlockCopolymer(:AB4A8, [A, B1, B2, B3, B4, A1, A2, A3, A4, A5, A6, A7, A8])
end

# ╔═╡ 7df1a45a-824f-11eb-1538-23b1c5782c75
starAB4A8()

# ╔═╡ 993a0e20-824f-11eb-2963-fd2dcbf63882
starg = BlockCopolymerGraph(starAB4A8());

# ╔═╡ 748a414a-7fc9-11eb-0250-6312741527d7
branchpoints(8)

# ╔═╡ f1b33da2-7fc9-11eb-1f75-bf245a82d0ea
freeends(6; prefix="B")

# ╔═╡ 469c1d94-7fb8-11eb-307b-a13577067b52
function branchAB()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb = branchpoints(8)
	fe = freeends(6; prefix="B")
	A1 = PolymerBlock(:A1, sA, 0.12, eb[1], eb[3])
	A2 = PolymerBlock(:A2, sA, 0.12, eb[2], eb[3])
	A3 = PolymerBlock(:A1, sA, 0.12, eb[6], eb[8])
	A4 = PolymerBlock(:A1, sA, 0.12, eb[6], eb[7])
	A5 = PolymerBlock(:A1, sA, 0.2, eb[4], eb[5])
	B1 = PolymerBlock(:B1, sB, 0.02, eb[1], fe[1])
	B2 = PolymerBlock(:B2, sB, 0.02, eb[1], fe[2])
	B3 = PolymerBlock(:B3, sB, 0.02, eb[2], fe[3])
	B4 = PolymerBlock(:B4, sB, 0.02, eb[8], fe[4])
	B5 = PolymerBlock(:B5, sB, 0.02, eb[7], fe[5])
	B6 = PolymerBlock(:B6, sB, 0.02, eb[7], fe[6])
	B7 = PolymerBlock(:B7, sB, 0.08, eb[3], eb[4])
	B8 = PolymerBlock(:B8, sB, 0.12, eb[5], eb[6])
	return BlockCopolymer(:AB3A3, [A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, B6, B7, B8])
end

# ╔═╡ 1dd95408-8247-11eb-2a29-09164b734f9b
function branchAB2()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb = branchpoints(11)
	fe = freeends(6; prefix="B")
	fe10 = FreeEnd(:B10)
	fe11 = FreeEnd(:B11)
	fe12 = FreeEnd(:B12)
	A1 = PolymerBlock(:A1, sA, 0.08, eb[1], eb[3])
	A2 = PolymerBlock(:A2, sA, 0.08, eb[2], eb[3])
	A3 = PolymerBlock(:A3, sA, 0.08, eb[6], eb[8])
	A4 = PolymerBlock(:A4, sA, 0.08, eb[6], eb[7])
	A5 = PolymerBlock(:A5, sA, 0.04, eb[4], eb[5])
	A6 = PolymerBlock(:A6, sA, 0.08, eb[9], eb[10])
	A7 = PolymerBlock(:A7, sA, 0.08, eb[9], eb[11])
	B1 = PolymerBlock(:B1, sB, 0.02, eb[1], fe[1])
	B2 = PolymerBlock(:B2, sB, 0.02, eb[1], fe[2])
	B3 = PolymerBlock(:B3, sB, 0.02, eb[2], fe[3])
	B4 = PolymerBlock(:B4, sB, 0.02, eb[8], fe[4])
	B5 = PolymerBlock(:B5, sB, 0.02, eb[7], fe[5])
	B6 = PolymerBlock(:B6, sB, 0.02, eb[7], fe[6])
	B7 = PolymerBlock(:B7, sB, 0.09, eb[3], eb[4])
	B8 = PolymerBlock(:B8, sB, 0.12, eb[5], eb[6])
	B9 = PolymerBlock(:B9, sB, 0.09, eb[5], eb[9])
	B10 = PolymerBlock(:B10, sB, 0.02, eb[10], fe10)
	B11 = PolymerBlock(:B11, sB, 0.02, eb[11], fe11)
	B12 = PolymerBlock(:B12, sB, 0.02, eb[11], fe12)
	return BlockCopolymer(:AB3A3, [A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, B8, B9, B10, B11, B12])
end

# ╔═╡ 53e6ef00-8248-11eb-1837-ad033d85c8f8
branchAB2()

# ╔═╡ 277a4454-7fb8-11eb-08f0-c308ad0fcd32
branchAB()

# ╔═╡ 8d194a3c-7f20-11eb-105a-0b91f6a69bc1
# chain = chainAB3A3()
chain = chainAB3A6()

# ╔═╡ 6d608ed6-7f49-11eb-25ef-272fdae43007
chaintype(chain)

# ╔═╡ cd4c77a2-824d-11eb-2f65-998feeb8cf5e
chaintype(branchAB())

# ╔═╡ dcc4cf54-824d-11eb-04ed-138f3813ab5e
chaintype(branchAB2())

# ╔═╡ ec653e1a-824d-11eb-290c-6dead284203d
chaintype(chainABA())

# ╔═╡ 9567dd66-7f20-11eb-0bda-315ad04a42ff
bcg = BlockCopolymerGraph(chain);

# ╔═╡ e31ec7ec-824d-11eb-0101-714661027e0a
chaintype(bcg)

# ╔═╡ a8d5721e-7f20-11eb-1c3f-fd4ffea1637e
source = argmax(degree(bcg.graph))

# ╔═╡ ba365686-7f20-11eb-22c6-0782aaa5897d
bfs_tree(bcg.graph, source) |> plot

# ╔═╡ f506bc06-800b-11eb-197a-fb31c6d4dd22
dfs_tree(bcg.graph, source) |> plot

# ╔═╡ 0278037c-800c-11eb-10c9-21130c78671c
(dfs1 = dfs_tree(bcg.graph, 1)) |> plot

# ╔═╡ 1bac3db8-800c-11eb-27ad-d10ce333659f
(dfs6 = dfs_tree(bcg.graph, 6)) |> plot

# ╔═╡ 28d4a0e6-800c-11eb-3010-492eada4bce8
(dfs7 = dfs_tree(bcg.graph, 7)) |> plot

# ╔═╡ 1989d0c2-7f21-11eb-13a6-2f87902db271
subg, nvmap = induced_subgraph(bcg.graph, [2,9,6,3])

# ╔═╡ 84157cf4-7f29-11eb-36cf-59dea7d7f5b6
nvmap

# ╔═╡ c3493360-7f20-11eb-0a22-637ae33e5833
plot(subg)

# ╔═╡ 5e0481d4-7f21-11eb-1b31-2f6dd8b2738c
subg2, nvmap2 = induced_subgraph(bcg.graph, [4,7])

# ╔═╡ 15274aa6-7f4d-11eb-09c8-5dda2485fe85
subg2 |> plot

# ╔═╡ e08c53ce-7f51-11eb-1e85-1386c9a2757d
"""
Find those neighbors and all their neighborhoods that are connected to vertex v, excluding the one vfrom and its all neighborhoods.
"""
function find_neighbors!(v, vfrom, vertices, graph)
	for vx in neighbors(graph, v)
		(vx == vfrom) && continue
		push!(vertices, vx)
		if degree(graph, vx) == 1
			continue
		else
			find_neighbors!(vx, v, vertices, graph)
		end
	end
	return nothing
end

# ╔═╡ c04a2b14-8002-11eb-38ce-218064340191
function find_neighbors!(v, vfrom, vertices, bcg::BlockCopolymerGraph)
	return find_neighbors!(v, vfrom, vertices, bcg.graph)
end

# ╔═╡ 5df84cf8-7f3e-11eb-2840-b98a8db5d97c
"""
Make two induced subgraphs of the input graph by cutting the edge e.
Note that the input edge is not included in the induced subgraphs.
"""
function induced_subtree(e::Edge, graph)
	v1, v2 = e.src, e.dst
	vertices1 = [v1]
	find_neighbors!(v1, v2, vertices1, graph)
	vertices2 = [v2]
	find_neighbors!(v2, v1, vertices2, graph)
	return induced_subgraph(graph, vertices1), induced_subgraph(graph, vertices2)
end

# ╔═╡ fed4c060-8002-11eb-3e86-e55799b0498d
function induced_subtree(e::Edge, bcg::BlockCopolymerGraph)
	return induced_subtree(e, bcg.graph)
end

# ╔═╡ ad015e6c-7fba-11eb-2a22-63b0cab369d3
"""
Make induced subgraphs of the input graph by disconnecting all branches from vertex v.
Note that vertex v is included in all induced subgraphs.
"""
function induced_subtree(v, graph)
	subgraphs = []
	for vx in neighbors(graph, v)
		vertices = [v, vx]
		find_neighbors!(vx, v, vertices, graph)
		push!(subgraphs, induced_subgraph(graph, vertices))
	end
	return subgraphs
end

# ╔═╡ 372af4ac-8003-11eb-0b48-b5372458c90b
function induced_subtree(v, bcg::BlockCopolymerGraph)
	return induced_subtree(v, bcg.graph)
end

# ╔═╡ 9791b2aa-7ff3-11eb-26b0-1ba7eaba36ed
function reverse_vmap(vmap)
	d = Dict()
	for i in 1:length(vmap)
		d[vmap[i]] = i
	end
	return d
end

# ╔═╡ 7fbcff20-7ff5-11eb-37a1-651728486d23
reverse_vmap(nvmap)

# ╔═╡ 98e92396-7ff2-11eb-3782-3f5616a2434c
"""
Cut the input graph at edges e1 and e2.
Return the subtree wich contains both edges' vertices, a vmap which maps the vertex of the subtree to the original graph, and v1 of e1 and v2 of e2 both present in the returned subtree.
"""
function induced_subtree(e1::Edge, e2::Edge, graph)
	# Make sure e1 and e2 has are two distinct edges
	length(unique([e1.src, e1.dst, e2.src, e2.dst])) < 3 && return nothing
	
	# Make two subgraphs by cutting edge e1
	subtrees1 = induced_subtree(e1, graph)
	subtree1_src, vmap1_src = subtrees1[1]
	subtree1_dst, vmap1_dst = subtrees1[2]

	# Find out edge e2 is in which subtree1 (src or dst).
	# index of vmap1 is vertex of subtree1
	# vlaue of vmap1 is vertex of graph
	issrc = e2.src ∈ vmap1_src
	subtree1 = issrc ? subtree1_src : subtree1_dst
	vmap1 = issrc ? vmap1_src : vmap1_dst
	v1 = issrc ? e1.src : e1.dst
	
	# Make two subgraphs of subtree by cutting edge e2
	# Note that we have to map e2 of bcg to subtree1 edge.
	# key of rvmap1 is vertex of graph
	# value of rvmap1 is vertex of subtree1
	rvmap1 = reverse_vmap(vmap1)
	subtrees2 = induced_subtree(Edge(rvmap1[e2.src], rvmap1[e2.dst]), subtree1)
	# index of vmap2_* is vertex of subgraph of subtree1
	# value of vmap2_* is vertex of subtree1
	subtree2_src, vmap2_src = subtrees2[1]
	subtree2_dst, vmap2_dst = subtrees2[2]
	
	# Find out v1 from edge e1 is in which subtree2 (src or dst)
	issrc = rvmap1[v1] ∈ vmap2_src
	subtree2 = issrc ? subtree2_src : subtree2_dst
	vmap2 = issrc ? vmap2_src : vmap2_dst
	v2 = issrc ? e2.src : e2.dst
	
	# Remaping from vertex of subtree1 to orginal graph
	vmap = similar(vmap2)
	for i in 1:length(vmap2)
		vmap[i] = vmap1[vmap2[i]]
	end
	
	return subtree2, vmap, v1, v2
end

# ╔═╡ 60e61312-8003-11eb-3084-f357cd092257
function induced_subtree(e1::Edge, e2::Edge, bcg::BlockCopolymerGraph)
	return induced_subtree(e1, e2, bcg.graph)
end

# ╔═╡ 4b1029e2-801c-11eb-04f6-135d3cf91f6b
function induced_subtree(e1::AbstractVector, e2::AbstractVector, bcg::BlockCopolymerGraph)
	return induced_subtree(Edge(e1[1],e1[2]), Edge(e2[1],e2[2]), bcg)
end

# ╔═╡ b49fadaa-801d-11eb-01af-1be437348d56
function induced_subtree(e1::Tuple{Int,Int}, e2::Tuple{Int,Int}, bcg::BlockCopolymerGraph)
	return induced_subtree(Edge(e1[1],e1[2]), Edge(e2[1],e2[2]), bcg)
end

# ╔═╡ 9ebd2494-801e-11eb-276f-93ae1af1b013
function find_edge_in_subtree(v, vmap, bcg::BlockCopolymerGraph)
	for vx in neighbors(bcg.graph, v)
		(vx ∈ vmap) && return Edge(v, vx)
	end
end

# ╔═╡ ca86285c-801c-11eb-2a1b-0ff444ed0721
"""
An induced subgraph of the input graph that contains vertices v1 and v2 but excluded all vertices in vmap1 and vmap2.
It is equivalent to say the operation split the original graph into three parts: an induced subgraph with vertices in vmap1, an induced subgraph with vertices in vmap2, and the returned induced subgraph.
"""
function induced_subtree(v1::Int, v2::Int, vmap1, vmap2, bcg::BlockCopolymerGraph)
	e1 = find_edge_in_subtree(v1, vmap1, bcg)
	e2 = find_edge_in_subtree(v2, vmap2, bcg)
	subtree, vmap, v1out, v2out = induced_subtree(e1, e2, bcg)
	# push!(vmap, v1)
	# push!(vmap, v2)
	return induced_subgraph(bcg.graph, vmap)
end

# ╔═╡ f51febe4-8598-11eb-3f19-7b2888b0723a
induced_subtree(2, bcg)

# ╔═╡ 8344b206-8085-11eb-39c7-113ca193e89d
function is_symmetric_tree(tree, v1, v2, bcg::BlockCopolymerGraph; vmap=collect(1:nv(bcg.graph)))
	return is_symmetric_tree(tree, v1, v2, vmap, bcg)
end

# ╔═╡ d4f01daa-8087-11eb-0dc2-974351888dd7
function is_symmetric_tree(v1, v2, bcg::BlockCopolymerGraph)
	vmap = collect(1:nv(bcg.graph))
	return is_symmetric_tree(bcg.graph, v1, v2, vmap, bcg)
end

# ╔═╡ a2c40156-7f52-11eb-27f2-cfee4296f1b2
subtree1, subtree2 = induced_subtree(Edge(2, 3), bcg)

# ╔═╡ edd548a2-7fbc-11eb-0956-cfb97b1c695d
subtree1 |> first |> plot

# ╔═╡ f0fce2f6-7fb7-11eb-1ece-b1e2494058db
subtree2 |> first |> plot

# ╔═╡ aaeba326-7fba-11eb-116e-296d4a438f04
subtrees = induced_subtree(argmax(degree(bcg.graph)), bcg)

# ╔═╡ 784a7912-7fd3-11eb-2591-17d4b1eaa347
branch = branchAB();

# ╔═╡ 8aa193c0-7fd3-11eb-156b-c12d122a18c3
branchg = BlockCopolymerGraph(branch);

# ╔═╡ f336d56a-80a0-11eb-2056-6d2cdf79d46a
induced_subtree(Edge(2,3), Edge(4,5), branchg)

# ╔═╡ b2984fbe-7ff3-11eb-0e1c-dbd158dd04ae
induced_subtree(Edge(2,7), Edge(4,8), branchg)

# ╔═╡ 3cd857c2-801d-11eb-3591-dd51a1c57699
induced_subtree([2,1], [4,5], branchg)

# ╔═╡ c63d92fc-801d-11eb-1a6b-fdd23bde076a
induced_subtree((2,1), (4,5), branchg)

# ╔═╡ d060a24a-801f-11eb-06ae-73cb7c729035
begin
	st1, vm1 = induced_subtree(Edge(2,7), branchg)[1]
	st2, vm2 = induced_subtree(Edge(4,8), branchg)[1]
	st3, vm3 = induced_subtree(2, 4, vm1, vm2, branchg)
	@show vm3
	plot(st3)
end

# ╔═╡ 555b5956-8082-11eb-0b1e-c9eef34c041e
(dfs_st3_3 = dfs_tree(st3, 3)) |> plot

# ╔═╡ a5898600-8082-11eb-24dc-df23129e37c4
(dfs_st3_4 = dfs_tree(st3, 4)) |> plot

# ╔═╡ b878ddaa-7fd4-11eb-06ef-3bac1ac83c38
branchg.node2joint

# ╔═╡ 83c1e328-7fe2-11eb-38cf-116d4d032ce9
branchg.node2free

# ╔═╡ 65aece8a-7fd3-11eb-3978-27fa9683e21b
subbranchtrees = induced_subtree(7, branchg)

# ╔═╡ 0e7d8880-7fe2-11eb-1f77-754c4263a0b5
subbranchtree1, _b1 = subbranchtrees[1]

# ╔═╡ 22a8a6b2-7fe2-11eb-2037-55a1a1e060f9
subbranchtree2, _b2 = subbranchtrees[2]

# ╔═╡ 83eba7e2-7f21-11eb-2685-d781ba1b10a3
difference(subg, subg2) |> plot

# ╔═╡ 97517fb4-7f21-11eb-0547-0963ab9ca934
star_graph(9) |> plot

# ╔═╡ 814c4f0c-7f22-11eb-13af-fb3eb53810b3
path_graph(4) |> plot

# ╔═╡ 053fe0b2-7f4d-11eb-09b8-41ba9a4d6387
blockdiag(subg, subg2) |> plot

# ╔═╡ 2c395a22-7f4d-11eb-2435-c9660719b17a
join(subg, subg2) |> plot

# ╔═╡ 319a9c80-7f2d-11eb-1a0c-bbc78da2540a
(3,2) |> Polymer._sort_tuple2

# ╔═╡ 016996e4-7f2e-11eb-3301-6da3c7b996ce
bcg.edge2block

# ╔═╡ 0161b5be-801b-11eb-2317-c539091f1e8d
function equivalent_blockend(v1, v2, bcg::BlockCopolymerGraph; vmap1=collect(1:nv(bcg.graph)), vmap2=collect(1:nv(bcg.graph)))
	# convert vertices to BlockCopolymerGraph vertices
	v1n = vmap1[v1]
	v2n = vmap2[v2]
	# check the degree
	(degree(bcg.graph, v1n) == degree(bcg.graph, v2n)) && return true
	# Otherwise
	return false
end

# ╔═╡ a9bd1960-7f38-11eb-1c7f-cd3105f49480
"""
v1 is in the input bcg.graph.
v2 is in the induced subgraph of the input bcg.graph.
"""
function equivalent_blockend(v1, v2, vmap, bcg::BlockCopolymerGraph)
	return equivalent_blockend(v1, v2, bcg; vmap2=vmap)
end

# ╔═╡ ae82e0d8-7fe4-11eb-3655-41da2101eb5a
function equivalent_blockend(v1, v2, vmap1, vmap2, bcg)
	return equivalent_blockend(v1, v2, bcg; vmap1=vmap1, vmap2=vmap2)
end

# ╔═╡ 89a27b48-801b-11eb-0a17-414e085d11d0
equivalent_blockend(4, 7, _b1, _b2, branchg)

# ╔═╡ cb84ef14-801b-11eb-0a4e-e5268c3a9965
# 4 => 9
# 7 => 13
equivalent_blockend(4, 7, branchg; vmap1=_b1, vmap2=_b2)

# ╔═╡ 33c8dfa8-8018-11eb-26a8-4b7c0122ee31
"""
Edges e1 and e2 can be in the induced subgraphs of input bcg::BlockCopolymerGraph.
The maps from subgraphs to original graph are provided by vmap1 and vmap2 for e1 and e2, respectively.
By default, vmap1 and vmap2 are identical mapping, which means e1 and e2 are both in the input bcg.graph.
"""
function equivalent_block(e1::Edge, e2::Edge, bcg::BlockCopolymerGraph; vmap1=collect(1:nv(bcg.graph)), vmap2=collect(1:nv(bcg.graph)))
	# convert edges to BlockCopolymerGraph edges
	v11, v12 = e1.src, e1.dst
	e1n = (vmap1[v11], vmap1[v12]) |> Polymer._sort_tuple2
	# convert to BlockCopolymerGraph
	v21, v22 = e2.src, e2.dst
	e2n = (vmap2[v21], vmap2[v22]) |> Polymer._sort_tuple2
	# compare if this two edges are blocks with identical properties
	block1 = bcg.edge2block[e1n]
	block2 = bcg.edge2block[e2n]
	# 1. Edges correspond to same block
	(block1 == block2) && return true
	# 2. Edges correspond to two blocks but with same specie and length
	(specie(block1) == specie(block2) && block1.f == block2.f) && return true
	# Otherwise
	return false
end

# ╔═╡ dd1e0226-8019-11eb-14e7-8767d9361a82
function equivalent_block(e1, e2, bcg::BlockCopolymerGraph; vmap1=collect(1:nv(bcg.graph)), vmap2=collect(1:nv(bcg.graph)))
	return equivalent_block(Edge(e1[1], e1[2]), Edge(e2[1], e2[2]), bcg; vmap1=vmap1, vmap2=vmap2)
end

# ╔═╡ 5025a718-7f2c-11eb-2779-a5bad86fd6a0
"""
e1 is in the input bcg.graph.
e2 is in the induced subgraph of the input bcg.graph.
"""
function equivalent_block(e1, e2, vmap, bcg::BlockCopolymerGraph)
	return equivalent_block(e1, e2, bcg; vmap2=vmap)
end

# ╔═╡ dc1bb0d8-7fe4-11eb-1778-5dd086ed216a
function equivalent_block(e1, e2, vmap1, vmap2, bcg::BlockCopolymerGraph)
	return equivalent_block(e1, e2, bcg; vmap1=vmap1, vmap2=vmap2)
end

# ╔═╡ 46ebdeb6-800c-11eb-203f-bb21e582625e
has_isomorph(dfs6, dfs7; 
			vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,bcg),
			edge_relation=(e1,e2)->equivalent_block(e1,e2,bcg))

# ╔═╡ 2b3fc368-8083-11eb-3504-9b5892365b45
"""
We say a tree is symmetric about two vertices v1 and v2, if its depth first trasversal path starting from v1 and starting from v2 are isomorphic, in the sense that all blocks are equivalent, i.e. the specie and length of the block are identical.

* tree can be bcg.graph or any induced subgraph of bcg.graph.
* v1, v2 are the vertices of tree.
* vmap is a map which maps the vertices of tree to vertices of bcg.graph.
* bcg is a BlockCopolymerGraph corresponding to a BlockCopolymer object.
"""
function is_symmetric_tree(tree, v1, v2, vmap, bcg::BlockCopolymerGraph)
	# If tree only has two vertices (one edge), then it is symmetric.
	(nv(tree) == 2) && return true
	
	dfs1 = dfs_tree(tree, v1)
	dfs2 = dfs_tree(tree, v2)
	return has_isomorph(dfs1, dfs2;
				vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,vmap,vmap,bcg),
				edge_relation=(e1,e2)->equivalent_block(e1,e2,vmap,vmap,bcg))
end

# ╔═╡ 711caa16-8085-11eb-1a2c-7133fd44aa62
is_symmetric_tree(bcg.graph, 6, 7, bcg)

# ╔═╡ 072042e6-8088-11eb-2c38-bd19e1578001
is_symmetric_tree(6, 7, bcg)

# ╔═╡ c18f4e4a-8085-11eb-013b-97978173c7cc
is_symmetric_tree(bcg.graph, 1, 6, bcg)

# ╔═╡ 0d6fec00-8088-11eb-36d5-719ab414513e
is_symmetric_tree(1, 6, bcg)

# ╔═╡ aba16426-8083-11eb-055a-6500d0de28de
is_symmetric_tree(st3, 3, 4, vm3, branchg)

# ╔═╡ d99eb46a-8082-11eb-2bf0-cdbe6f0aacc9
has_isomorph(dfs_st3_3, dfs_st3_4; 
			vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,vm3,vm3,branchg),
			edge_relation=(e1,e2)->equivalent_block(e1,e2,vm3,vm3,branchg))

# ╔═╡ a738558a-8019-11eb-35a1-991259834d46
# (1, 2) => (7, 2)
# (2, 3) => (8, 4)
equivalent_block((1,2), (2,3), branchg; vmap1=_b1, vmap2=_b2)

# ╔═╡ 1db19018-801b-11eb-152f-01608867a266
equivalent_block((1,2), (2,3), _b1, _b2, branchg)

# ╔═╡ 54d6310a-8018-11eb-0468-2b1e98c5a6a1
nv(bcg.graph)

# ╔═╡ d783467e-7f2d-11eb-3266-498d9ed39e5f
equivalent_block(Edge(4,7), Edge(3,6), bcg)

# ╔═╡ 86bde198-801a-11eb-3412-cbc6fb494aec
equivalent_block((4,7), (3,6), bcg)

# ╔═╡ 90a74d98-801a-11eb-2ecd-19f4b2ed9e72
equivalent_block([4,7], [3,6], bcg)

# ╔═╡ b8ea5c08-7f22-11eb-002d-e923402835ef
LightGraphs.Experimental.all_induced_subgraphisomorph(bcg.graph, subg;
			vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,nvmap,bcg),
			edge_relation=(e1,e2)->equivalent_block(e1,e2,nvmap,bcg)) |> collect

# ╔═╡ f66615e4-7f28-11eb-3bd7-8f123f3b8837
LightGraphs.Experimental.count_induced_subgraphisomorph(bcg.graph, subg)

# ╔═╡ 4b24db38-7f29-11eb-0164-7fa310866165
LightGraphs.Experimental.all_induced_subgraphisomorph(branchg.graph, subbranchtree1;
			vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,_b1,branchg),
			edge_relation=(e1,e2)->equivalent_block(e1,e2,_b1,branchg)) |> collect

# ╔═╡ 61362204-7f2b-11eb-00cc-dfa70199b912
has_isomorph(branchg.graph, branchg.graph)

# ╔═╡ f3d9ee22-80ba-11eb-3ecc-a54d00a89045
function find_leaf_neighbor(v, bcg::BlockCopolymerGraph)
	return first(neighbors(bcg.graph, v))
end

# ╔═╡ cbaca3e4-80be-11eb-3938-7133e143f92e
function count_isomorphic_subtree(subtree::AbstractGraph, vmap, bcg::BlockCopolymerGraph)
	return count_induced_subgraphisomorph(bcg.graph, subtree;
			vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,vmap,bcg), 
			edge_relation=(e1,e2)->equivalent_block(e1,e2,vmap,bcg))
end

# ╔═╡ c0de8e22-80bf-11eb-039e-054f2df66301
"""
Find all isomorphic subtrees for input subtree in bcg.graph.
v: the front vertex of the input subtree (its index is in the bcg.graph).
vmap: mapping vertex in subtree to bcg.graph. key (subtree) => value (bcg.graph).
"""
function all_isomorphic_subtree(subtree::AbstractGraph, v, vmap, bcg::BlockCopolymerGraph)
	ilist = all_induced_subgraphisomorph(bcg.graph, subtree;
			vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,vmap,bcg), 
			edge_relation=(e1,e2)->equivalent_block(e1,e2,vmap,bcg)) |> collect
	subtrees = []
	front_vertex = Int[]
	for subtree in ilist
		vertices = Int[]
		for (v1, v2) in subtree
			push!(vertices, v1)
			(vmap[v2] == v) && push!(front_vertex, v1)
		end
		push!(subtrees, vertices)
	end
	return subtrees, front_vertex
end

# ╔═╡ 1845fcee-823d-11eb-2550-015a785a6453
all_isomorphic_subtree(induced_subgraph(bcg.graph, [6,3])[1], 3, [6,3], bcg)

# ╔═╡ 71ed8db8-823d-11eb-1486-c32cd5de0c38
induced_subgraph(bcg.graph, [6,3])

# ╔═╡ e012cbd2-8191-11eb-3ff1-69c1364ad19c
function is_equivalent_subtree(subtree1, subtree2, v1, v2, bcg::BlockCopolymerGraph)
	subtree, vmap = induced_subtree(v1, v2, subtree1, subtree2, bcg)
	rvmap = reverse_vmap(vmap)
	return is_symmetric_tree(subtree, rvmap[v1], rvmap[v2], vmap, bcg)
end

# ╔═╡ 7c4d201e-818b-11eb-2d02-218a6ed92a49
function process_isolated_isomorphic_subtrees(subtrees, front_vertices, bcg::BlockCopolymerGraph)
	equivalent_subtrees = []
	equiv_front_vertices = []
	semi_equivalent_subtrees = []
	semi_equiv_front_vertices = []
	
	m = length(subtrees)
	
	if m == 1
		# Here subtrees = [vertices] and front_vertices = [v]
		push!(semi_equivalent_subtress, subtrees)
		push!(semiequiv_front_vertices, front_vertices)
	end
	
	if m == 2
		vertices1, vertices2 = subtrees
		v1, v2 = front_vertices
		issym = is_equivalent_subtree(vertices1, vertices2, v1, v2, bcg)
		if issym
			push!(equivalent_subtrees, [vertices1, vertices2])
			push!(equiv_front_vertices, [v1, v2])
		else
			push!(semi_equivalent_subtrees, [vertices1, vertices2])
			push!(semi_equiv_front_vertices, [v1, v2])
		end
	end
	
	processed = [false for _ in 1:m]
	if m > 2
		for i in 1:m
			equiv_subtrees_subset = []
			equiv_vertex_subset = []
			semiequiv_subtrees_subset = []
			semiequiv_vertex_subset = []
			processed[i] && continue
			vertices1 = subtrees[i]
			v1 = front_vertices[i]
			for j in (i+1):m
				processed[j] && continue
				vertices2 = subtrees[j]
				v2 = front_vertices[j]
				issym = is_equivalent_subtree(vertices1, vertices2, v1, v2, bcg)
				if issym
					push!(equiv_subtrees_subset, vertices2)
					push!(equiv_vertex_subset, v2)
					processed[j] = true
				end
			end
			if length(equiv_subtrees_subset) > 0
				push!(equiv_subtrees_subset, vertices1)
				push!(equiv_vertex_subset, v1)
				push!(equivalent_subtrees, equiv_subtrees_subset)
				push!(equiv_front_vertices, equiv_vertex_subset)
				processed[i] = true
			end
		end
		for i in 1:m
			if !processed[i]
				push!(semi_equivalent_subtrees, subtrees[i])
				push!(semi_equiv_front_vertices, front_vertices[i])
				processed[i] = true
			end
		end
	end
	
	return equivalent_subtrees, equiv_front_vertices, semi_equivalent_subtrees, semi_equiv_front_vertices
end

# ╔═╡ 3c62ebe6-8198-11eb-3255-dfbaa392145f
"""
For each equivalent subtree, there should be only one branch (emanating from the front vertex of that subtree) contains all other equivalent subtrees.
Therefore, we can always propagate one edge.

The special case is when there are only two equivalent subtrees, and their front vertices form an edge. In such case, we can still propagate one edge, and the front vertices of resulted two equivalent subtrees exchange their positions.
"""
function process_equivalent_subtrees(subtrees, front_vertices, bcg::BlockCopolymerGraph)
	equivalent_subtrees = []
	equivalent_front_vertices = []
	isomorphic_subtrees = []
	isomorphic_front_vertices = []
	
	m = length(subtrees)
	
	(m == 1) && error("At least two equivalent subtrees must be provided. Only got one.")
	
	# Finalize equivalent trees when there are only two equivalent trees and they are connected by an edge. No more isomorphic subtrees are returned.
	if m == 2
		subtree1, subtree2 = subtrees
		v1, v2 = front_vertices
		if has_edge(bcg.graph, v1, v2)
			# Note we have to include all other branches emanating from v1 and v2 into our two equivalent trees.
			# In this case, the equivalent trees are final, meaning that they can not be expanded further.
			(g1, vmap1), (g2, vmap2) = induced_subtree(Edge(v1, v2), bcg)
			push!(vmap1, v2)
			push!(vmap2, v1)
			push!(equivalent_subtrees, vmap1)
			push!(equivalent_subtrees, vmap2)
			push!(equivalent_front_vertices, v2)
			push!(equivalent_front_vertices, v1)
			return equivalent_subtrees, equivalent_front_vertices, isomorphic_subtrees, isomorphic_front_vertices
		end
	end
	
	for i in 1:m
		subtree = subtrees[i]
		v = front_vertices[i]
		isomorphic_vertices = copy(subtree)
		new_front_vertex = nothing
		for vx in neighbors(bcg.graph, v)
			(vx ∈ subtree) && continue
			(g, vmap), _ = induced_subtree(Edge(vx, v), bcg)
			tomerge = true
			for vi in front_vertices
				# we are in the branch containing all other equivalent subtrees
				(vi ∈ vmap) && (tomerge = false; break) 
			end
			tomerge ? append!(isomorphic_vertices, vmap) : new_front_vertex = vx
		end
		push!(isomorphic_vertices, new_front_vertex)
		push!(isomorphic_subtrees, isomorphic_vertices)
		push!(isomorphic_front_vertices, new_front_vertex)
	end

	return equivalent_subtrees, equivalent_front_vertices, isomorphic_subtrees, isomorphic_front_vertices
end

# ╔═╡ a719a518-81ae-11eb-3935-a3ef2903b1cf
function group_isomorphic_subtrees(subtrees, front_vertices, bcg::BlockCopolymerGraph)
	unique_vertices = front_vertices |> unique
	subtree_groups = []
	vertex_groups = []
	ntrees = []
	for v in unique_vertices
		subtree_group = []
		for i in 1:length(front_vertices)
			subtree = subtrees[i]
			vertex = front_vertices[i]
			(vertex == v) && push!(subtree_group, subtree)
		end
		push!(subtree_groups, subtree_group)
		push!(vertex_groups, v)
		push!(ntrees, length(subtree_group))
	end
	
	return subtree_groups, vertex_groups, ntrees
end

# ╔═╡ 7668845c-8217-11eb-1c04-a1cb72e69815
function process_isomorphic_subtree_groups(groups, front_vertices, ntrees)
	equivalent_subtrees = []
	equivalent_front_vertices = []
	semi_equivalent_subtrees = []
	semi_equivalent_front_vertices = []
	isolated_isomorphic_subtrees = []
	isolated_isomorphic_vertices = []
	
	sum(ntrees) < 2 && error("At least two isomorphic subtrees in the groups.")
	
	unique_ntrees = ntrees |> unique
	for n in unique_ntrees
		Qgroup = []
		Qv = []
		for i in 1:length(groups)
			if ntrees[i] == n
				push!(Qgroup, groups[i])
				push!(Qv, front_vertices[i])
			end
		end
		if length(Qgroup) == 1
			if n == 1
				push!(semi_equivalent_subtrees, Qgroup)
				push!(semi_equivalent_front_vertices, Qv)
			else
				push!(equivalent_subtrees, Qgroup)
				push!(equivalent_front_vertices, Qv)
			end
		else
			Qgroup_new = []
			for j in 1:length(Qgroup)
				new_subtree = []
				for subtree in Qgroup[j]
					append!(new_subtree, subtree)
				end
				push!(Qgroup_new, unique(new_subtree))
			end
			push!(isolated_isomorphic_subtrees, Qgroup_new)
			push!(isolated_isomorphic_vertices, Qv)
		end
	end
	
	return equivalent_subtrees, equivalent_front_vertices, semi_equivalent_subtrees, semi_equivalent_front_vertices, isolated_isomorphic_subtrees, isolated_isomorphic_vertices
end

# ╔═╡ 4a94a37e-7f2b-11eb-110a-09cefcf5da72
function process_leaf(v, bcg::BlockCopolymerGraph)
	equivalent_subtrees = []
	equivalent_front_vertices = []
	semi_equivalent_subtrees = []
	semi_equivalent_front_vertices = []
	
	vnext = find_leaf_neighbor(v, bcg)
	subtree, vmap = induced_subgraph(bcg.graph, [v, vnext])

	n = count_isomorphic_subtree(subtree, vmap, bcg)
	(n == 1) && return ([], [], [], [])
	
	subtrees, front_vertices = all_isomorphic_subtree(subtree, vnext, vmap, bcg)
	isomorphic_subtrees_list = [subtrees]
	isomorphic_vertices_list = [front_vertices]
	
# 	subtree_groups, vertex_groups, ntrees = group_isomorphic_subtrees(subtrees, front_vertices, bcg)
# 	es, ev, ss, sv, is, iv = process_isomorphic_subtree_groups(subtree_groups, vertex_groups, ntrees)
# 	es, ev, ss, sv = process_isolated_isomorphic_subtrees(first(is), first(iv), bcg)
# 	es, ev, is, iv = process_equivalent_subtrees(first(es), first(ev), bcg)
	
# 	subtree_groups, vertex_groups, ntrees = group_isomorphic_subtrees(is, iv, bcg)
# 	es, ev, ss, sv, is, iv = process_isomorphic_subtree_groups(subtree_groups, vertex_groups, ntrees)
	
	
	while length(isomorphic_subtrees_list) > 0
		equiv_subtrees_list = []
		equiv_vertices_list = []
		for i in 1:length(isomorphic_subtrees_list)
			isomorphic_subtrees = isomorphic_subtrees_list[i]
			isomorphic_vertices = isomorphic_vertices_list[i]
			@show isomorphic_subtrees
			@show isomorphic_vertices
			sg, vg, ntrees = group_isomorphic_subtrees(isomorphic_subtrees, isomorphic_vertices, bcg)
			@show sg
			@show vg
			@show ntrees
			es, ev, ss, sv, is, iv = process_isomorphic_subtree_groups(sg, vg, ntrees)
			@show is
			@show iv
			append!(equivalent_subtrees, es)
			append!(equivalent_front_vertices, ev)
			append!(semi_equivalent_subtrees, ss)
			append!(semi_equivalent_front_vertices, sv)
			for j in 1:length(is)
				isolated_subtrees = is[j]
				isolated_vertices = iv[j]
				@show isolated_subtrees
				@show isolated_vertices
				es, ev, ss, sv = process_isolated_isomorphic_subtrees(
									isolated_subtrees, 
									isolated_vertices,
									bcg)
				append!(semi_equivalent_subtrees, ss)
				append!(semi_equivalent_front_vertices, sv)
				append!(equiv_subtrees_list, es)
				append!(equiv_vertices_list, ev)
			end
		end

		isomorphic_subtrees_list = []
		isomorphic_vertices_list = []
		for i in 1:length(equiv_subtrees_list)
			equiv_subtrees = equiv_subtrees_list[i]
			equiv_vertices = equiv_vertices_list[i]
			es, ev, is, iv = process_equivalent_subtrees(equiv_subtrees, equiv_vertices, bcg)
			(length(es) > 0) && push!(equivalent_subtrees, es)
			(length(ev) > 0) && push!(equivalent_front_vertices, ev)
			(length(is) > 0) && push!(isomorphic_subtrees_list, is)
			(length(iv) > 0) && push!(isomorphic_vertices_list, iv)
		end
		@show equivalent_subtrees
		@show equivalent_front_vertices
		@show isomorphic_subtrees_list
		@show isomorphic_vertices_list
	end
	
	return equivalent_subtrees, equivalent_front_vertices, semi_equivalent_subtrees, semi_equivalent_front_vertices
end

# ╔═╡ aad7f0fa-824f-11eb-01e0-bb56fc01cdb3
process_leaf(7, starg)

# ╔═╡ f6989fd0-8236-11eb-1bbe-a7c5e8a3ef62
process_leaf(9, BlockCopolymerGraph(chainAB3A6()))

# ╔═╡ 85b88436-81ae-11eb-0f96-1d4851bbbcfe
process_leaf(11, BlockCopolymerGraph(branchAB()))

# ╔═╡ 55d90414-823a-11eb-3523-d1683c46b6fe
process_leaf(9, BlockCopolymerGraph(branchAB()))

# ╔═╡ 6ae737b0-823a-11eb-25f7-8f85d648122f
process_leaf(6, BlockCopolymerGraph(chainAB3A3()))

# ╔═╡ 7669f144-8248-11eb-0a39-c7fe01ec3062
branchg2 = BlockCopolymerGraph(branchAB2());

# ╔═╡ 87e6de28-8248-11eb-3a1a-ed215a26b75a
process_leaf(16, branchg2)

# ╔═╡ c4109672-8249-11eb-0463-dd7319bded51
process_leaf(18, branchg2)

# ╔═╡ cf27af50-8249-11eb-378c-2dbe6b11846d
process_leaf(14, branchg2)

# ╔═╡ b879a606-824a-11eb-355e-4b4d61722144
process_leaf(2, chainABAg)

# ╔═╡ c11b46e2-8163-11eb-1030-abb9389714c4
function process_leaf_first_step(v, bcg::BlockCopolymerGraph)
	vnext = find_leaf_neighbor(v, bcg)
	subtree, vmap = induced_subgraph(bcg.graph, [v, vnext])
	n = count_isomorphic_subtree(subtree, vmap, bcg)
	(n == 1) && return nothing
	
	subtrees, front_vertex = all_isomorphic_subtree(subtree, vnext, vmap, bcg)
	equiv_subtrees = []
	for (vertices, vf) in zip(subtrees, front_vertex)
		(v ∈ vertices || vnext ∈ vertices) && (push!(equiv_subtrees, vertices); continue)
		subtreex, vmapx = induced_subgraph(bcg.graph, vertices)
		subtreey, vmapy = induced_subtree(vnext, vf, vmap, vmapx, bcg)
		rvmapy = reverse_vmap(vmapy)
		is_symmetric_tree(subtreey, rvmapy[vnext], rvmapy[vf], vmapy, bcg) && push!(equiv_subtrees, vertices)
	end
	
	return equiv_subtrees, front_vertex
end

# ╔═╡ 068b3182-80c0-11eb-280d-e3dd9abefec4
process_leaf_first_step(10, branchg)

# ╔═╡ 92e78932-80cf-11eb-1289-a5133bf1eea6
process_leaf_first_step(1, bcg)

# ╔═╡ cb4ddaca-80dd-11eb-0037-1f93de1ebab6
process_leaf_first_step(8, bcg)

# ╔═╡ 87e9c31a-8629-11eb-3959-970814e44072
# process_semi_equivalent_subtrees_one_step([[12,1,13,2],[11,19,20,9]], [2, 9], branchg2)

# ╔═╡ 453beccc-8639-11eb-0c99-1153564bd10c
process_leaf(14, branchg2)

# ╔═╡ 96f09382-8632-11eb-29b0-af9c02561b18
# process_semi_equivalent_subtrees_one_step([[11,3,2],[12,5,4]], [2, 4], branchg)

# ╔═╡ baf2a2d4-8572-11eb-36f3-d9603492bc5b
"""
Assumptions:

subtree1, subtree2: a list of vertices indexed in the original bcg.graph
v1, v2: front vertices of subtree1 and subtree2 indexed in the original bcg.graph
"""
function is_isomorphic_tree(subtree1, v1, subtree2, v2, bcg::BlockCopolymerGraph)
	tree1, vmap1 = induced_subgraph(bcg.graph, subtree1)
	tree2, vmap2 = induced_subgraph(bcg.graph, subtree2)
	rvmap1 = reverse_vmap(vmap1)
	rvmap2 = reverse_vmap(vmap2)
	dfs1 = dfs_tree(tree1, rvmap1[v1])
	dfs2 = dfs_tree(tree2, rvmap2[v2])
	return has_isomorph(dfs1, dfs2;
				vertex_relation=(v1,v2)->equivalent_blockend(v1,v2,vmap1,vmap2,bcg),
				edge_relation=(e1,e2)->equivalent_block(e1,e2,vmap1,vmap2,bcg))
end

# ╔═╡ 25e58756-8603-11eb-371f-99d5463fe797
is_isomorphic_tree([2,1,9,10], 2, [4,6,13,14], 4, branchg)

# ╔═╡ a2a13f06-858b-11eb-2a14-df05785b3b08
"""
Group according to number of branches emanating from the front vertices.
"""
function group_semi_equivalent_subtrees(subtrees, front_vertices, bcg::BlockCopolymerGraph)
	tree_groups = []
	vertex_groups = []
	
	nbranch = [degree(bcg.graph, i) for i in front_vertices]
	nbranch_unique = nbranch |> unique |> sort
	
	for n in nbranch_unique
		tgroup = []
		vgroup = []
		for i in 1:length(subtrees)
			(nbranch[i] == n) && push!(tgroup, subtrees[i])
			(nbranch[i] == n) && push!(vgroup, front_vertices[i])
		end
		push!(tree_groups, tgroup)
		push!(vertex_groups, vgroup)
	end
	
	return tree_groups, vertex_groups
end

# ╔═╡ 2db76df4-864d-11eb-0e6f-dd7c28703d70
function process_semi_equivalent_subtree_group(::Val{1}, subtrees, front_vertices, bcg::BlockCopolymerGraph)
	return [subtrees], [front_vertices], [], []
end

# ╔═╡ 34ee6a3c-864f-11eb-2775-fbd353ce8544
"""
This is for N-element group where N > 2.
"""
function process_semi_equivalent_subtree_group(::Val{N}, subtrees, front_vertices, bcg::BlockCopolymerGraph) where N
	semi_equivalent_subtrees = []
	semi_equivalent_vertices = []
	isomorphic_subtrees = []
	isomorphic_front_vertices = []
	
	for i in 1:length(subtrees)
		subtree1 = subtrees[i]
		v1 = front_vertices[i]
		for j in (i+1):length(subtrees)
			subtree2 = subtrees[j]
			v2 = front_vertices[j]
			ss, sv, is, iv = process_semi_equivalent_subtree_group(Val(2), [subtree1, subtree2], [v1, v2], bcg)
			(length(ss) > 0) && append!(semi_equivalent_subtrees, ss)
			(length(sv) > 0) && append!(semi_equivalent_vertices, sv)
			(length(is) > 0) && append!(isomorphic_subtrees, is)
			(length(iv) > 0) && append!(isomorphic_front_vertices, iv)
		end
	end
	
	return semi_equivalent_subtrees, semi_equivalent_vertices, isomorphic_subtrees, isomorphic_front_vertices
end

# ╔═╡ 636a4996-858c-11eb-3a29-63ef42d57113
function process_semi_equivalent_subtree_group(subtrees, front_vertices, bcg::BlockCopolymerGraph)
	m = length(subtrees)
	return process_semi_equivalent_subtree_group(Val(m), subtrees, front_vertices, bcg)
end

# ╔═╡ 4948c96c-8601-11eb-2d56-ede3eb0c6966
function group_isomorphic_branches(branches, front_vertex, bcg::BlockCopolymerGraph)
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
			if is_isomorphic_tree(subtree1, front_vertex, subtree2, front_vertex, bcg)
				push!(bgroup, subtree2)
				grouped[j] = true
			end
		end
		grouped[i] = true
		push!(branch_groups, bgroup)
	end
	
	return branch_groups
end

# ╔═╡ 20f996ba-864e-11eb-0802-41f97dc221d4
function process_semi_equivalent_subtree_group(::Val{2}, subtrees, front_vertices, bcg::BlockCopolymerGraph)
	semi_equivalent_subtrees = []
	semi_equivalent_vertices = []
	isomorphic_subtrees = []
	isomorphic_front_vertices = []
	
	nbranch = degree(bcg.graph, first(front_vertices))
	
	subtree1, subtree2 = deepcopy(subtrees)
	v1, v2 = front_vertices
	v1pre = nothing
	for vx in neighbors(bcg.graph, v1)
		(vx ∈ subtree1) && (v1pre = vx; break)
	end
	v2pre = nothing
	for vx in neighbors(bcg.graph, v2)
		(vx ∈ subtree2) && (v2pre = vx; break)
	end

	e1, e2 = nothing, nothing
	if nbranch == 2
		for vx in neighbors(bcg.graph, v1)
			(vx ∉ subtree1) && (e1 = Edge(vx, v1))
		end
		for vx in neighbors(bcg.graph, v2)
			(vx ∉ subtree2) && (e2 = Edge(vx, v2))
		end
		if equivalent_block(e1, e2, bcg)
			push!(subtree1, src(e1))
			push!(subtree2, src(e2))
			push!(isomorphic_subtrees, [subtree1, subtree2])
			push!(isomorphic_front_vertices, [src(e1), src(e2)])
		else
			push!(semi_equivalent_subtrees, subtrees)
			push!(semi_equivalent_vertices, front_vertices)
		end
		return semi_equivalent_subtrees, semi_equivalent_vertices, isomorphic_subtrees, isomorphic_front_vertices
	end

	# when nbranch > 2
	isomorphic_branches1 = []
	isomorphic_vertices1 = []
	isomorphic_branches2 = []
	isomorphic_vertices2 = []
	branch12 = nothing  # branch detached from v1 containing v2
	branch21 = nothing  # branch detached from v2 containing v1
	branches1 = induced_subtree(v1, bcg)
	newbranches1 = []
	for (branchg, vmap) in branches1
		(v1pre ∈ vmap) && continue
		(v2 ∈ vmap) && (branch12 = vmap; continue)
		push!(newbranches1, vmap)
	end
	branches2 = induced_subtree(v2, bcg)
	newbranches2 = []
	for (branchg, vmap) in branches2
		(v2pre ∈ vmap) && continue
		(v1 ∈ vmap) && (branch21 = vmap; continue)
		push!(newbranches2, vmap)
	end
	bgroups1 = group_isomorphic_branches(newbranches1, v1, bcg)
	bgroups2 = group_isomorphic_branches(newbranches2, v2, bcg)
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
				append!(isomorphic_vertices1, fill(v1, p))
				append!(isomorphic_branches2, bgroup2[1:p])
				append!(isomorphic_vertices2, fill(v2, p))
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
	return semi_equivalent_subtrees, semi_equivalent_vertices, isomorphic_subtrees, isomorphic_front_vertices
end

# ╔═╡ 82f617a6-84c5-11eb-10ce-e7b668973a01
function process_semi_equivalent_subtrees_one_step(subtrees, front_vertices, bcg::BlockCopolymerGraph)
	semi_equivalent_subtrees = []
	semi_equivalent_vertices = []
	isomorphic_subtrees = []
	isomorphic_front_vertices = []

	tree_groups, vertex_groups = group_semi_equivalent_subtrees(subtrees, front_vertices, bcg)
	
	for i in length(tree_groups)
		ss, sv, is, iv = process_semi_equivalent_subtree_group(tree_groups[i], vertex_groups[i], bcg)
		append!(semi_equivalent_subtrees, ss)
		append!(semi_equivalent_vertices, sv)
		append!(isomorphic_subtrees, is)
		append!(isomorphic_front_vertices, iv)
	end
	
	return semi_equivalent_subtrees, semi_equivalent_vertices, isomorphic_subtrees, isomorphic_front_vertices
end

# ╔═╡ 0da9b6f0-8633-11eb-15de-957a9d39260c
function process_semi_equivalent_subtrees(subtrees, front_vertices, bcg::BlockCopolymerGraph)
	semi_equivalent_subtrees = []
	semi_equivalent_vertices = []
	
	isomorphic_subtrees = [subtrees]
	isomorphic_front_vertices = [front_vertices]
	
	while length(isomorphic_subtrees) > 0
		tree_groups = []
		vertex_groups = []
		for i in 1:length(isomorphic_subtrees)
			trees = isomorphic_subtrees[i]
			vertices = isomorphic_front_vertices[i]
			tgroups, vgroups = group_semi_equivalent_subtrees(trees, vertices, bcg)
			append!(tree_groups, tgroups)
			append!(vertex_groups, vgroups)
		end
		
		isomorphic_subtrees = []
		isomorphic_front_vertices = []
		for i in 1:length(tree_groups)
			ss, sv, is, iv = process_semi_equivalent_subtree_group(tree_groups[i], vertex_groups[i], bcg)
			@show ss
			(length(ss) > 0) && append!(semi_equivalent_subtrees, ss)
			(length(sv) > 0) && append!(semi_equivalent_vertices, sv)
			(length(is) > 0) && append!(isomorphic_subtrees, is)
			(length(iv) > 0) && append!(isomorphic_front_vertices, iv)
		end
	end
	
	return semi_equivalent_subtrees, semi_equivalent_vertices
end

# ╔═╡ 84d66802-8635-11eb-3540-63c9f4ca8b79
process_semi_equivalent_subtrees([[14,3], [18,10]], [3, 10], branchg2)

# ╔═╡ 5f32ea4e-8635-11eb-29fe-655dc625f746
process_semi_equivalent_subtrees([[11,3],[12,5]], [3, 5], branchg)

# ╔═╡ 4ad3ced0-8638-11eb-3d3b-316e2459f55a
process_semi_equivalent_subtrees([[10,1],[13,6]], [1, 6], branchg)

# ╔═╡ a7490832-865c-11eb-1a62-c78bd869534a
process_semi_equivalent_subtrees([[14,3], [18,10], [15,5]], [3, 10, 5], branchg2)

# ╔═╡ b7339410-8605-11eb-1d20-6f1080487157
group_isomorphic_branches([[2,1], [2,3,7,11], [2,4,8,12], [2,5,9,13], [2,6,10,14]], 2, starg)

# ╔═╡ 0c65112c-8606-11eb-1dd1-7b6a45f44523
starg

# ╔═╡ fa935b60-8572-11eb-39fd-f3616623f688
"""
`subtrees`: a set of semi-equivalent subtrees found by the `process_leaf` method. Note that all subtrees in this set are isomorphic. At present, each subtree is merely a set of vertices in the original graph.
`front_vertices`: corresponding front vertices
"""
function maximize_semi_equivalent_subtree(subtrees, front_vertices, bcg::BlockCopolymerGraph)
	semi_equivalent_subtrees = []
	semi_equivalent_vertices = []
	for i in 1:length(subtrees)
		for j in (i+1):length(subtrees)
			subtree1 = subtrees[i]
			v1 = front_vertices[i]
			subtree2 = subtrees[j]
			v2 = front_vertices[j]
			
		end
	end
end

# ╔═╡ cc9293a4-865e-11eb-1685-fff5a99dce6c
function starABCDO()
	sA = KuhnSegment(:A)
	sB = KuhnSegment(:B)
	eb = branchpoints(4)
	
	O1 = PolymerBlock(:O1, sA, 0.05, eb[1], FreeEnd(:O1))
	O2 = PolymerBlock(:O2, sA, 0.05, eb[2], FreeEnd(:O2))
	O3 = PolymerBlock(:O3, sA, 0.05, eb[3], FreeEnd(:O3))
	O4 = PolymerBlock(:O4, sA, 0.05, eb[4], FreeEnd(:O4))
	
	A2 = PolymerBlock(:A2, sB, 0.04, eb[2], FreeEnd(:A2))
	A3 = PolymerBlock(:A3, sB, 0.04, eb[3], FreeEnd(:A3))
	A4 = PolymerBlock(:A4, sB, 0.04, eb[4], FreeEnd(:A4))
	
	B2 = PolymerBlock(:B2, sB, 0.03, eb[2], FreeEnd(:B2))
	B3 = PolymerBlock(:B3, sB, 0.03, eb[3], FreeEnd(:B3))
	
	C1 = PolymerBlock(:C1, sA, 0.06, eb[1], FreeEnd(:C1))
	C2 = PolymerBlock(:C2, sA, 0.06, eb[2], FreeEnd(:C2))
	C4 = PolymerBlock(:C4, sA, 0.06, eb[4], FreeEnd(:C4))
	
	D3 = PolymerBlock(:D3, sB, 0.1, eb[3], FreeEnd(:D3))
	D4 = PolymerBlock(:D4, sB, 0.1, eb[4], FreeEnd(:D4))
	
	X12 = PolymerBlock(:X12, sA, 0.08, eb[1], eb[2])
	X13 = PolymerBlock(:X13, sB, 0.08, eb[1], eb[3])
	X14 = PolymerBlock(:X14, sA, 0.08, eb[1], eb[4])
	
	return BlockCopolymer(:ABCDO, [O1, O2, O3, O4, A2, A3, A4, B2, B3, C1, C2, C4, D3, D4, X12, X13, X14])
end

# ╔═╡ 0cfb898a-8661-11eb-01fd-13397e2d299f
starABCDOg = BlockCopolymerGraph(starABCDO())

# ╔═╡ 2f8b1880-8661-11eb-2234-39294f5b958d
process_leaf(2, starABCDOg)

# ╔═╡ 80ea4a02-8661-11eb-34ab-ed1846867e84
process_semi_equivalent_subtrees([[2,1], [4,3], [6,5], [8,7]], [1, 3, 5, 7], starABCDOg)

# ╔═╡ 87ee1990-8662-11eb-3739-cd3e623d4602
process_leaf(12, starABCDOg)

# ╔═╡ a0f7deee-8662-11eb-11a5-3f5f82adfa3f
process_semi_equivalent_subtrees([[12,3], [13,5]], [3, 5], starABCDOg)

# ╔═╡ 3e293a52-8661-11eb-30f6-73dceaa28049
starABCDOg.node2free

# ╔═╡ Cell order:
# ╠═5b913a92-7f20-11eb-1f10-b3059d07c947
# ╠═71fa9742-7f20-11eb-369a-ef11d25edb97
# ╠═088b1254-7f21-11eb-0ab9-234397a05e9d
# ╠═79fac120-7f2e-11eb-2007-a3bde9fad3f7
# ╠═f20867c8-8017-11eb-2b31-9b24a0235101
# ╠═a3376b3a-7f2e-11eb-3c1b-b1fdf40461a5
# ╠═36cfce44-824a-11eb-367e-9b18ea4df6bf
# ╠═997c616a-824a-11eb-0a7a-e7dcb8a66df0
# ╠═110d5194-824b-11eb-19c4-4991b80b89bc
# ╠═77040be2-7f20-11eb-2132-7d12bbe3751e
# ╠═88554b48-7f36-11eb-15b2-6d5cb42186d1
# ╠═8c2f0162-824e-11eb-3a4d-7110b677cbea
# ╠═7df1a45a-824f-11eb-1538-23b1c5782c75
# ╠═993a0e20-824f-11eb-2963-fd2dcbf63882
# ╠═aad7f0fa-824f-11eb-01e0-bb56fc01cdb3
# ╠═0e1c3dee-7fc9-11eb-0b44-ad794b13a6ca
# ╠═f0e9ac34-7fc8-11eb-27a3-b95c83218a18
# ╠═a166207e-7fc9-11eb-2ada-77f2740283b4
# ╠═748a414a-7fc9-11eb-0250-6312741527d7
# ╠═f1b33da2-7fc9-11eb-1f75-bf245a82d0ea
# ╠═469c1d94-7fb8-11eb-307b-a13577067b52
# ╠═1dd95408-8247-11eb-2a29-09164b734f9b
# ╠═53e6ef00-8248-11eb-1837-ad033d85c8f8
# ╠═277a4454-7fb8-11eb-08f0-c308ad0fcd32
# ╠═8d194a3c-7f20-11eb-105a-0b91f6a69bc1
# ╠═6d608ed6-7f49-11eb-25ef-272fdae43007
# ╠═cd4c77a2-824d-11eb-2f65-998feeb8cf5e
# ╠═dcc4cf54-824d-11eb-04ed-138f3813ab5e
# ╠═e31ec7ec-824d-11eb-0101-714661027e0a
# ╠═ec653e1a-824d-11eb-290c-6dead284203d
# ╠═9567dd66-7f20-11eb-0bda-315ad04a42ff
# ╠═a8d5721e-7f20-11eb-1c3f-fd4ffea1637e
# ╠═ba365686-7f20-11eb-22c6-0782aaa5897d
# ╠═f506bc06-800b-11eb-197a-fb31c6d4dd22
# ╠═0278037c-800c-11eb-10c9-21130c78671c
# ╠═1bac3db8-800c-11eb-27ad-d10ce333659f
# ╠═28d4a0e6-800c-11eb-3010-492eada4bce8
# ╠═46ebdeb6-800c-11eb-203f-bb21e582625e
# ╠═711caa16-8085-11eb-1a2c-7133fd44aa62
# ╠═072042e6-8088-11eb-2c38-bd19e1578001
# ╠═c18f4e4a-8085-11eb-013b-97978173c7cc
# ╠═0d6fec00-8088-11eb-36d5-719ab414513e
# ╠═1989d0c2-7f21-11eb-13a6-2f87902db271
# ╠═84157cf4-7f29-11eb-36cf-59dea7d7f5b6
# ╠═c3493360-7f20-11eb-0a22-637ae33e5833
# ╠═5e0481d4-7f21-11eb-1b31-2f6dd8b2738c
# ╠═15274aa6-7f4d-11eb-09c8-5dda2485fe85
# ╠═e08c53ce-7f51-11eb-1e85-1386c9a2757d
# ╠═c04a2b14-8002-11eb-38ce-218064340191
# ╠═5df84cf8-7f3e-11eb-2840-b98a8db5d97c
# ╠═fed4c060-8002-11eb-3e86-e55799b0498d
# ╠═ad015e6c-7fba-11eb-2a22-63b0cab369d3
# ╠═372af4ac-8003-11eb-0b48-b5372458c90b
# ╠═f51febe4-8598-11eb-3f19-7b2888b0723a
# ╠═9791b2aa-7ff3-11eb-26b0-1ba7eaba36ed
# ╠═7fbcff20-7ff5-11eb-37a1-651728486d23
# ╠═98e92396-7ff2-11eb-3782-3f5616a2434c
# ╠═60e61312-8003-11eb-3084-f357cd092257
# ╠═4b1029e2-801c-11eb-04f6-135d3cf91f6b
# ╠═b49fadaa-801d-11eb-01af-1be437348d56
# ╠═9ebd2494-801e-11eb-276f-93ae1af1b013
# ╠═ca86285c-801c-11eb-2a1b-0ff444ed0721
# ╠═f336d56a-80a0-11eb-2056-6d2cdf79d46a
# ╠═b2984fbe-7ff3-11eb-0e1c-dbd158dd04ae
# ╠═3cd857c2-801d-11eb-3591-dd51a1c57699
# ╠═c63d92fc-801d-11eb-1a6b-fdd23bde076a
# ╠═d060a24a-801f-11eb-06ae-73cb7c729035
# ╠═555b5956-8082-11eb-0b1e-c9eef34c041e
# ╠═a5898600-8082-11eb-24dc-df23129e37c4
# ╠═2b3fc368-8083-11eb-3504-9b5892365b45
# ╠═8344b206-8085-11eb-39c7-113ca193e89d
# ╠═d4f01daa-8087-11eb-0dc2-974351888dd7
# ╠═aba16426-8083-11eb-055a-6500d0de28de
# ╠═d99eb46a-8082-11eb-2bf0-cdbe6f0aacc9
# ╠═a2c40156-7f52-11eb-27f2-cfee4296f1b2
# ╠═edd548a2-7fbc-11eb-0956-cfb97b1c695d
# ╠═f0fce2f6-7fb7-11eb-1ece-b1e2494058db
# ╠═aaeba326-7fba-11eb-116e-296d4a438f04
# ╠═784a7912-7fd3-11eb-2591-17d4b1eaa347
# ╠═8aa193c0-7fd3-11eb-156b-c12d122a18c3
# ╠═b878ddaa-7fd4-11eb-06ef-3bac1ac83c38
# ╠═83c1e328-7fe2-11eb-38cf-116d4d032ce9
# ╠═65aece8a-7fd3-11eb-3978-27fa9683e21b
# ╠═0e7d8880-7fe2-11eb-1f77-754c4263a0b5
# ╠═22a8a6b2-7fe2-11eb-2037-55a1a1e060f9
# ╠═a738558a-8019-11eb-35a1-991259834d46
# ╠═1db19018-801b-11eb-152f-01608867a266
# ╠═89a27b48-801b-11eb-0a17-414e085d11d0
# ╠═cb84ef14-801b-11eb-0a4e-e5268c3a9965
# ╠═83eba7e2-7f21-11eb-2685-d781ba1b10a3
# ╠═97517fb4-7f21-11eb-0547-0963ab9ca934
# ╠═814c4f0c-7f22-11eb-13af-fb3eb53810b3
# ╠═053fe0b2-7f4d-11eb-09b8-41ba9a4d6387
# ╠═2c395a22-7f4d-11eb-2435-c9660719b17a
# ╠═319a9c80-7f2d-11eb-1a0c-bbc78da2540a
# ╠═016996e4-7f2e-11eb-3301-6da3c7b996ce
# ╠═0161b5be-801b-11eb-2317-c539091f1e8d
# ╠═a9bd1960-7f38-11eb-1c7f-cd3105f49480
# ╠═ae82e0d8-7fe4-11eb-3655-41da2101eb5a
# ╠═33c8dfa8-8018-11eb-26a8-4b7c0122ee31
# ╠═dd1e0226-8019-11eb-14e7-8767d9361a82
# ╠═5025a718-7f2c-11eb-2779-a5bad86fd6a0
# ╠═dc1bb0d8-7fe4-11eb-1778-5dd086ed216a
# ╠═54d6310a-8018-11eb-0468-2b1e98c5a6a1
# ╠═d783467e-7f2d-11eb-3266-498d9ed39e5f
# ╠═86bde198-801a-11eb-3412-cbc6fb494aec
# ╠═90a74d98-801a-11eb-2ecd-19f4b2ed9e72
# ╠═b8ea5c08-7f22-11eb-002d-e923402835ef
# ╠═f66615e4-7f28-11eb-3bd7-8f123f3b8837
# ╠═4b24db38-7f29-11eb-0164-7fa310866165
# ╠═61362204-7f2b-11eb-00cc-dfa70199b912
# ╠═f3d9ee22-80ba-11eb-3ecc-a54d00a89045
# ╠═cbaca3e4-80be-11eb-3938-7133e143f92e
# ╠═c0de8e22-80bf-11eb-039e-054f2df66301
# ╠═1845fcee-823d-11eb-2550-015a785a6453
# ╠═71ed8db8-823d-11eb-1486-c32cd5de0c38
# ╠═e012cbd2-8191-11eb-3ff1-69c1364ad19c
# ╠═7c4d201e-818b-11eb-2d02-218a6ed92a49
# ╠═3c62ebe6-8198-11eb-3255-dfbaa392145f
# ╠═a719a518-81ae-11eb-3935-a3ef2903b1cf
# ╠═7668845c-8217-11eb-1c04-a1cb72e69815
# ╠═4a94a37e-7f2b-11eb-110a-09cefcf5da72
# ╠═f6989fd0-8236-11eb-1bbe-a7c5e8a3ef62
# ╠═85b88436-81ae-11eb-0f96-1d4851bbbcfe
# ╠═55d90414-823a-11eb-3523-d1683c46b6fe
# ╠═6ae737b0-823a-11eb-25f7-8f85d648122f
# ╠═7669f144-8248-11eb-0a39-c7fe01ec3062
# ╠═87e6de28-8248-11eb-3a1a-ed215a26b75a
# ╠═c4109672-8249-11eb-0463-dd7319bded51
# ╠═cf27af50-8249-11eb-378c-2dbe6b11846d
# ╠═b879a606-824a-11eb-355e-4b4d61722144
# ╠═c11b46e2-8163-11eb-1030-abb9389714c4
# ╠═068b3182-80c0-11eb-280d-e3dd9abefec4
# ╠═92e78932-80cf-11eb-1289-a5133bf1eea6
# ╠═cb4ddaca-80dd-11eb-0037-1f93de1ebab6
# ╠═82f617a6-84c5-11eb-10ce-e7b668973a01
# ╠═0da9b6f0-8633-11eb-15de-957a9d39260c
# ╠═87e9c31a-8629-11eb-3959-970814e44072
# ╠═84d66802-8635-11eb-3540-63c9f4ca8b79
# ╠═453beccc-8639-11eb-0c99-1153564bd10c
# ╠═96f09382-8632-11eb-29b0-af9c02561b18
# ╠═5f32ea4e-8635-11eb-29fe-655dc625f746
# ╠═4ad3ced0-8638-11eb-3d3b-316e2459f55a
# ╠═baf2a2d4-8572-11eb-36f3-d9603492bc5b
# ╠═25e58756-8603-11eb-371f-99d5463fe797
# ╠═a2a13f06-858b-11eb-2a14-df05785b3b08
# ╠═2db76df4-864d-11eb-0e6f-dd7c28703d70
# ╠═20f996ba-864e-11eb-0802-41f97dc221d4
# ╠═a7490832-865c-11eb-1a62-c78bd869534a
# ╠═34ee6a3c-864f-11eb-2775-fbd353ce8544
# ╠═636a4996-858c-11eb-3a29-63ef42d57113
# ╠═4948c96c-8601-11eb-2d56-ede3eb0c6966
# ╠═b7339410-8605-11eb-1d20-6f1080487157
# ╠═0c65112c-8606-11eb-1dd1-7b6a45f44523
# ╠═fa935b60-8572-11eb-39fd-f3616623f688
# ╠═cc9293a4-865e-11eb-1685-fff5a99dce6c
# ╠═0cfb898a-8661-11eb-01fd-13397e2d299f
# ╠═2f8b1880-8661-11eb-2234-39294f5b958d
# ╠═80ea4a02-8661-11eb-34ab-ed1846867e84
# ╠═87ee1990-8662-11eb-3739-cd3e623d4602
# ╠═a0f7deee-8662-11eb-11a5-3f5f82adfa3f
# ╠═3e293a52-8661-11eb-30f6-73dceaa28049
