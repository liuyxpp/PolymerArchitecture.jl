@testset "subtree.jl: Constructor" begin
    subtree = Subtree(abg.graph, [2,3], 3)
    @test subtree.vmap == [2,3]
    @test subtree.vi == 2
    @test subtree.v == 3

    subtree = Subtree(abg, [2,3], 3)
    @test subtree.vmap == [2,3]
    @test subtree.vi == 2
    @test subtree.v == 3
end

@testset "subtree.jl: front_edges" begin
    subtree1, subtree2 = induced_subtree(branchg, Edge(2,7))
    @test front_vertex(subtree1) == 2
    @test front_vertex(subtree2) == 7
    @test front_edges(subtree1) == [Edge(2,1), Edge(2,3)]
    @test front_edges(subtree2) == [Edge(7,8)]
end

@testset "subtree.jl: find_neighbors" begin
    vs = find_neighbors(abg.graph, 2, 1)
    @test vs == [3]
    vs = find_neighbors(abg.graph, 1, 2)
    @test length(vs) == 0
end

@testset "subtree.jl: induced_subtree - edge" begin
    subtree1, subtree2 = induced_subtree(abg, Edge(1,2))
    @test subtree1.v == 1
    @test subtree1.vmap == [1]
    @test subtree2.v == 2
    @test subtree2.vmap == [3,2]
end

@testset "subtree.jl: induced_subtree - vertex" begin
    subtree1, subtree2 = induced_subtree(abg, 2)
    @test subtree1.v == 2
    @test subtree1.vmap == [2,1]
    @test subtree2.v == 2
    @test subtree2.vmap == [2,3]

    subtree1, = induced_subtree(abg, 1)
    @test subtree1.v == 1
    @test subtree1.vmap == [3,1,2]
end

@testset "subtree.jl: induced_subtree - edge1, edge2" begin
    subtree, v1, v2 = induced_subtree(abg, Edge(1,2), Edge(2,3))
    @test subtree.v == 2
    @test subtree.vmap == [2]
    @test v1 == 2
    @test v2 == 2
end

@testset "subtree.jl: induced_subtree - subtree1, subtree2" begin
    subtree1 = Subtree(branchg, [10,9,1,2,3,11], 2)
    subtree2 = Subtree(branchg, [4,5,12,6,13,14], 4)
    subtree = induced_subtree(branchg, subtree1, subtree2)
    @test Set(subtree.vmap) == Set([7, 8, 2, 4])
end

@testset "subtree.jl: merge_subtrees" begin
    subtree1, subtree2 = induced_subtree(branchg, Edge(2,7))
    subtree = merge_subtrees(branchg, [subtree1, subtree2])
    @test Set(subtree.vmap) == Set(1:nv(branchg.graph))
end