@testset "equivalent.jl: equivalent_blockend" begin
    @test equivalent_blockend(branchg, 2, 4)
    @test !equivalent_blockend(branchg, 2, 12)

    subtree1, subtree2 = induced_subtree(branchg, Edge(2,7))
    @test equivalent_blockend(branchg, 6, 2, subtree1.vmap, subtree2.vmap)
    @test equivalent_blockend(branchg, 1, 5, subtree1.vmap, subtree2.vmap)
    @test !equivalent_blockend(branchg, 1, 5)
end

@testset "equivalent.jl: equivalent_block" begin
    @test equivalent_block(branchg, Edge(2,7), Edge(4,8))
    @test !equivalent_block(branchg, Edge(2,7),Edge(4,6))

    subtree1, subtree2 = induced_subtree(branchg, Edge(2,7))
    @test equivalent_block(branchg, Edge(4,5), Edge(3,4), subtree1.vmap, subtree2.vmap)
    @test !equivalent_block(branchg, Edge(4,5), Edge(2,3), subtree1.vmap, subtree2.vmap)
end