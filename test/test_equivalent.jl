@testset "equivalent.jl: equivalent_blockend" begin
    @test equivalent_blockend(branchg, 2, 4)
    @test !equivalent_blockend(branchg, 2, 12)

    subtree1, subtree2 = induced_subtree(branchg, Edge(2,7))
    @test equivalent_blockend(branchg, 6, 2, subtree1.vmap, subtree2.vmap)
    @test equivalent_blockend(branchg, 1, 5, subtree1.vmap, subtree2.vmap)
    @test !equivalent_blockend(branchg, 1, 5)
end