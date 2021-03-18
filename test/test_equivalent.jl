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

@testset "equivalent.jl: is_symmetric_subtree" begin
    subtree1 = Subtree(branchg, [10,9,1,2], 2)
    subtree2 = Subtree(branchg, [4,6,13,14], 4)
    subtree = induced_subtree(branchg, subtree1, subtree2)
    @test is_symmetric_subtree(branchg, subtree, 2, 4)

    subtree1 = Subtree(semibranchg, [10,9,1,2], 2)
    subtree2 = Subtree(semibranchg, [4,6,13,14], 4)
    subtree = induced_subtree(semibranchg, subtree1, subtree2)
    @test !is_symmetric_subtree(semibranchg, subtree, 2, 4)
end

@testset "equivalent.jl: is_isomorphic_subtree" begin
    subtree1 = Subtree(branchg, [10,9,1,2], 2)
    subtree2 = Subtree(branchg, [4,6,13,14], 4)
    @test is_isomorphic_subtree(branchg, subtree1, subtree2)

    subtree1 = Subtree(branchg, [10,9,1,2,3,11], 2)
    subtree2 = Subtree(branchg, [4,6,13,14,5,12], 4)
    @test is_isomorphic_subtree(branchg, subtree1, subtree2)
end

@testset "equivalent.jl: is_equivalent_subtree" begin
    subtree1 = Subtree(branchg, [10,9,1,2], 2)
    subtree2 = Subtree(branchg, [4,6,13,14], 4)
    @test is_equivalent_subtree(branchg, subtree1, subtree2)

    subtree1 = Subtree(semibranchg, [10,9,1,2], 2)
    subtree2 = Subtree(semibranchg, [4,6,13,14], 4)
    @test !is_equivalent_subtree(semibranchg, subtree1, subtree2)
end