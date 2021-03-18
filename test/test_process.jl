@testset "equivalent.jl: group_isomorphic_subtrees" begin
    subtree = Subtree(branchg, [1,9], 1)
    isotrees = all_isomorphic_subtree(branchg, subtree)
    groups = group_isomorphic_subtrees(isotrees)
    @test length(groups) == 2
    group1, group2 = groups
    @test group1[1].vmap == [1, 9]
    @test group1[2].vmap == [1, 10]
    @test group2[1].vmap == [6, 13]
    @test group2[2].vmap == [6, 14]
end