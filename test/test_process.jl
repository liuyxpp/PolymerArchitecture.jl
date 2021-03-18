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

@testset "equivalent.jl: process_isomorphic_subtree_groups" begin
    subtree = Subtree(branchg, [1,9], 1)
    isotrees = all_isomorphic_subtree(branchg, subtree)
    groups = group_isomorphic_subtrees(isotrees)
    es, ss, is = process_isomorphic_subtree_groups(branchg, groups)
    @test length(es) == 0
    @test length(ss) == 0
    @test length(is) == 1
    new_isotrees = first(is)
    @test new_isotrees[1].vmap == [1, 9, 10]
    @test new_isotrees[2].vmap == [6, 13, 14]
end