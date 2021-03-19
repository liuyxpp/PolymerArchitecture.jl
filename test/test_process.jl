@testset "process.jl: group_isomorphic_subtrees" begin
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

@testset "process.jl: process_isomorphic_subtree_groups" begin
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

@testset "process.jl: process_isolated_isomorphic_subtrees" begin
    subtree = Subtree(branchg, [1,9], 1)
    isotrees = all_isomorphic_subtree(branchg, subtree)
    groups = group_isomorphic_subtrees(isotrees)
    es, ss, is = process_isomorphic_subtree_groups(branchg, groups)
    es, ss = process_isolated_isomorphic_subtrees(branchg, first(is))
    @test length(es) == 1
    @test length(ss) == 0
    etree1, etree2 = first(es)
    @test Set(etree1.vmap) == Set([1, 9, 10])
    @test Set(etree2.vmap) == Set([6, 13, 14])
end

@testset "process.jl: process_equivalent_subtrees" begin
    subtree = Subtree(branchg, [1,9], 1)
    isotrees = all_isomorphic_subtree(branchg, subtree)
    groups = group_isomorphic_subtrees(isotrees)
    es, ss, is = process_isomorphic_subtree_groups(branchg, groups)
    es, ss = process_isolated_isomorphic_subtrees(branchg, first(is))
    es, is = process_equivalent_subtrees(branchg, first(es))
    @test length(es) == 0
    @test length(is) == 2
    itree1, itree2 = is
    @test Set(itree1.vmap) == Set([1, 9, 10, 2])
    @test Set(itree2.vmap) == Set([6, 13, 14, 4])
end

@testset "process.jl: process_leaf" begin
    es, ss = process_leaf(branchg, 13)
    @test length(es) == 1
    @test length(ss) == 0
    etree1, etree2 = first(es)
    @test Set(etree1.vmap) == Set([4, 5, 12, 6, 13, 14, 8, 7])
    @test etree1.v == 7
    @test Set(etree2.vmap) == Set([2, 1, 9, 10, 3, 11, 7, 8])
    @test etree2.v == 8
end