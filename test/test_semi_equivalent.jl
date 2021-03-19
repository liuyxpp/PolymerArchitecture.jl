@testset "semi_equivalent.jl: group_isomorphic_branches" begin
    groups = group_isomorphic_branches(starAB4A8g, induced_subtree(starAB4A8g, 2))
    @test length(groups) == 3
    g1, g2, g3 = groups
    @test Set(g1[1].vmap) == Set([2, 1])
    @test Set(g2[1].vmap) == Set([7, 11, 2, 3])
    @test Set(g2[2].vmap) == Set([9, 13, 2, 5])
    @test Set(g3[1].vmap) == Set([8, 12, 2, 4])
    @test Set(g3[2].vmap) == Set([10, 14, 2, 6])
end

@testset "semi_equivalent.jl: group_semi_equivalent_subtrees" begin
    es, ss = process_leaf(branch2g, 13)
    groups = group_semi_equivalent_subtrees(branch2g, ss)
    @test length(groups) == 1
    group = first(groups)
    @test length(group) == 3
    @test Set(group[1].vmap) == Set([13, 1, 12])
    @test Set(group[2].vmap) == Set([16, 6, 17])
    @test Set(group[3].vmap) == Set([19, 11, 20])
end