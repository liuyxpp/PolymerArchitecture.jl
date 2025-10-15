using Test
using PolymerArchitecture: group_equivalent_blocks_traversal

# chains.jl has been included in runtests.jl

canon_group(g) = sort([(p[1], p[2]) for p in g])
canon_groups(v) = sort([canon_group(g) for g in v])

@testset "group_equivalent_blocks_traversal: branchg" begin
    got = group_equivalent_blocks_traversal(branchg)
    # expected from user description; order of entries not important
    # Order of pairs in each group is also not important
    expected = [
        [9 => 1, 10 => 1, 11 => 3, 12 => 5, 13 => 6, 14 => 6],
        [1 => 9, 1 => 10, 6 => 13, 6 => 14],
        [2 => 1, 4 => 6],
        [2 => 3, 4 => 5],
        [3 => 11, 5 => 12],
        [3 => 2, 5 => 4],
        [1 => 2, 6 => 4],
        [2 => 7, 4 => 8],
        [7 => 2, 8 => 4],
        [7 => 8, 8 => 7]
    ]
    @test canon_groups(got) == canon_groups(expected)
end

@testset "group_equivalent_blocks_traversal: semibranchg" begin
    got = group_equivalent_blocks_traversal(semibranchg)
    # expected from user description; order of entries not important
    # Order of pairs in each group is also not important
    expected = [
        [9=>1, 10=>1, 11=>3, 12=>5, 13=>6, 14=>6],
        [3=>2, 5=>4],
        [1=>9, 1=>10],
        [6=>13, 6=>14],
        [1=>2, 6=>4],
        [2=>1],
        [2=>3],
        [4=>5],
        [4=>6],
        [3=>11],
        [5=>12],
        [2=>7],
        [7=>2],
        [4=>8],
        [8=>4],
        [7=>8],
        [8=>7]
    ]
    @test canon_groups(got) == canon_groups(expected)
end
