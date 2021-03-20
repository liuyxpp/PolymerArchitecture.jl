branchpoints(n, prefix="EB") = [BranchPoint(Symbol(prefix*string(i))) for i in 1:n]

freeends(n, prefix="A") = [FreeEnd(Symbol(prefix*string(i))) for i in 1:n]

function branchAB(fB7=0.1, fB8=0.1)
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb = branchpoints(8)
    fe = freeends(6, "B")
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
    B7 = PolymerBlock(:B7, sB, fB7, eb[3], eb[4])
    B8 = PolymerBlock(:B8, sB, fB8, eb[5], eb[6])
    return BlockCopolymer(:AB3A3, [A1, A2, A3, A4, A5, B1, B2, B3, B4, B5, B6, B7, B8])
end

function branchAB2()
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb = branchpoints(11)
    fe = freeends(6, "B")
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

function starAB3A3()
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

function starAB3A6()
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

function starABCDO4()
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    sC = KuhnSegment(:C)
    sD = KuhnSegment(:D)
    sO = KuhnSegment(:O)
    eb = branchpoints(4)

    O1 = PolymerBlock(:O1, sO, 0.05, eb[1], FreeEnd(:O1))
    O2 = PolymerBlock(:O2, sO, 0.05, eb[2], FreeEnd(:O2))
    O3 = PolymerBlock(:O3, sO, 0.05, eb[3], FreeEnd(:O3))
    O4 = PolymerBlock(:O4, sO, 0.05, eb[4], FreeEnd(:O4))

    A2 = PolymerBlock(:A2, sA, 0.04, eb[2], FreeEnd(:A2))
    A3 = PolymerBlock(:A3, sA, 0.04, eb[3], FreeEnd(:A3))
    A4 = PolymerBlock(:A4, sA, 0.04, eb[4], FreeEnd(:A4))

    B2 = PolymerBlock(:B2, sB, 0.03, eb[2], FreeEnd(:B2))
    B3 = PolymerBlock(:B3, sB, 0.03, eb[3], FreeEnd(:B3))

    C1 = PolymerBlock(:C1, sC, 0.06, eb[1], FreeEnd(:C1))
    C2 = PolymerBlock(:C2, sC, 0.06, eb[2], FreeEnd(:C2))
    C4 = PolymerBlock(:C4, sC, 0.06, eb[4], FreeEnd(:C4))

    D3 = PolymerBlock(:D3, sD, 0.1, eb[3], FreeEnd(:D3))
    D4 = PolymerBlock(:D4, sD, 0.1, eb[4], FreeEnd(:D4))

    X12 = PolymerBlock(:X12, sA, 0.08, eb[1], eb[2])
    X13 = PolymerBlock(:X13, sB, 0.08, eb[1], eb[3])
    X14 = PolymerBlock(:X14, sA, 0.08, eb[1], eb[4])

    return BlockCopolymer(:ABCDO, [O1, O2, O3, O4, A2, A3, A4, B2, B3, C1, C2, C4, D3, D4, X12, X13, X14])
end

abg = diblock_chain() |> BlockCopolymerGraph

branchg = branchAB() |> BlockCopolymerGraph

semibranchg = branchAB(0.08, 0.12) |> BlockCopolymerGraph

branch2g = branchAB2() |> BlockCopolymerGraph

starAB3A3g = starAB3A3() |> BlockCopolymerGraph

starAB3A6g = starAB3A6() |> BlockCopolymerGraph

starAB4A8g = starAB4A8() |> BlockCopolymerGraph

starABCDOg = starABCDO() |> BlockCopolymerGraph

starABCDO4g = starABCDO4() |> BlockCopolymerGraph