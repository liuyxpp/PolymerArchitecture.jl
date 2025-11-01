using Polymer
using Polymer: branchpoints, freeends
using PolymerArchitecture

function linearABA(fA1=0.3, fA2=0.3, fB=0.4)
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb = branchpoints(2)
    fe = freeends(2)
    A1 = PolymerBlock(:A1, sA, fA1, eb[1], fe[1])
    A2 = PolymerBlock(:A2, sA, fA2, eb[2], fe[2])
    B = PolymerBlock(:B, sB, fB, eb[1], eb[2])
    return BlockCopolymer(:ABA, [A1, A2, B])
end

function branchAB(fB7=0.1, fB8=0.1)
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb = branchpoints(8)
    fe = freeends(6, "B")
    A1 = PolymerBlock(:A1, sA, 0.12, eb[1], eb[3])
    A2 = PolymerBlock(:A2, sA, 0.12, eb[2], eb[3])
    A3 = PolymerBlock(:A3, sA, 0.12, eb[6], eb[8])
    A4 = PolymerBlock(:A4, sA, 0.12, eb[6], eb[7])
    A5 = PolymerBlock(:A5, sA, 0.2, eb[4], eb[5])
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

"""
A - B - A - C1
    |   |
    C2  D
"""
function chainABCACD()
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    sC = KuhnSegment(:C)
    sD = KuhnSegment(:D)
    eb1 = BranchPoint(:EB1)
    eb2 = BranchPoint(:EB2)
    eb3 = BranchPoint(:EB3)
    A1 = PolymerBlock(:A1, sA, 0.2, FreeEnd(:A1), eb1)
    B = PolymerBlock(:B, sB, 0.3, eb1, eb2)
    A2 = PolymerBlock(:A2, sA, 0.3, eb2, eb3)
    C1 = PolymerBlock(:C1, sC, 0.1, eb3, FreeEnd(:C1))
    C2 = PolymerBlock(:C2, sC, 0.05, eb2, FreeEnd(:C2))
    D = PolymerBlock(:D, sD, 0.05, eb3, FreeEnd(:D))
    return BlockCopolymer(:ABCACD, [A1, B, A2, C1, C2, D])
end

function chainA6B6()
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb0 = BranchPoint(:EB0)
    eb1 = BranchPoint(:EB1)
    eb2 = BranchPoint(:EB2)
    eb3 = BranchPoint(:EB3)
    eb4 = BranchPoint(:EB4)
    A1 = PolymerBlock(:A1, sA, 0.1, eb0, FreeEnd(:A1))
    A2 = PolymerBlock(:A2, sA, 0.1, eb0, FreeEnd(:A2))
    A3 = PolymerBlock(:A3, sA, 0.05, eb1, FreeEnd(:A3))
    A4 = PolymerBlock(:A4, sA, 0.1, eb1, eb2)
    A5 = PolymerBlock(:A5, sA, 0.05, eb3, FreeEnd(:A5))
    A6 = PolymerBlock(:A6, sA, 0.1, eb3, eb4)
    B1 = PolymerBlock(:B1, sB, 0.1, eb0, eb1)
    B2 = PolymerBlock(:B2, sB, 0.05, eb1, FreeEnd(:B2))
    B3 = PolymerBlock(:B3, sB, 0.1, eb2, FreeEnd(:B3))
    B4 = PolymerBlock(:B4, sB, 0.1, eb2, eb3)
    B5 = PolymerBlock(:B5, sB, 0.05, eb2, FreeEnd(:B5))
    B6 = PolymerBlock(:B6, sB, 0.1, eb4, FreeEnd(:B6))
    return BlockCopolymer(:A6B6, [A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6])
end

function starA2B2()
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb0 = BranchPoint(:EB0)
    A1 = PolymerBlock(:A1, sA, 0.25, eb0, FreeEnd(:A1))
    A2 = PolymerBlock(:A2, sA, 0.25, eb0, FreeEnd(:A2))
    B1 = PolymerBlock(:B1, sB, 0.25, eb0, FreeEnd(:B1))
    B2 = PolymerBlock(:B2, sB, 0.25, eb0, FreeEnd(:B2))
    return BlockCopolymer(:A2B2, [A1, A2, B1, B2])
end

function starA3()
    sA = KuhnSegment(:A)
    eb0 = BranchPoint(:EB0)
    A1 = PolymerBlock(:A1, sA, 1/3, eb0, FreeEnd(:A1))
    A2 = PolymerBlock(:A2, sA, 1/3, eb0, FreeEnd(:A2))
    A3 = PolymerBlock(:A3, sA, 1/3, eb0, FreeEnd(:A3))
    return BlockCopolymer(:A3, [A1, A2, A3])
end

function starA4()
    sA = KuhnSegment(:A)
    eb0 = BranchPoint(:EB0)
    A1 = PolymerBlock(:A1, sA, 1/4, eb0, FreeEnd(:A1))
    A2 = PolymerBlock(:A2, sA, 1/4, eb0, FreeEnd(:A2))
    A3 = PolymerBlock(:A3, sA, 1/4, eb0, FreeEnd(:A3))
    A4 = PolymerBlock(:A4, sA, 1/4, eb0, FreeEnd(:A4))
    return BlockCopolymer(:A4, [A1, A2, A3, A4])
end

function starA5()
    sA = KuhnSegment(:A)
    eb0 = BranchPoint(:EB0)
    A1 = PolymerBlock(:A1, sA, 1/5, eb0, FreeEnd(:A1))
    A2 = PolymerBlock(:A2, sA, 1/5, eb0, FreeEnd(:A2))
    A3 = PolymerBlock(:A3, sA, 1/5, eb0, FreeEnd(:A3))
    A4 = PolymerBlock(:A4, sA, 1/5, eb0, FreeEnd(:A4))
    A5 = PolymerBlock(:A5, sA, 1/5, eb0, FreeEnd(:A5))
    return BlockCopolymer(:A5, [A1, A2, A3, A4, A5])
end

function starA6()
    sA = KuhnSegment(:A)
    eb0 = BranchPoint(:EB0)
    A1 = PolymerBlock(:A1, sA, 1/6, eb0, FreeEnd(:A1))
    A2 = PolymerBlock(:A2, sA, 1/6, eb0, FreeEnd(:A2))
    A3 = PolymerBlock(:A3, sA, 1/6, eb0, FreeEnd(:A3))
    A4 = PolymerBlock(:A4, sA, 1/6, eb0, FreeEnd(:A4))
    A5 = PolymerBlock(:A5, sA, 1/6, eb0, FreeEnd(:A5))
    A6 = PolymerBlock(:A6, sA, 1/6, eb0, FreeEnd(:A6))
    return BlockCopolymer(:A6, [A1, A2, A3, A4, A5, A6])
end

function starAB3(; fA=0.4)
    fB = (1 - fA) / 3
    sA = KuhnSegment(:A)
    sB = KuhnSegment(:B)
    eb0 = BranchPoint(:EB0)
    A = PolymerBlock(:A, sA, fA, FreeEnd(:A), eb0)
    B1 = PolymerBlock(:B1, sB, fB, FreeEnd(:B1), eb0)
    B2 = PolymerBlock(:B2, sB, fB, FreeEnd(:B2), eb0)
    B3 = PolymerBlock(:B3, sB, fB, FreeEnd(:B3), eb0)
    return BlockCopolymer(:AB3, [A, B1, B2, B3])
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
    return BlockCopolymer(:AB3A6, [A, B1, B2, B3, A1, A2, A3, A4, A5, A6])
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
    B4 = PolymerBlock(:B4, sB, 0.04, eb[5], eb[1])
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

    return BlockCopolymer(:ABCDO4, [O1, O2, O3, O4, A2, A3, A4, B2, B3, C1, C2, C4, D3, D4, X12, X13, X14])
end

function chainM2(fA, τ)
    fA1 = fA * (τ) / (τ + 14)
    fA2 = fA * (1) / (τ + 14)
    fB = (1 - fA) / 16
    A1 = PolymerBlock(:A1, KuhnSegment(:A), fA1, FreeEnd(:A1), BranchPoint(:EB1))

    A21 = PolymerBlock(:A21, KuhnSegment(:A), fA2, BranchPoint(:EB1), BranchPoint(:EB21))
    A22 = PolymerBlock(:A22, KuhnSegment(:A), fA2, BranchPoint(:EB1), BranchPoint(:EB22))

    A31 = PolymerBlock(:A31, KuhnSegment(:A), fA2, BranchPoint(:EB21), BranchPoint(:EB31))
    A32 = PolymerBlock(:A32, KuhnSegment(:A), fA2, BranchPoint(:EB21), BranchPoint(:EB32))
    A33 = PolymerBlock(:A33, KuhnSegment(:A), fA2, BranchPoint(:EB22), BranchPoint(:EB33))
    A34 = PolymerBlock(:A34, KuhnSegment(:A), fA2, BranchPoint(:EB22), BranchPoint(:EB34))

    A41 = PolymerBlock(:A41, KuhnSegment(:A), fA2, BranchPoint(:EB31), BranchPoint(:EB41))
    A42 = PolymerBlock(:A42, KuhnSegment(:A), fA2, BranchPoint(:EB31), BranchPoint(:EB42))
    A43 = PolymerBlock(:A43, KuhnSegment(:A), fA2, BranchPoint(:EB32), BranchPoint(:EB43))
    A44 = PolymerBlock(:A44, KuhnSegment(:A), fA2, BranchPoint(:EB32), BranchPoint(:EB44))
    A45 = PolymerBlock(:A45, KuhnSegment(:A), fA2, BranchPoint(:EB33), BranchPoint(:EB45))
    A46 = PolymerBlock(:A46, KuhnSegment(:A), fA2, BranchPoint(:EB33), BranchPoint(:EB46))
    A47 = PolymerBlock(:A47, KuhnSegment(:A), fA2, BranchPoint(:EB34), BranchPoint(:EB47))
    A48 = PolymerBlock(:A48, KuhnSegment(:A), fA2, BranchPoint(:EB34), BranchPoint(:EB48))

    B1 = PolymerBlock(:B1, KuhnSegment(:B), fB, BranchPoint(:EB41), BranchPoint(:EB51))
    B2 = PolymerBlock(:B2, KuhnSegment(:B), fB, BranchPoint(:EB41), BranchPoint(:EB52))
    B3 = PolymerBlock(:B3, KuhnSegment(:B), fB, BranchPoint(:EB42), BranchPoint(:EB53))
    B4 = PolymerBlock(:B4, KuhnSegment(:B), fB, BranchPoint(:EB42), BranchPoint(:EB54))
    B5 = PolymerBlock(:B5, KuhnSegment(:B), fB, BranchPoint(:EB43), BranchPoint(:EB55))
    B6 = PolymerBlock(:B6, KuhnSegment(:B), fB, BranchPoint(:EB43), BranchPoint(:EB56))
    B7 = PolymerBlock(:B7, KuhnSegment(:B), fB, BranchPoint(:EB44), BranchPoint(:EB57))
    B8 = PolymerBlock(:B8, KuhnSegment(:B), fB, BranchPoint(:EB44), BranchPoint(:EB58))
    B9 = PolymerBlock(:B9, KuhnSegment(:B), fB, BranchPoint(:EB45), BranchPoint(:EB59))
    B10 = PolymerBlock(:B10, KuhnSegment(:B), fB, BranchPoint(:EB45), BranchPoint(:EB60))
    B11 = PolymerBlock(:B11, KuhnSegment(:B), fB, BranchPoint(:EB46), BranchPoint(:EB61))
    B12 = PolymerBlock(:B12, KuhnSegment(:B), fB, BranchPoint(:EB46), BranchPoint(:EB62))
    B13 = PolymerBlock(:B13, KuhnSegment(:B), fB, BranchPoint(:EB47), BranchPoint(:EB63))
    B14 = PolymerBlock(:B14, KuhnSegment(:B), fB, BranchPoint(:EB47), BranchPoint(:EB64))
    B15 = PolymerBlock(:B15, KuhnSegment(:B), fB, BranchPoint(:EB48), BranchPoint(:EB65))
    B16 = PolymerBlock(:B16, KuhnSegment(:B), fB, BranchPoint(:EB48), BranchPoint(:EB66))

    return BlockCopolymer(:dendron_AB, [A1,A21,A22,A31,A32,A33,A34,A41,A42,A43,A44,A45,A46,A47,A48,
    B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16])
end

function chain_b9(fA::Float64)
    fA_segment = fA / 10
    fB_segment = (1 - fA) / 9

    A1 = PolymerBlock(:A1, KuhnSegment(:A), 2*fA_segment, FreeEnd(:A1), BranchPoint(:EB1))

    A2 = PolymerBlock(:A2, KuhnSegment(:A), fA_segment, BranchPoint(:EB1), BranchPoint(:EB2))
    A3 = PolymerBlock(:A3, KuhnSegment(:A), fA_segment, BranchPoint(:EB2), BranchPoint(:EB3))
    A4 = PolymerBlock(:A4, KuhnSegment(:A), fA_segment, BranchPoint(:EB3), BranchPoint(:EB4))
    A5 = PolymerBlock(:A5, KuhnSegment(:A), fA_segment, BranchPoint(:EB4), BranchPoint(:EB5))
    A6 = PolymerBlock(:A6, KuhnSegment(:A), fA_segment, BranchPoint(:EB5), BranchPoint(:EB6))
    A7 = PolymerBlock(:A7, KuhnSegment(:A), fA_segment, BranchPoint(:EB6), BranchPoint(:EB7))
    A8 = PolymerBlock(:A8, KuhnSegment(:A), fA_segment, BranchPoint(:EB7), BranchPoint(:EB8))
    A9 = PolymerBlock(:A9, KuhnSegment(:A), fA_segment, BranchPoint(:EB8), BranchPoint(:EB9))

    B1 = PolymerBlock(:B1, KuhnSegment(:B), fB_segment, FreeEnd(:B1), BranchPoint(:EB1))
    B2 = PolymerBlock(:B2, KuhnSegment(:B), fB_segment, FreeEnd(:B2), BranchPoint(:EB2))
    B3 = PolymerBlock(:B3, KuhnSegment(:B), fB_segment, FreeEnd(:B3), BranchPoint(:EB3))
    B4 = PolymerBlock(:B4, KuhnSegment(:B), fB_segment, FreeEnd(:B4), BranchPoint(:EB4))
    B5 = PolymerBlock(:B5, KuhnSegment(:B), fB_segment, FreeEnd(:B5), BranchPoint(:EB5))
    B6 = PolymerBlock(:B6, KuhnSegment(:B), fB_segment, FreeEnd(:B6), BranchPoint(:EB6))
    B7 = PolymerBlock(:B7, KuhnSegment(:B), fB_segment, FreeEnd(:B7), BranchPoint(:EB7))
    B8 = PolymerBlock(:B8, KuhnSegment(:B), fB_segment, FreeEnd(:B8), BranchPoint(:EB8))
    B9 = PolymerBlock(:B9, KuhnSegment(:B), fB_segment, FreeEnd(:B9), BranchPoint(:EB9))

    return BlockCopolymer(:comb_AB, [A1, A2, A3, A4, A5, A6, A7, A8, A9,
                                     B1, B2, B3, B4, B5, B6, B7, B8, B9])
end

ab = diblock_chain()
abg = ab |> BlockCopolymerGraph

aba = linearABA()
abag = aba |> BlockCopolymerGraph

branch = branchAB()
branchg = branch |> BlockCopolymerGraph

semibranch = branchAB(0.08, 0.12)
semibranchg = semibranch |> BlockCopolymerGraph

branch2 = branchAB2()
branch2g = branch2 |> BlockCopolymerGraph

ABCACD = chainABCACD()
chainABCACDg = ABCACD |> BlockCopolymerGraph

A6B6 = chainA6B6()
chainA6B6g = A6B6 |> BlockCopolymerGraph

A2B2 = starA2B2()
starA2B2g = A2B2 |> BlockCopolymerGraph

AB3A3 = starAB3A3()
starAB3A3g = AB3A3 |> BlockCopolymerGraph

AB3A6 = starAB3A6()
starAB3A6g = AB3A6 |> BlockCopolymerGraph

AB4A8 = starAB4A8()
starAB4A8g = AB4A8 |> BlockCopolymerGraph

ABCDO = starABCDO()
starABCDOg = ABCDO |> BlockCopolymerGraph

ABCDO4 = starABCDO4()
starABCDO4g = ABCDO4 |> BlockCopolymerGraph

M2 = chainM2(0.25, 10)
M2g = M2 |> BlockCopolymerGraph

chainb9 = chain_b9(0.4)
chainb9g = chainb9 |> BlockCopolymerGraph

# Prevent displaying last line in the REPL
nothing